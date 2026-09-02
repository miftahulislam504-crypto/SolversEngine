"""
Analysis Orchestration — Frontend Model → Solver Call
=========================================================

এই মডিউল frontend থেকে আসা পূর্ণাঙ্গ structural model (Firestore থেকে
পড়া elements, materials, sections, loads — সবগুলো একসাথে) নিয়ে
C++ solver এর জন্য প্রয়োজনীয় node graph বানায়, এবং ফলাফল আবার
frontend-বোধগম্য shape এ ফেরত দেয়।

গুরুত্বপূর্ণ সীমাবদ্ধতা (Phase 4a, আপডেট করা হয়েছে uniform load conversion যোগ হওয়ার পর):
  - Line Element (Beam, Column, Brace, Pile) সম্পূর্ণ সমর্থিত। Slab/
    Wall/Shear-Wall/Core-Wall shell element হিসেবে mesh+solve হয়
    (mesh_generation.py, shell.cpp — membrane+plate bending stiffness),
    কিন্তু shell এর জন্য এখনো internal force/moment recovery নেই (শুধু
    displacement/reaction পাওয়া যায়, stress contour একটা honest
    displacement-magnitude proxy)। Footing/Combined-Footing/Strip-
    Footing/Mat-Foundation/Pile-Cap/Pile-Group — এগুলো foundation-soil
    interaction element, এখনো সমর্থিত না (ভবিষ্যৎ Phase এর কাজ), skip
    হয়ে explicit warning সহ রিপোর্ট হয়, silent-drop না।
  - "pin" connectionType (Brace) এখন প্রয়োগ করা হয় (static condensation,
    cpp/src/stiffness.cpp এর applyEndReleases()), কিন্তু শুধু both-end
    release হিসেবে — single-end release সমর্থিত না (বিস্তারিত নিচে
    element assembly অংশে)।
  - Uniform-line ও uniform-area load এখন equivalent nodal load এ
    রূপান্তরিত হয় (নিচে ধাপ ৪ দেখুন) — uniform-line কে প্রতিটা
    sub-element এ w·L_sub/2 করে দুই প্রান্তে lump করা হয় (simply-
    supported approximation, exact fixed-end-moment consistent load
    না)। uniform-area কে shell mesh এর প্রতিটা quad এ area×intensity÷4
    করে ৪টা corner এ lump করা হয়। দুটোই total force সংরক্ষণ করে
    (নির্ভুল), কিন্তু span/quad এর মাঝামাঝি internal force diagram
    সামান্য approximate — ছোট mesh এ প্রভাব নগণ্য, বড় uniform-area
    load এর ক্ষেত্রে shell এর internal force recovery না থাকার কারণে
    stress/moment ফলাফল এমনিতেও পাওয়া যায় না।
  - Mid-span Point Load এখন সঠিকভাবে হ্যান্ডল হয় (নিচে দেখুন —
    element split করে intermediate node বসানো হয়, "nearest-endpoint
    snap" আচরণ আর নেই)।
"""

import math
from typing import Any

from app.model_conversion import (
    convert_section_to_solver_units,
    convert_material_to_solver_units,
    convert_unit_weight_to_density,
)
from app.mesh_generation import generate_quad_mesh, MeshGenerationError

LINE_ELEMENT_CATEGORIES = {"beam", "column", "brace", "pile"}
# Slab/Wall/Shear-Wall/Core-Wall এখন shell element দিয়ে solve করা যায়
# (mesh_generation.py এর polygon-to-quad-mesh, ও shell.cpp এর 4-node
# quadrilateral shell stiffness — membrane+plate bending)। Footing
# এখনো সমর্থিত না (এটা conceptually একটা point/pad support element,
# slab এর মতো "floating" area element না — foundation-soil interaction
# মডেলিং একটা ভিন্ন সমস্যা, Phase 7 Foundation Module এর কাজ)।
SHELL_ELEMENT_CATEGORIES = {"slab", "wall", "shear-wall", "core-wall"}
UNSUPPORTED_ELEMENT_CATEGORIES = {"footing", "combined-footing", "strip-footing", "mat-foundation", "pile-cap", "pile-group"}

# positionRatio এই মান থেকে endpoint (0 বা 1) এর কতটা কাছে হলে সেটাকে
# "কার্যত endpoint" ধরা হবে (নতুন split node না বসিয়ে সরাসরি existing
# start/end node ব্যবহার করা হবে)। এটা floating-point noise থেকে
# অপ্রয়োজনীয় অতি-ছোট sub-element তৈরি হওয়া ঠেকায়।
ENDPOINT_SNAP_TOLERANCE = 1e-4


class ModelParsingError(Exception):
    """Element/material/section রেফারেন্স broken হলে (যেমন materialId যা library তে নেই)।"""


def _find_by_id(items: list[dict[str, Any]], id_field: str, id_value: str) -> dict[str, Any]:
    for item in items:
        if item.get(id_field) == id_value:
            return item
    raise ModelParsingError(f"{id_field}='{id_value}' খুঁজে পাওয়া যায়নি — library তে অনুপস্থিত")


def _interpolate_point(
    start: dict[str, float], end: dict[str, float], ratio: float
) -> dict[str, float]:
    """start থেকে end পর্যন্ত সরলরেখা বরাবর ratio (0..1) অবস্থানের coordinate।"""
    return {
        "x": start["x"] + (end["x"] - start["x"]) * ratio,
        "y": start["y"] + (end["y"] - start["y"]) * ratio,
        "z": start["z"] + (end["z"] - start["z"]) * ratio,
    }


def _vector_length(start: dict[str, float], end: dict[str, float]) -> float:
    """দুইটা 3D point এর মধ্যে সরলরেখার দৈর্ঘ্য (m, ইনপুট coordinate এর একক অনুযায়ী)।"""
    return math.sqrt(
        (end["x"] - start["x"]) ** 2 + (end["y"] - start["y"]) ** 2 + (end["z"] - start["z"]) ** 2
    )


def _cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _quad_area_3d(corners: list[dict[str, float]]) -> float:
    """
    একটা planar quad (৪টা 3D corner, counter-clockwise) এর ক্ষেত্রফল —
    দুইটা ত্রিভুজে ভেঙে (0-1-2, 0-2-3) প্রতিটার cross-product magnitude/2
    যোগ করে (shoelace formula এর 3D সাধারণীকরণ)। mesh_generation.py এর
    generate_quad_mesh() সবসময় একটা planar polygon থেকে quad তৈরি করে
    (2D projection → unproject), তাই non-planar warping নিয়ে চিন্তা
    করার দরকার নেই।
    """
    p0 = (corners[0]["x"], corners[0]["y"], corners[0]["z"])
    p1 = (corners[1]["x"], corners[1]["y"], corners[1]["z"])
    p2 = (corners[2]["x"], corners[2]["y"], corners[2]["z"])
    p3 = (corners[3]["x"], corners[3]["y"], corners[3]["z"])

    def tri_area(a: tuple[float, float, float], b: tuple[float, float, float], c: tuple[float, float, float]) -> float:
        ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        cx, cy, cz = _cross(ab, ac)
        return 0.5 * math.sqrt(cx ** 2 + cy ** 2 + cz ** 2)

    return tri_area(p0, p1, p2) + tri_area(p0, p2, p3)


class NodeGraph:
    """
    Node তালিকা ও coordinate→index লুকআপ একসাথে বহন করে, যাতে
    build_solver_model এর একাধিক ধাপ (element split point বসানো, load
    application point বসানো) একই get_or_create_node() যুক্তি পুনরায়
    ব্যবহার করতে পারে কোনো ডুপ্লিকেট node তৈরি না করে।
    """

    def __init__(self) -> None:
        self.node_list: list[dict[str, float]] = []
        self._coordinate_to_index: dict[tuple[float, float, float], int] = {}

    def get_or_create_node(self, x: float, y: float, z: float) -> int:
        # ১ মিলিমিটার precision এ রাউন্ড করে key বানানো হচ্ছে, যাতে
        # floating-point এর সামান্য পার্থক্য (যেমন 5.000000001 বনাম
        # 4.999999998) ভুলবশত দুটো আলাদা node হিসেবে গণ্য না হয়।
        key = (round(x, 3), round(y, 3), round(z, 3))
        if key in self._coordinate_to_index:
            return self._coordinate_to_index[key]

        index = len(self.node_list)
        self.node_list.append({"nodeId": f"node-{index}", "x": x, "y": y, "z": z})
        self._coordinate_to_index[key] = index
        return index

    def index_of(self, x: float, y: float, z: float) -> int | None:
        key = (round(x, 3), round(y, 3), round(z, 3))
        return self._coordinate_to_index.get(key)


def build_node_graph(
    elements: list[dict[str, Any]],
) -> tuple[list[dict[str, float]], dict[tuple[float, float, float], int]]:
    """
    Line element গুলোর start/end point থেকে একটা unique node তালিকা
    বানায়। দুটো element যদি একই coordinate এ মিলিত হয় (যেমন Beam এর
    end point = Column এর top point), তাদের একই node হিসেবে merge
    করা হয়।

    নোট: এই ফাংশন শুধু element endpoint (start/end) থেকে node বানায় —
    mid-span load split point এখানে অন্তর্ভুক্ত না (সেটা
    build_solver_model এ আলাদাভাবে হ্যান্ডল হয়, কারণ split point শুধু
    load-application-driven, geometry-driven না)। এই ফাংশন backward-
    compatible রাখা হয়েছে (পুরনো signature/return shape অপরিবর্তিত)।
    """
    graph = NodeGraph()
    for element in elements:
        if element["category"] not in LINE_ELEMENT_CATEGORIES:
            continue
        start = element["startPoint"]
        end = element["endPoint"]
        graph.get_or_create_node(start["x"], start["y"], start["z"])
        graph.get_or_create_node(end["x"], end["y"], end["z"])

    return graph.node_list, graph._coordinate_to_index


def build_solver_model(
    elements: list[dict[str, Any]],
    materials: list[dict[str, Any]],
    sections: list[dict[str, Any]],
    load_cases: list[dict[str, Any]],
    supported_node_ids: set[str],
    support_overrides: list[dict[str, Any]] | None = None,
) -> tuple[dict[str, Any], list[str], list[dict[str, Any]]]:
    """
    পূর্ণাঙ্গ solver input dict বানায়, একটা warnings তালিকা, এবং একটা
    sub_element_registry রিটার্ন করে।

    sub_element_registry: mid-span load এর কারণে element split হলে
    কোন solver_elements[i] কোন original elementId থেকে এসেছে তার
    ম্যাপিং — [{"originalElementId": str, "solverElementIndex": int,
    "subStartRatio": float, "subEndRatio": float}, ...]। Split না
    হলে প্রতিটা original element এর জন্য একটাই entry থাকে (ratio
    0.0→1.0)। main.py এই তালিকা ব্যবহার করে elementEndForces কে আবার
    original elementId এর সাথে যুক্ত করে ফেরত দেয় (C++ solver নিজে
    শুধু positional index জানে, elementId-aware না)।

    supported_node_ids: যেসব node এ boundary condition বসাতে হবে
    (Phase 4a তে একটা সরলীকৃত heuristic ব্যবহার করা হচ্ছে: base story
    এ থাকা সব node কে automatically fully-fixed support ধরা হয়)।

    support_overrides (Phase 5): ঐচ্ছিক — coordinate-ভিত্তিক support
    override, Y≈0 heuristic কে override/supplement করে। প্রতিটা entry:
    {"x", "y", "z", "supportType": "fixed"|"pinned"|"free"|"custom",
    "restrainX"..."restrainRz" (শুধু supportType="custom" হলে পড়া হয়)}।
    "fixed" মানে সব ৬ DOF restrain, "pinned" মানে শুধু ৩টা translation
    restrain (rotation free — geotechnically একটা সরলীকরণ, প্রকৃত pile/
    footing rotational stiffness এখানে model করা হয় না), "free" মানে
    কোনো DOF-ই restrain না (ওই coordinate এ heuristic থেকে support বাদ
    দিতে, যেমন ভুলভাবে Y≈0 এ ধরা পড়া non-support node)। কোনো override
    coordinate কোনো node এর সাথে না মিললে warning এ জানানো হয়, silent
    ignore হয় না। override প্রয়োগের ক্রম: প্রথমে Y≈0 heuristic দিয়ে base
    boundary_conditions dict বানানো হয় (নিচে), তারপর override গুলো সেই
    dict এ coordinate-key দিয়ে বসানো/প্রতিস্থাপন করা হয় — অর্থাৎ override
    সবসময় heuristic কে জেতে, override না থাকা node গুলো heuristic
    অনুযায়ীই থাকে।
    """
    warnings: list[str] = []

    # ইনপুট validation: প্রতিটা element-এ 'category' ও 'elementId' key
    # থাকতেই হবে — নিচের পুরো ফাংশন জুড়ে এই দুটো key ধরে নিয়েই কাজ করা
    # হয়। এই আপফ্রন্ট চেক না থাকলে malformed/incomplete element (যেমন
    # কোনো diagnostic/placeholder payload, বা future client bug) সরাসরি
    # raw KeyError তুলে পুরো সার্ভিসকে unhandled 500 crash করাত। এখন এর
    # বদলে একটা readable ModelParsingError ওঠে, যা main.py 422 এ ম্যাপ
    # করে দেয়।
    for idx, e in enumerate(elements):
        missing = [k for k in ("elementId", "category") if k not in e]
        if missing:
            raise ModelParsingError(
                f"elements[{idx}] এ আবশ্যক key অনুপস্থিত: {', '.join(missing)}। "
                f"প্রতিটা element এ কমপক্ষে 'elementId' ও 'category' থাকতে হবে।"
            )

    skipped_elements = [e for e in elements if e["category"] in UNSUPPORTED_ELEMENT_CATEGORIES]
    if skipped_elements:
        categories = sorted({e["category"] for e in skipped_elements})
        warnings.append(
            f"⚠️ {len(skipped_elements)}টা element ({', '.join(categories)}) এই Phase এ solve করা "
            f"যায়নি — এগুলো foundation-soil interaction element, যা এখনো implement করা হয়নি "
            f"(ভবিষ্যৎ Phase এর কাজ)। ফলাফলে এই element গুলোর প্রভাব ধরা পড়েনি।"
        )

    line_elements = [e for e in elements if e["category"] in LINE_ELEMENT_CATEGORIES]
    shell_source_elements = [e for e in elements if e["category"] in SHELL_ELEMENT_CATEGORIES]

    if not line_elements and not shell_source_elements:
        raise ModelParsingError(
            "কোনো solve-যোগ্য (Beam/Column/Brace/Pile/Slab/Wall/Shear-Wall/Core-Wall) element "
            "পাওয়া যায়নি"
        )

    brace_elements = [e for e in line_elements if e["category"] == "brace"]
    if brace_elements:
        warnings.append(
            f"ℹ️ {len(brace_elements)}টা Brace element আছে — pin-connection release প্রয়োগ করা "
            f"হয়েছে (উভয় প্রান্তে bending moment DOF release, static condensation দিয়ে)। "
            f"⚠️ সীমাবদ্ধতা: এই মুহূর্তে single-end release সমর্থিত না — connectionType='pin' "
            f"হলে element-এর উভয় প্রান্তই release হয়, একদিকে moment/অন্যদিকে pin এমন আংশিক "
            f"release মডেল করা যায় না।"
        )

    # ধাপ ১: প্রতিটা line element-এ কোন positionRatio-তে mid-span point
    # load আছে তা আগে থেকে জোগাড় করা, যাতে element split (ধাপ ৩) একবারেই
    # সব split point জেনে করা যায়।
    split_ratios_by_element: dict[str, set[float]] = {}
    point_loads_only = [lc for lc in load_cases if lc["applicationType"] == "point"]
    shell_source_ids = {e["elementId"] for e in shell_source_elements}
    unsupported_shell_point_loads = 0
    for load_case in point_loads_only:
        if load_case["elementId"] in shell_source_ids:
            # Shell element এ Point Load এখনো সমর্থিত না (এই ধাপে শুধু
            # geometry/mesh solve করা হচ্ছে, shell এ load application
            # ভবিষ্যতের কাজ) — silently ignore না করে গণনা রাখা হচ্ছে,
            # নিচে একটা warning এ জানানো হবে।
            unsupported_shell_point_loads += 1
            continue
        target_element = _find_by_id(line_elements, "elementId", load_case["elementId"])
        ratio = load_case.get("positionRatio", 0.5)
        if ENDPOINT_SNAP_TOLERANCE < ratio < (1 - ENDPOINT_SNAP_TOLERANCE):
            split_ratios_by_element.setdefault(target_element["elementId"], set()).add(ratio)

    if unsupported_shell_point_loads:
        warnings.append(
            f"⚠️ {unsupported_shell_point_loads}টা Point Load একটা Slab/Wall/Shear-Wall/Core-Wall "
            f"element কে target করছে — shell element এ load application এখনো সমর্থিত না এই ধাপে "
            f"(শুধু geometry/stiffness solve করা হচ্ছে), তাই এই load গুলো বাদ দেওয়া হয়েছে।"
        )

    # ধাপ ১.৫: প্রতিটা Slab/Wall/Shear-Wall/Core-Wall element কে quad
    # mesh এ রূপান্তর করা (mesh_generation.py) — node graph বানানোর
    # আগে করা হচ্ছে, যাতে mesh-generated quad corner point গুলোও
    # (যার মধ্যে নতুন centroid/midpoint node থাকবে যা কোনো original
    # element এর endpoint না) node graph এ অন্তর্ভুক্ত হয়।
    #
    # shell_meshes: elementId → list of quad (প্রতিটা quad = ৪টা
    # {"x":..,"y":..,"z":..} dict এর তালিকা, counter-clockwise)
    shell_meshes: dict[str, list[list[dict[str, float]]]] = {}
    for shell_source in shell_source_elements:
        try:
            quads_3d = generate_quad_mesh(shell_source["vertices"])
        except MeshGenerationError as mesh_error:
            raise ModelParsingError(
                f"'{shell_source.get('name', shell_source['elementId'])}' ({shell_source['category']}) "
                f"এর mesh generation ব্যর্থ হয়েছে: {mesh_error}"
            )
        shell_meshes[shell_source["elementId"]] = [
            [{"x": p[0], "y": p[1], "z": p[2]} for p in quad] for quad in quads_3d
        ]

    if shell_source_elements:
        total_quads = sum(len(qs) for qs in shell_meshes.values())
        warnings.append(
            f"ℹ️ {len(shell_source_elements)}টা Slab/Wall/Shear-Wall/Core-Wall element মোট "
            f"{total_quads}টা shell (quad) element এ mesh করা হয়েছে (ear-clipping triangulation + "
            f"triangle-to-quad conversion)। ⚠️ সীমাবদ্ধতা: mesh density uniform/adaptive না — "
            f"বড় বা জটিল-আকৃতির element এ mesh coarse হতে পারে, যা displacement/stress নির্ভুলতা "
            f"প্রভাবিত করতে পারে। Element internal force (moment/shear per unit width) এই ধাপে "
            f"এখনো রিপোর্ট করা হয় না — শুধু displacement ফলাফল পাওয়া যাবে।"
        )

    # ধাপ ২: Node graph — element endpoint (যথারীতি) + প্রতিটা split-
    # প্রয়োজনীয় element-এর intermediate node(গুলো) + shell mesh corner।
    graph = NodeGraph()
    for element in line_elements:
        start = element["startPoint"]
        end = element["endPoint"]
        graph.get_or_create_node(start["x"], start["y"], start["z"])
        graph.get_or_create_node(end["x"], end["y"], end["z"])

    for element in line_elements:
        ratios = split_ratios_by_element.get(element["elementId"])
        if not ratios:
            continue
        start = element["startPoint"]
        end = element["endPoint"]
        for ratio in ratios:
            p = _interpolate_point(start, end, ratio)
            graph.get_or_create_node(p["x"], p["y"], p["z"])

    for quads in shell_meshes.values():
        for quad in quads:
            for corner in quad:
                graph.get_or_create_node(corner["x"], corner["y"], corner["z"])

    node_list = graph.node_list

    def point_index(point: dict[str, float]) -> int:
        idx = graph.index_of(point["x"], point["y"], point["z"])
        if idx is None:
            raise ModelParsingError(
                f"Node coordinate ({point['x']}, {point['y']}, {point['z']}) node graph এ পাওয়া "
                f"যায়নি — এটা একটা internal bug।"
            )
        return idx

    # ধাপ ৩: প্রতিটা line element কে এক বা একাধিক sub-element এ রূপান্তর।
    # Split প্রয়োজন না হলে element অপরিবর্তিত (১টা sub-element)। Split
    # হলে sorted ratio অনুযায়ী ধারাবাহিক sub-element — প্রতিটার section/
    # material মূল element থেকে অভিন্ন।
    #
    # connectionType handling: internal cut point এ সবসময় rigid/moment
    # continuity — mid-span লোড শুধু numerical split, physical hinge না।
    sub_element_registry: list[dict[str, Any]] = []
    solver_elements: list[dict[str, Any]] = []

    for element in line_elements:
        material = _find_by_id(materials, "materialId", element["materialId"])
        section = _find_by_id(sections, "sectionId", element["sectionId"])

        if "properties" not in section:
            raise ModelParsingError(
                f"Section '{section.get('name', element['sectionId'])}' এর geometric properties "
                f"পাওয়া যায়নি — request payload এ precomputed section properties অন্তর্ভুক্ত করা "
                f"আছে কিনা যাচাই করুন (frontend এর computeSectionProperties() থেকে)।"
            )

        # yieldMomentMzKNm/yieldMomentMyKNm (Nonlinear Static Analysis)
        # frontend এর BaseSection টাইপে top-level field হিসেবে থাকে
        # (src/lib/types/section.ts), কিন্তু section["properties"] হলো
        # computeSectionProperties() এর আউটপুট (শুধু geometry:
        # area/ixx/iyy/j/centroidY) — সেখানে yield capacity থাকে না।
        # তাই এখানে দুটো উৎস merge করা হচ্ছে convert_section_to_solver_units()
        # এ পাঠানোর আগে।
        section_properties_with_yield = {
            **section["properties"],
            "yieldMomentMzKNm": section.get("yieldMomentMzKNm", 0.0),
            "yieldMomentMyKNm": section.get("yieldMomentMyKNm", 0.0),
        }
        section_m = convert_section_to_solver_units(section_properties_with_yield)

        if material["type"] == "concrete":
            fc = material["fc"]
            ec_override = material.get("ec")
            elastic_modulus_mpa = ec_override if ec_override is not None else 4700 * (fc ** 0.5)
            poissons_ratio = material.get("poissonsRatio", 0.2)
        elif material["type"] == "steel":
            elastic_modulus_mpa = material["es"]
            poissons_ratio = material.get("poissonsRatio", 0.3)
        else:
            raise ModelParsingError(
                f"Material type '{material['type']}' এখনো solver এ সমর্থিত না (শুধু concrete/steel)।"
            )

        shear_modulus_mpa = elastic_modulus_mpa / (2 * (1 + poissons_ratio))
        material_converted = convert_material_to_solver_units(elastic_modulus_mpa, shear_modulus_mpa)
        # Modal Analysis এর mass matrix গণনার জন্য material density যোগ
        # করা হচ্ছে — Linear Static Analysis এ এই key ব্যবহৃত হয় না
        # (bindings.cpp এ optional হিসেবে handle করা আছে), তাই এই যোগ
        # করা Linear Static এর কোনো existing আচরণ বদলায় না।
        material_converted["density"] = convert_unit_weight_to_density(material["unitWeight"])

        connection_type = element.get("connectionType", "moment")
        ratios = sorted(split_ratios_by_element.get(element["elementId"], set()))
        boundary_ratios = [0.0] + ratios + [1.0]

        start_pt = element["startPoint"]
        end_pt = element["endPoint"]

        for i in range(len(boundary_ratios) - 1):
            r0 = boundary_ratios[i]
            r1 = boundary_ratios[i + 1]
            p0 = start_pt if r0 == 0.0 else _interpolate_point(start_pt, end_pt, r0)
            p1 = end_pt if r1 == 1.0 else _interpolate_point(start_pt, end_pt, r1)

            is_first_sub = (i == 0)
            is_last_sub = (i == len(boundary_ratios) - 2)
            sub_connection_type = connection_type if (is_first_sub or is_last_sub) else "moment"
            # নোট: বর্তমান solver "pin" কে সবসময় both-end release হিসেবে
            # treat করে (getEffectiveLocalStiffness, stiffness.h দেখুন)।
            # তাই একটা multi-sub-element pin member এ sub_connection_type=
            # "pin" থাকলে সেই sub-element এর *দুই* প্রান্তই release হবে,
            # শুধু বাইরের প্রান্ত না — এটা split pin element এর জন্য
            # conservative approximation, নিচে warning এ জানানো আছে।

            # Plastic hinge assignment (Nonlinear Static Analysis) — শুধু
            # physical outer প্রান্তে বসে (connectionType এর মতোই
            # যুক্তিতে), internal split (mid-span load cut) point এ না।
            # একটা physical hinge মূল element এর প্রকৃত প্রান্তে থাকে
            # (যেমন column base) — mid-span split সম্পূর্ণ numerical
            # কারণে, সেই কৃত্রিম cut point এ hinge বসানো ভুল হবে (একটা
            # continuous element কে ভুলভাবে দুইটা যৌগিক hinge-behavior
            # সহ element বানিয়ে ফেলবে)।
            sub_hinge_at_start = element.get("hingeAtStart", False) if is_first_sub else False
            sub_hinge_at_end = element.get("hingeAtEnd", False) if is_last_sub else False

            solver_element_index = len(solver_elements)
            solver_elements.append({
                "elementId": f"{element['elementId']}#{i}" if len(boundary_ratios) > 2 else element["elementId"],
                "startNodeIndex": point_index(p0),
                "endNodeIndex": point_index(p1),
                "connectionType": sub_connection_type,
                "section": section_m,
                "material": material_converted,
                "hingeAtStart": sub_hinge_at_start,
                "hingeAtEnd": sub_hinge_at_end,
            })
            sub_element_registry.append({
                "originalElementId": element["elementId"],
                "solverElementIndex": solver_element_index,
                "subStartRatio": r0,
                "subEndRatio": r1,
            })

    split_pin_elements = [
        e["elementId"] for e in line_elements
        if e.get("connectionType") == "pin" and e["elementId"] in split_ratios_by_element
    ]
    if split_pin_elements:
        warnings.append(
            f"⚠️ {len(split_pin_elements)}টা pin-connected element ({', '.join(split_pin_elements)}) "
            f"এ mid-span point load থাকায় সেগুলো split করা হয়েছে — কিন্তু বর্তমান সলভার 'pin' কে "
            f"সবসময় sub-element এর উভয় প্রান্তে release করে, তাই split point (internal cut) এও "
            f"অনিচ্ছাকৃতভাবে moment release হয়ে যাচ্ছে। ছোট মডেলে প্রভাব সীমিত হতে পারে, কিন্তু "
            f"নির্ভুল ফলাফলের জন্য এই সীমাবদ্ধতা মাথায় রাখুন।"
        )

    # -----------------------------------------------------------------
    # ধাপ ৩.৫: Shell (Slab/Wall/Shear-Wall/Core-Wall) element assembly —
    # প্রতিটা mesh-generated quad কে একটা solver shellElement এ রূপান্তর
    # করা, node graph এ ইতিমধ্যে যোগ হওয়া corner coordinate ব্যবহার করে
    # (point_index() দিয়ে lookup)।
    # -----------------------------------------------------------------
    solver_shell_elements: list[dict[str, Any]] = []
    for shell_source in shell_source_elements:
        material = _find_by_id(materials, "materialId", shell_source["materialId"])

        if material["type"] == "concrete":
            fc = material["fc"]
            ec_override = material.get("ec")
            elastic_modulus_mpa = ec_override if ec_override is not None else 4700 * (fc ** 0.5)
            poissons_ratio = material.get("poissonsRatio", 0.2)
        elif material["type"] == "steel":
            elastic_modulus_mpa = material["es"]
            poissons_ratio = material.get("poissonsRatio", 0.3)
        else:
            raise ModelParsingError(
                f"Material type '{material['type']}' এখনো shell element এ সমর্থিত না (শুধু "
                f"concrete/steel)।"
            )

        shell_material_converted = convert_material_to_solver_units(elastic_modulus_mpa, elastic_modulus_mpa)
        # নোট: convert_material_to_solver_units() দুইটা আর্গুমেন্ট নেয়
        # (E ও G), কিন্তু shell element G সরাসরি ব্যবহার করে না (G,
        # poissonsRatio থেকে shell.cpp এর ভেতরেই derive হয় shear
        # constitutive matrix এর জন্য) — তাই এখানে শুধু elasticModulus
        # ব্যবহার করা হবে dict থেকে, shearModulus entry থাকলেও ignore
        # হবে (bindings.cpp এর parseModelFromDict এ shell এর জন্য
        # shearModulus ঐচ্ছিক)।
        shell_material_converted["poissonsRatio"] = poissons_ratio
        shell_material_converted["density"] = convert_unit_weight_to_density(material["unitWeight"])

        thickness_m = shell_source["thickness"] / 1000.0  # mm → m (element.ts এ thickness: number // mm)

        for quad_idx, quad in enumerate(shell_meshes[shell_source["elementId"]]):
            node_indices = [point_index(corner) for corner in quad]
            solver_shell_elements.append({
                "elementId": f"{shell_source['elementId']}#{quad_idx}",
                "nodeIndices": node_indices,
                "thickness": thickness_m,
                "material": shell_material_converted,
            })

    # Boundary conditions — base-level heuristic
    boundary_conditions_by_node_index: dict[int, dict[str, Any]] = {}
    for i, node in enumerate(node_list):
        if node["nodeId"] in supported_node_ids or node["y"] <= 1e-3:
            boundary_conditions_by_node_index[i] = {
                "nodeIndex": i,
                "restrainX": True, "restrainY": True, "restrainZ": True,
                "restrainRx": True, "restrainRy": True, "restrainRz": True,
            }

    heuristic_count = len(boundary_conditions_by_node_index)
    if heuristic_count:
        warnings.append(
            f"ℹ️ {heuristic_count}টা node কে স্বয়ংক্রিয়ভাবে fully-fixed support ধরা "
            f"হয়েছে (elevation Y≈0 হওয়ার ভিত্তিতে) — এটা একটা preliminary heuristic, প্রকৃত "
            f"support condition define করার UI এখনো নেই। এই মুহূর্তে ফলাফল সেই সতর্কতার সাথে "
            f"ব্যবহার করুন।"
        )
    elif not support_overrides:
        warnings.append(
            "⚠️ কোনো node base level (Y≈0) এ পাওয়া যায়নি — কোনো boundary condition বসানো "
            "হয়নি, সলভার ব্যর্থ হবে (unstable structure)।"
        )

    # Phase 5 — support_overrides প্রয়োগ। heuristic এর ফলাফল উপরে
    # ইতিমধ্যে তৈরি হয়ে গেছে, override সেটাকে coordinate-key দিয়ে
    # বসিয়ে/প্রতিস্থাপন করে — override সবসময় heuristic কে জেতে।
    SUPPORT_TYPE_DOF_PRESETS = {
        "fixed": {"restrainX": True, "restrainY": True, "restrainZ": True, "restrainRx": True, "restrainRy": True, "restrainRz": True},
        "pinned": {"restrainX": True, "restrainY": True, "restrainZ": True, "restrainRx": False, "restrainRy": False, "restrainRz": False},
        "free": {"restrainX": False, "restrainY": False, "restrainZ": False, "restrainRx": False, "restrainRy": False, "restrainRz": False},
    }
    if support_overrides:
        applied_count = 0
        unmatched_overrides = 0
        for override in support_overrides:
            node_index = graph.index_of(override["x"], override["y"], override["z"])
            if node_index is None:
                unmatched_overrides += 1
                continue

            support_type = override.get("supportType", "fixed")
            if support_type == "custom":
                dof_flags = {
                    "restrainX": bool(override.get("restrainX", False)),
                    "restrainY": bool(override.get("restrainY", False)),
                    "restrainZ": bool(override.get("restrainZ", False)),
                    "restrainRx": bool(override.get("restrainRx", False)),
                    "restrainRy": bool(override.get("restrainRy", False)),
                    "restrainRz": bool(override.get("restrainRz", False)),
                }
            elif support_type in SUPPORT_TYPE_DOF_PRESETS:
                dof_flags = SUPPORT_TYPE_DOF_PRESETS[support_type]
            else:
                warnings.append(
                    f"⚠️ support_overrides এ অজানা supportType (\"{support_type}\") — এই override "
                    f"({override['x']}, {override['y']}, {override['z']}) বাদ দেওয়া হয়েছে।"
                )
                continue

            if support_type == "free":
                boundary_conditions_by_node_index.pop(node_index, None)
            else:
                boundary_conditions_by_node_index[node_index] = {"nodeIndex": node_index, **dof_flags}
            applied_count += 1

        if applied_count:
            warnings.append(
                f"ℹ️ {applied_count}টা support_override প্রয়োগ করা হয়েছে (Y≈0 heuristic কে "
                f"override/supplement করে) — উপরের fully-fixed heuristic count এ এই পরিবর্তন "
                f"প্রতিফলিত না, নিচের চূড়ান্ত boundary condition তালিকাতেই আসল ফলাফল।"
            )
        if unmatched_overrides:
            warnings.append(
                f"⚠️ {unmatched_overrides}টা support_override এর coordinate মডেলের কোনো node এর "
                f"সাথে মেলেনি (নীরবে বাদ দেওয়া হয়েছে) — coordinate ভুল, বা সেই বিন্দুতে কোনো element "
                f"endpoint নেই কিনা যাচাই করুন।"
            )
        if not boundary_conditions_by_node_index:
            warnings.append(
                "⚠️ সব boundary condition override দিয়ে সরিয়ে ফেলা হয়েছে (supportType='free') — "
                "কোনো boundary condition বসানো হয়নি, সলভার ব্যর্থ হবে (unstable structure)।"
            )

    boundary_conditions = list(boundary_conditions_by_node_index.values())

    # ধাপ ৪: Load conversion — point load সরাসরি actual position-এর node-এ
    # apply হয় (mid-span split node, নিচে)। uniform-line/uniform-area load
    # কে "equivalent nodal load" এ রূপান্তর করা হয় (নিচে বিস্তারিত) — এটা
    # standard FE practice, কোনো নতুন element type বা solver পরিবর্তন
    # লাগে না, C++ সলভার শুধু NodalLoad-ই চেনে (assembleGlobalLoadVector()
    # প্রতিটা node-এ একাধিক load entry থাকলে += দিয়ে সরাসরি sum করে, তাই
    # একই node-এ একাধিক sub-element/quad থেকে contribution আলাদা entry
    # হিসেবে append করাই যথেষ্ট, আগে থেকে accumulate করার দরকার নেই)।
    solver_loads = []
    unsupported_load_types = set()
    uniform_load_applied = False

    # elementId → তার সব sub-element registry entry (একাধিক হতে পারে যদি
    # mid-span point load এর কারণে split হয়ে থাকে) — uniform-line load কে
    # প্রতিটা sub-segment এ আলাদাভাবে distribute করার জন্য।
    sub_elements_by_original: dict[str, list[dict[str, Any]]] = {}
    for entry in sub_element_registry:
        sub_elements_by_original.setdefault(entry["originalElementId"], []).append(entry)

    for load_case in load_cases:
        application_type = load_case["applicationType"]

        if application_type == "point":
            if load_case["elementId"] in shell_source_ids:
                continue  # উপরে (split-ratio scanning ধাপে) ইতিমধ্যে warning যোগ করা হয়েছে

            target_element = _find_by_id(line_elements, "elementId", load_case["elementId"])
            ratio = load_case.get("positionRatio", 0.5)

            if ratio <= ENDPOINT_SNAP_TOLERANCE:
                point = target_element["startPoint"]
            elif ratio >= (1 - ENDPOINT_SNAP_TOLERANCE):
                point = target_element["endPoint"]
            else:
                point = _interpolate_point(target_element["startPoint"], target_element["endPoint"], ratio)

            node_idx = point_index(point)

            solver_loads.append({
                "nodeIndex": node_idx,
                "fx": load_case["forceX"], "fy": load_case["forceY"], "fz": load_case["forceZ"],
                "mx": 0, "my": 0, "mz": 0,
            })

        elif application_type == "uniform-line":
            if load_case["elementId"] in shell_source_ids:
                continue

            target_element = _find_by_id(line_elements, "elementId", load_case["elementId"])
            total_length = _vector_length(target_element["startPoint"], target_element["endPoint"])
            sub_entries = sub_elements_by_original.get(load_case["elementId"], [])

            if not sub_entries:
                # mid-span point load না থাকলে element split হয়নি, কিন্তু
                # sub_element_registry তে তবুও একটা single "0.0→1.0"
                # sub-element থাকার কথা (build_solver_model এর element
                # assembly loop সবসময় অন্তত একটা sub-element রেজিস্টার
                # করে) — তাই খালি থাকলে সেটা একটা প্রকৃত অসঙ্গতি, নীরবে
                # বাদ না দিয়ে warning দেওয়া হচ্ছে।
                warnings.append(
                    f"⚠️ Uniform Line Load এর element '{load_case['elementId']}' এর জন্য কোনো "
                    f"sub-element registry entry পাওয়া যায়নি — এই লোড প্রয়োগ করা যায়নি।"
                )
                continue

            # frontend এর intensityX/Y/Z পুরো original element এর length
            # বরাবর ধ্রুবক (kN/m, module-data.ts এর UniformLineLoadCase
            # দেখুন) — প্রতিটা sub-segment এ তার নিজস্ব দৈর্ঘ্য অনুপাতে
            # total force ভাগ হবে, তারপর সেই sub-segment এর দুই প্রান্তে
            # w·L_sub/2 করে lump হবে (simply-supported equivalent nodal
            # load — standard consistent/lumped load conversion, exact
            # fixed-end-moment consistent load না, কিন্তু ছোট mesh এ
            # যথেষ্ট নির্ভুল approximation, warning এ জানানো আছে)।
            intensity = (
                load_case.get("forceX", 0.0),
                load_case.get("forceY", 0.0),
                load_case.get("forceZ", 0.0),
            )
            for sub in sub_entries:
                sub_length = (sub["subEndRatio"] - sub["subStartRatio"]) * total_length
                solver_element = solver_elements[sub["solverElementIndex"]]
                half_force = tuple(v * sub_length / 2.0 for v in intensity)

                solver_loads.append({
                    "nodeIndex": solver_element["startNodeIndex"],
                    "fx": half_force[0], "fy": half_force[1], "fz": half_force[2],
                    "mx": 0, "my": 0, "mz": 0,
                })
                solver_loads.append({
                    "nodeIndex": solver_element["endNodeIndex"],
                    "fx": half_force[0], "fy": half_force[1], "fz": half_force[2],
                    "mx": 0, "my": 0, "mz": 0,
                })
            uniform_load_applied = True

        elif application_type == "uniform-area":
            if load_case["elementId"] not in shell_source_ids:
                warnings.append(
                    f"⚠️ Uniform Area Load এর element '{load_case['elementId']}' একটা shell "
                    f"(Slab/Wall) না — বাদ দেওয়া হয়েছে।"
                )
                continue

            quads = shell_meshes.get(load_case["elementId"], [])
            if not quads:
                warnings.append(
                    f"⚠️ Uniform Area Load এর element '{load_case['elementId']}' এর কোনো mesh quad "
                    f"পাওয়া যায়নি (mesh generation ব্যর্থ হয়ে থাকতে পারে) — এই লোড প্রয়োগ করা যায়নি।"
                )
                continue

            # intensity সবসময় gravity (global -Y) দিক ধরা হয় (frontend এর
            # UniformAreaLoadCase এ শুধু "intensity" (kN/m²) থাকে, কোনো
            # axis component না — module-data.ts এর কমেন্ট "সাধারণত gravity
            # load" নিশ্চিত করে)। প্রতিটা quad এর area × intensity ÷ 4 করে
            # তার ৪টা corner node এ সমানভাবে lump করা হয় (consistent load
            # না, কিন্তু uniform pressure এর জন্য uniform mesh এ ভালো
            # approximation)।
            intensity_kn_m2 = load_case.get("intensity", 0.0)
            for quad in quads:
                quad_area = _quad_area_3d(quad)
                quad_total_force = intensity_kn_m2 * quad_area
                force_per_corner = quad_total_force / 4.0

                for corner in quad:
                    node_idx = point_index(corner)
                    solver_loads.append({
                        "nodeIndex": node_idx,
                        "fx": 0.0, "fy": force_per_corner, "fz": 0.0,
                        "mx": 0, "my": 0, "mz": 0,
                    })
            uniform_load_applied = True

        else:
            unsupported_load_types.add(application_type)

    if unsupported_load_types:
        warnings.append(
            f"⚠️ এই লোড টাইপ গুলো এই Phase এ সমর্থিত না, বাদ দেওয়া হয়েছে: "
            f"{', '.join(sorted(unsupported_load_types))}।"
        )

    if uniform_load_applied:
        warnings.append(
            "ℹ️ Uniform Line/Area Load কে equivalent nodal load এ রূপান্তর করা হয়েছে "
            "(simply-supported lumped approximation, exact fixed-end-moment consistent load না) — "
            "তাই এই লোড থেকে আসা reaction/displacement মোটামুটি নির্ভুল কিন্তু sub-element "
            "internal force diagram এ (span এর মাঝামাঝি) সামান্য পার্থক্য থাকতে পারে প্রকৃত "
            "distributed-load internal force curve এর তুলনায়। Uniform Area Load এর ক্ষেত্রে "
            "shell element এ এখনো internal force/moment recovery নেই (শুধু displacement solve "
            "হয়), তাই সেই লোড থেকে stress/moment ফলাফল পাওয়া যাবে না।"
        )

    if not solver_loads and not solver_shell_elements:
        warnings.append("⚠️ কোনো ব্যবহারযোগ্য load পাওয়া যায়নি — সলভার শূন্য-লোড অবস্থায় চলবে।")
    elif not solver_loads and solver_shell_elements:
        warnings.append("ℹ️ কোনো load প্রয়োগযোগ্য পাওয়া যায়নি — শুধু shell (Slab/Wall) geometry/stiffness solve হয়েছে।")

    solver_input = {
        "nodes": node_list,
        "elements": solver_elements,
        "shellElements": solver_shell_elements,
        "boundaryConditions": boundary_conditions,
        "loads": solver_loads,
    }

    return solver_input, warnings, sub_element_registry
