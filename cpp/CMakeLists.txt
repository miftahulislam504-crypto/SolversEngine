cmake_minimum_required(VERSION 3.15)
project(civilos_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# pybind11 কে pip দিয়ে ইনস্টল করা হয় (requirements.txt এ), CMake config
# pybind11 এর নিজস্ব প্যাকেজের ভিতরেই আছে — pip install pybind11 করলে
# `pybind11 --includes`/cmake dir স্বয়ংক্রিয়ভাবে খুঁজে পাওয়া যায়।
execute_process(
    COMMAND python3 -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_package(pybind11 REQUIRED)

# --- Core solver library (stiffness + assembly + solve) ---
# Object library হিসেবে রাখা হয়েছে যাতে Python module এবং test
# executable দুটোই একই compiled object পুনর্ব্যবহার করতে পারে,
# একই সোর্স দুইবার আলাদাভাবে কম্পাইল না করে।
#
# POSITION_INDEPENDENT_CODE অত্যাবশ্যক এখানে: pybind11_add_module
# একটা shared library (.so) বানায় Python module এর জন্য, এবং shared
# library তে link হওয়া যেকোনো object code position-independent
# (-fPIC) হতে হবে — object library ডিফল্টভাবে এটা inherit করে না,
# explicitly সেট করতে হয় (এই লাইন ছাড়া build link-time এ ব্যর্থ
# হয়, যা এই sandbox এ প্রথম CMake build attempt এ ধরা পড়েছিল)।
add_library(civilos_solver_core OBJECT
    src/stiffness.cpp
    src/solver.cpp
    src/shell.cpp
)
set_target_properties(civilos_solver_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(civilos_solver_core PUBLIC include)
target_link_libraries(civilos_solver_core PUBLIC Eigen3::Eigen)

# --- Python extension module (এটাই FastAPI import করবে) ---
pybind11_add_module(civilos_solver src/bindings.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(civilos_solver PRIVATE include)
target_link_libraries(civilos_solver PRIVATE Eigen3::Eigen)

# --- Test executables (CI/build-time verification এর জন্য) ---
enable_testing()

add_executable(test_stiffness tests/test_stiffness.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_stiffness PRIVATE include)
target_link_libraries(test_stiffness PRIVATE Eigen3::Eigen)
add_test(NAME test_stiffness COMMAND test_stiffness)

add_executable(test_solver_integration tests/test_solver_integration.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_solver_integration PRIVATE include)
target_link_libraries(test_solver_integration PRIVATE Eigen3::Eigen)
add_test(NAME test_solver_integration COMMAND test_solver_integration)

add_executable(test_pin_release tests/test_pin_release.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_pin_release PRIVATE include)
target_link_libraries(test_pin_release PRIVATE Eigen3::Eigen)
add_test(NAME test_pin_release COMMAND test_pin_release)

add_executable(test_modal_analysis tests/test_modal_analysis.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_modal_analysis PRIVATE include)
target_link_libraries(test_modal_analysis PRIVATE Eigen3::Eigen)
add_test(NAME test_modal_analysis COMMAND test_modal_analysis)

add_executable(test_buckling_analysis tests/test_buckling_analysis.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_buckling_analysis PRIVATE include)
target_link_libraries(test_buckling_analysis PRIVATE Eigen3::Eigen)
add_test(NAME test_buckling_analysis COMMAND test_buckling_analysis)

add_executable(test_pdelta_analysis tests/test_pdelta_analysis.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_pdelta_analysis PRIVATE include)
target_link_libraries(test_pdelta_analysis PRIVATE Eigen3::Eigen)
add_test(NAME test_pdelta_analysis COMMAND test_pdelta_analysis)

add_executable(test_shell_element tests/test_shell_element.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_shell_element PRIVATE include)
target_link_libraries(test_shell_element PRIVATE Eigen3::Eigen)
add_test(NAME test_shell_element COMMAND test_shell_element)

add_executable(test_shell_integration tests/test_shell_integration.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_shell_integration PRIVATE include)
target_link_libraries(test_shell_integration PRIVATE Eigen3::Eigen)
add_test(NAME test_shell_integration COMMAND test_shell_integration)

add_executable(test_response_spectrum_analysis tests/test_response_spectrum_analysis.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_response_spectrum_analysis PRIVATE include)
target_link_libraries(test_response_spectrum_analysis PRIVATE Eigen3::Eigen)
add_test(NAME test_response_spectrum_analysis COMMAND test_response_spectrum_analysis)

add_executable(test_nonlinear_static_analysis tests/test_nonlinear_static_analysis.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_nonlinear_static_analysis PRIVATE include)
target_link_libraries(test_nonlinear_static_analysis PRIVATE Eigen3::Eigen)
add_test(NAME test_nonlinear_static_analysis COMMAND test_nonlinear_static_analysis)

add_executable(test_pushover_analysis tests/test_pushover_analysis.cpp $<TARGET_OBJECTS:civilos_solver_core>)
target_include_directories(test_pushover_analysis PRIVATE include)
target_link_libraries(test_pushover_analysis PRIVATE Eigen3::Eigen)
add_test(NAME test_pushover_analysis COMMAND test_pushover_analysis)
