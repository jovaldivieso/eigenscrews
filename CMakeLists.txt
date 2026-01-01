cmake_minimum_required(VERSION 3.10)
project(eigenscrews VERSION 1.0.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
find_package(Eigen3 REQUIRED NO_MODULE)

# --- Library Target ---
add_library(${PROJECT_NAME} src/screws.cpp)

# Define include directories for the library itself and consumers
target_include_directories(${PROJECT_NAME} PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${EIGEN3_INCLUDE_DIRS}
)

# Link Eigen to the library
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

# --- Options ---
option(BUILD_EXAMPLES "Build example applications" ON)
option(BUILD_TESTS "Build unit tests" OFF)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
