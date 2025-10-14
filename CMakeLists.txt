cmake_minimum_required(VERSION 3.10)
project(graphene VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Compiler flags
add_compile_options(-Wall -Wextra -O2)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
file(GLOB SOURCES "src/*.c")

# Create static library
add_library(graphene STATIC ${SOURCES})

# Installation
install(TARGETS graphene
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(DIRECTORY include/
    DESTINATION include/graphene
    FILES_MATCHING PATTERN "*.h"
)

# Tests
enable_testing()

file(GLOB TEST_SOURCES "tests/*.c")
foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} graphene)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

# Example app (optional)
# add_subdirectory(apps)

