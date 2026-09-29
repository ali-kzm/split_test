cmake_minimum_required(VERSION 3.16)
project(split_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(split_mesher
    src/quad_levelset_mesher.cpp
)
target_include_directories(split_mesher PUBLIC include)
target_compile_options(split_mesher PRIVATE -Wall -Wextra -Wpedantic)

add_executable(split_demo src/main.cpp)
target_link_libraries(split_demo PRIVATE split_mesher)

include(CTest)
if(BUILD_TESTING)
    add_executable(test_mesher tests/test_mesher.cpp)
    target_link_libraries(test_mesher PRIVATE split_mesher)
    add_test(NAME mesher_tests COMMAND test_mesher)
endif()
