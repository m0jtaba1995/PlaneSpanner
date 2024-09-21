cmake_minimum_required(VERSION 3.18)
project(PlaneSpanner)

set(CMAKE_CXX_STANDARD 20)
set(CGAL_DATA_DIR "$HOME/CGAL-5.6/data")
find_package(CGAL REQUIRED)

add_executable(PlaneSpanner main.cpp)
target_link_libraries(PlaneSpanner PRIVATE CGAL::CGAL)