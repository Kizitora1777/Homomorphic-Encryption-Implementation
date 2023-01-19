set(CMAKE_CXX_STANDARD 17)

find_package(helib)

add_executable(debug debug.cpp)

target_link_libraries(debug helib)
