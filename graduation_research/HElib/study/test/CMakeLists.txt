set(CMAKE_CXX_STANDARD 17)

find_package(helib)

add_executable(test test.cpp)

target_link_libraries(test helib)
