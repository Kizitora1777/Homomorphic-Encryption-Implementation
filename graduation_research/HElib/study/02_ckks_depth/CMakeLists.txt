set(CMAKE_CXX_STANDARD 17)

find_package(helib)

add_executable(02_ckks_depth 02_ckks_depth.cpp)

target_link_libraries(02_ckks_depth helib)
