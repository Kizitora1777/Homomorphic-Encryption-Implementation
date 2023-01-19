set(CMAKE_CXX_STANDARD 17)

find_package(helib)

add_executable(main main.cpp)

target_link_libraries(main helib)
