cmake_minimum_required(VERSION 3.16.3)
set(CMAKE_CXX_STANDARD 17)
project(Thicken2)


include_directories("/opt/homebrew/include")

include_directories(/usr/local/include/osqp/)

include_directories(/opt/homebrew/include/eigen3)

include_directories(/opt/osqp-eigen/include/)

LINK_DIRECTORIES(/usr/local/lib/)

find_package(OsqpEigen REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})



file (GLOB CppFile
        MeshKernel/*
        GLKLib/*
        Eigen/*
        )
add_executable(Thicken2 main.cpp lib_impl.cpp DSU.cpp ${CppFile})

target_link_libraries(Thicken2
        /opt/homebrew/lib/libgmp.dylib
        -lpthread
        /opt/homebrew/lib/libmpfr.dylib
        ${Boost_LIBRARY_DIRS}
        OsqpEigen::OsqpEigen
        )


