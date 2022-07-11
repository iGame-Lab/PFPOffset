cmake_minimum_required(VERSION 3.16.3)
project(Thicken2)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "linux")
    include_directories("/usr/include")
#    SET(Boost_INCLUDE_DIRS "D:\\boost_1_79_0\\")
#    SET(Boost_LIBRARY_DIRS "D:\\boost_1_79_0\\lib64-msvc-14.3")
#
#    include_directories("C:\\Program Files (x86)\\CGAL\\include")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    SET(Boost_INCLUDE_DIRS "D:\\boost_1_79_0\\")
    SET(Boost_LIBRARY_DIRS "D:\\boost_1_79_0\\lib64-msvc-14.3")

    include_directories("C:\\Program Files (x86)\\CGAL\\include")
endif()


include_directories(/usr/local/include/osqp/)

include_directories(/opt/osqp-eigen/include/)

LINK_DIRECTORIES(/usr/local/lib/)

find_package(OsqpEigen REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})


set(CMAKE_CXX_STANDARD 14)
file (GLOB CppFile
        MeshKernel/*
        GLKLib/*
        Eigen/*
        )
add_executable(Thicken2 main.cpp lib_impl.cpp DSU.cpp ${CppFile})

target_link_libraries(Thicken2
        -lgmp
        -lpthread
        ${Boost_LIBRARY_DIRS}
        OsqpEigen::OsqpEigen
        )


