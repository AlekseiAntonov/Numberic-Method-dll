cmake_minimum_required(VERSION 3.15)

project(MATRIXDLL)

set(CMAKE_C++_STANDARD 11)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(Matrix SHARED Matrix.h Matrix_Protected.cpp Matrix_Public.cpp
        Vectors/Vector_Private.cpp Vectors/Vector_Public.cpp Vectors/Vector.h
        Matrix_Wrap.h Matrix_Wrap.cpp)
add_executable(MATRIXDLL main.cpp Matrix.h Matrix_Protected.cpp Matrix_Public.cpp
        Vectors/Vector_Private.cpp Vectors/Vector_Public.cpp Vectors/Vector.h)