cmake_minimum_required(VERSION 2.6)

project(Aberration_Correction)

add_executable(chromaberrat abberation.h centers.h correction.h distortion.cpp distortion.h image.cpp image.h LM.h main_centering.cpp matrix.cpp matrix.h misc.h misc.cpp numerics.cpp numerics.h pgm_io.cpp pgm_io.h spline.cpp spline.h vector.cpp)
