cmake_minimum_required(VERSION 2.6)

project(Aberration_Correction)

add_executable(chromaberrat abberation.h centers.h correction.h demo_lib_sift.cpp demo_lib_sift.h distortion.cpp distortion.h domain.cpp domain.h filter.cpp filter.h image.cpp image.h library.cpp library.h LM.h LWImage.cpp LWImage.h main_centering.cpp match.cpp match.h matrix.cpp matrix.h misc.h misc.cpp numerics.cpp numerics.h pgm_io.cpp pgm_io.h spline.cpp spline.h vector.cpp)