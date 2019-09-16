cmake_minimum_required(VERSION 3.14)

project(ImageProj)

#list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(OpenCV REQUIRED)
# find_package(Boost REQUIRED)
#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs REQUIRED)
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "-pthread -lpython2.7 -lboost_system")

add_executable(ImageProj
        src/main.cpp src/RAWImage.cpp src/CPU/DCT.cpp src/GPU/GPU_DCT.cpp src/DCT_Matrices.h)

set_source_files_properties(src/GPU/GPU_DCT.cpp PROPERTIES LANGUAGE CUDA)

target_link_libraries(ImageProj ${OpenCV_LIBS})

project(Clear_DCT)

#list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(OpenCV REQUIRED)
# find_package(Boost REQUIRED)
#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs REQUIRED)
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "-pthread -lpython2.7 -lboost_system")

add_executable(Clear_DCT
        src/CPU/main.cpp src/CPU/DCT.cpp src/CPU/ImageQuaolity.cpp)

target_link_libraries(Clear_DCT ${OpenCV_LIBS})



set(CMAKE_CUDA_COMPILER /usr/local/cuda-10.1/bin/nvcc)

project(cuda_ CXX CUDA)

find_package(OpenCV REQUIRED)
set(CMAKE_CXX_STANDARD 11)

#file(GLOB app_stc "tiff-4.0.9/libtiff/*.c")

add_executable(cuda_
        src/GPU/matrixMultiply.cpp src/GPU/cuda_heder.h src/GPU/testDCT.cpp src/RAWImage.cpp src/GPU/GPU_Functions.h src/GPU/GPU_DCT.cpp src/GPU/GPU_DCT_stub.cpp)

target_link_libraries(cuda_ ${OpenCV_LIBS})
set_source_files_properties(src/GPU/matrixMultiply.cpp PROPERTIES LANGUAGE CUDA)

target_compile_features(cuda_ PUBLIC cxx_std_11)

#target_include_directories(ImageProcessing_exe PRIVATE
#        /usr/include/python2.7/
#        ${Boost_LIBRARIES}
#        ${Boost_LIBRARY_DIR})