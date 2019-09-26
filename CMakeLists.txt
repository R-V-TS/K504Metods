cmake_minimum_required(VERSION 3.14)

project(ImageProjGPU)

#list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(OpenCV REQUIRED)
# find_package(Boost REQUIRED)
#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs REQUIRED)
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "-pthread -lpython2.7 -lboost_system")

add_executable(ImageProjGPU
        src/main.cpp src/RAWImage.cpp src/CPU/DCT.cpp src/GPU/GPU_DCT.cpp src/DCT_Matrices.h  src/CPU/ImageQuaolity.cpp )

set_source_files_properties(src/GPU/GPU_DCT.cpp PROPERTIES LANGUAGE CUDA)

target_link_libraries(ImageProjGPU ${OpenCV_LIBS})

project(ImageProjCPU)

#list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(OpenCV REQUIRED)
# find_package(Boost REQUIRED)
#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs REQUIRED)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
#set(CMAKE_CXX_FLAGS "-pthread -lpython2.7 -lboost_system")

add_executable(ImageProjCPU
        src/main.cpp src/RAWImage.cpp src/CPU/DCT.cpp src/GPU/GPU_DCT_stub.cpp src/CPU/ImageQuaolity.cpp src/CPU/BM3D.cpp)


target_link_libraries(ImageProjCPU ${OpenCV_LIBS})