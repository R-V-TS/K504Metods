cmake_minimum_required(VERSION 3.14)
project(Clear_DCT)

#list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(OpenCV REQUIRED)
# find_package(Boost REQUIRED)
#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs REQUIRED)
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "-pthread -lpython2.7 -lboost_system")

add_executable(Clear_DCT
        src/CPU/main.cpp src/CPU/DCT_FILT.cpp src/CPU/ImageQuaolity.cpp)

target_link_libraries(Clear_DCT ${OpenCV_LIBS})

#target_include_directories(ImageProcessing_exe PRIVATE
#        /usr/include/python2.7/
#        ${Boost_LIBRARIES}
#        ${Boost_LIBRARY_DIR})