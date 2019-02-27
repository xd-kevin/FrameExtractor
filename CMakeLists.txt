cmake_minimum_required(VERSION 3.13)
project(FrameExtractor)

set(CMAKE_CXX_STANDARD 14)

add_definitions(-std=c++11)

set(INC_DIR /usr/local/include)
set(LINK_DIR /usr/local/lib)

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${INC_DIR}")
LINK_DIRECTORIES(${LINK_DIR})


add_executable(FrameExtractor src/extractor.cpp src/video_parser.cpp)
target_link_libraries(FrameExtractor opencv_core opencv_highgui opencv_imgproc opencv_imgcodecs opencv_videoio)