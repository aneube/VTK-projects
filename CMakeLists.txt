cmake_minimum_required(VERSION 2.6)
project(VTK-projects)
set(VTK_DIR "/home/neu016/Dev/VTK6.3/build")
find_package(VTK REQUIRED)

SUBDIRS(DistancePolyData)
