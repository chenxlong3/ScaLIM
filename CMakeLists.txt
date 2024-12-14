cmake_minimum_required(VERSION 3.5)		# Cmake Version

project(ScaLIM) # Project Name
include_directories(include)
find_package(Boost COMPONENTS system filesystem)

set(CMAKE_BUILD_TYPE Release)
add_compile_options(-O3)
# add_executable(format_graph exec_files/main.cpp include/SFMT/dSFMT/dSFMT.c)
add_executable(eval_by_mc exec_files/evaluate.cpp include/SFMT/dSFMT/dSFMT.c)
add_executable(ScaLIM exec_files/ScaLIM.cpp include/SFMT/dSFMT/dSFMT.c)

# target_link_libraries(MCGreedy_only ${Boost_LIBRARIES})

target_compile_features(ScaLIM PRIVATE cxx_std_11)