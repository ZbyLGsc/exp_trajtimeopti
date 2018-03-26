cmake_minimum_required(VERSION 2.8.3)
project(pcd_trajectory)

set(CMAKE_VERBOSE_MAKEFILE "true")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
# set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -DEIGEN_DONT_PARALLELIZE")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall")

find_package(catkin REQUIRED COMPONENTS
    roscpp
    std_msgs
    nav_msgs
    sensor_msgs
    quadrotor_msgs
    tf
    kdtree
    arc_utilities
    grad_traj_optimization
    )

find_package(PCL 1.7 REQUIRED)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Eigen3  REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS kdtree
)

link_directories( /home/zby/Downloads/software/mosek/8/tools/platform/linux64x86/bin)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS}
    )

include_directories(/home/zby/workspaces/traj_ws/src/grad_traj_optimization/third_party/sdf_tools/include)
include_directories(/home/zby/workspaces/traj_ws/src/grad_traj_optimization/include)


link_directories(${PCL_LIBRARY_DIRS})


add_executable (random_map 
                src/random_map.cpp )
target_link_libraries (random_map
                       ${catkin_LIBRARIES}
                       ${PCL_LIBRARIES})


add_executable( time_optimal_traj 
                src/time_optimal_node.cpp
                src/trajectory_generator.cpp 
                src/trajectory_generator_lite.cpp 
                src/trajectory_generator_socp_lite.cpp 
                src/bezier_base.cpp
                src/rrgPathFinder.cpp
                src/time_optimizer.cpp
                src/qp_generator.cpp
                /home/zby/workspaces/traj_ws/src/grad_traj_optimization/src/grad_traj_optimizer.cpp
)
target_link_libraries( time_optimal_traj
                       mosek64
                       ${catkin_LIBRARIES}
                       ${PCL_LIBRARIES}
                       /home/zby/workspaces/traj_ws/src/grad_traj_optimization/lib/libnlopt.so.0.9.0
                       /home/zby/workspaces/traj_ws/src/grad_traj_optimization/third_party/sdf_tools/lib/libsdf_tools.so
)


# add_executable( tan src/test_any.cpp src/qp_generator.cpp)
# target_link_libraries(tan ${catkin_LIBRARIES}
# ${PCL_LIBRARIES} /home/zby/workspaces/traj_ws/src/grad_traj_optimization/third_party/sdf_tools/lib/libsdf_tools.so
# /home/zby/workspaces/traj_ws/src/grad_traj_optimization/lib/libnlopt.so.0.9.0)

set(PROJECT_SOURCE_DIR)

# add_executable(wp_pub src/waypoint_publisher.cpp)
# target_link_libraries(wp_pub ${catkin_LIBRARIES})

# add_executable( time_optimal 
#                 src/test_time_opti.cpp
#                 src/trajectory_generator.cpp 
#                 src/trajectory_generator_lite.cpp 
#                 src/trajectory_generator_socp_lite.cpp 
#                 src/bezier_base.cpp
#                 src/rrgPathFinder.cpp
#                 src/time_optimizer.cpp
#                 src/qp_generator.cpp
# )
# target_link_libraries( time_optimal
#                        mosek64
#                        ${catkin_LIBRARIES}
#                        ${PCL_LIBRARIES}
#                        /usr/local/lib/libnlopt.so.0.9.0
#                        /home/zby/workspaces/traj_ws/src/grad_traj_optimization/third_party/sdf_tools/lib/libsdf_tools.so
# )