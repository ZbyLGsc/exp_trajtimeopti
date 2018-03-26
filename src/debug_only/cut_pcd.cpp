#include <iostream>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/radius_outlier_removal.h>

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/registration/transforms.h>
#include <pcl/search/kdtree.h>
#include <pcl/io/pcd_io.h>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <fstream>

#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <math.h>

#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <random>

#include "visualization_msgs/MarkerArray.h"
#include "visualization_msgs/Marker.h"

using namespace std;

using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;
using namespace std;

pcl::PointCloud<pcl::PointXYZ>::Ptr CloudIn(new pcl::PointCloud<pcl::PointXYZ>());
int main (int argc, char** argv) {
        
      ros::init (argc, argv, "test_frame");
      ros::NodeHandle n( "~" );

      if (pcl::io::loadPCDFile<pcl::PointXYZ> ("input.pcd", *CloudIn) == -1) //* load the file input_less.pcd
        {
              PCL_ERROR ("Couldn't read file input.pcd \n");
              return (-1);
        }
      

      pcl::PointCloud<pcl::PointXYZ> Cloud_Less;
      pcl::PointXYZ pt_select;
      for (int i = 0; i < int(CloudIn->points.size()); i++) 
      {
          pt_select = CloudIn->points[i];
              
          if(fabs(pt_select.x) > 70.0 || fabs(pt_select.y) > 70.0 || pt_select.z < -7.5 )
                  continue;
                  
          Cloud_Less.points.push_back( pt_select);
      }
      
      pcl::PointCloud<pcl::PointXYZ> Cloud_DS;
      pcl::VoxelGrid<pcl::PointXYZ>  VoxelSampler;
      
      VoxelSampler.setLeafSize(0.1f, 0.1f, 0.1f);
      VoxelSampler.setInputCloud( Cloud_Less.makeShared() );      
      VoxelSampler.filter( Cloud_DS );       

      ROS_WARN("Check points number of pointCloud : %lu", Cloud_DS.points.size());

      Cloud_DS.width = Cloud_DS.points.size();
      Cloud_DS.height = 1;
      Cloud_DS.is_dense = true;
      pcl::io::savePCDFileASCII ("input_less.pcd", Cloud_DS);    
      


      ros::spin ();
}