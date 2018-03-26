#include <iostream>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>

#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <math.h>

#include <time.h>
#include <sys/time.h>
#include <random>
#include <nav_msgs/Path.h>
#include "visualization_msgs/MarkerArray.h"
#include "visualization_msgs/Marker.h"

using namespace std;

ros::Publisher _path_new;
ros::Publisher _traj_new;
ros::Publisher _path_ekf_new;

ros::Subscriber _path_sub;
ros::Subscriber _traj_sub;
ros::Subscriber _path_ekf;

visualization_msgs::MarkerArray path_new;

nav_msgs::Path path_ekf_new;

void pathCallBack(visualization_msgs::MarkerArray path_vis)
{     
       for (auto & mk: path_new.markers) 
        mk.action = visualization_msgs::Marker::DELETE;

      _path_new.publish(path_new);
      path_new.markers.clear();

      visualization_msgs::Marker mk;
      mk.header.frame_id = "map";
      mk.header.stamp = ros::Time::now();
      mk.ns = "convert_plot/path";
      mk.type = visualization_msgs::Marker::SPHERE;
      mk.action = visualization_msgs::Marker::ADD;
      mk.pose.orientation.x = 0.0;
      mk.pose.orientation.y = 0.0;
      mk.pose.orientation.z = 0.0;
      mk.pose.orientation.w = 1.0;
      mk.color.a = 0.4;
      mk.color.r = 0.0;
      mk.color.g = 1.0;
      mk.color.b = 0.0;

      for(int i = 0; i < int(path_vis.markers.size()); i++){
          mk.id = i;
          visualization_msgs::Marker ori;
          ori = path_vis.markers[i];

          mk.pose = ori.pose; 
          mk.scale = ori.scale;
          path_new.markers.push_back(mk);
    }

    _path_new.publish(path_new);

}

void trajCallBack(visualization_msgs::Marker traj_vis)
{
         visualization_msgs::Marker traj_new;

        traj_new.header.stamp       = ros::Time::now();
        traj_new.header.frame_id    = "map";

        traj_new.ns = "convert_plot/new_trajectory";
        traj_new.id = 0;
        traj_new.type = visualization_msgs::Marker::SPHERE_LIST;
        traj_new.action = visualization_msgs::Marker::ADD;
        traj_new.scale.x = 0.15;
        traj_new.scale.y = 0.15;
        traj_new.scale.z = 0.15;
        traj_new.pose.orientation.x = 0.0;
        traj_new.pose.orientation.y = 0.0;
        traj_new.pose.orientation.z = 0.0;
        traj_new.pose.orientation.w = 1.0;
        traj_new.color.r = 1.0;
        traj_new.color.g = 0.0;
        traj_new.color.b = 0.0;
        traj_new.color.a = 1.0;

        geometry_msgs::Point pt_new;

        for(int i = 0; i < int(traj_vis.points.size()); i++ ){
            geometry_msgs::Point pt = traj_vis.points[i];
            pt_new.x = pt.x;
            pt_new.y = pt.y;
            pt_new.z = pt.z;
                
            traj_new.points.push_back(pt_new);
        }

        _traj_new.publish(traj_new);
}

/*void pathekfCallBack(nav_msgs::Path path_ekf)
{   
    cout<<path_ekf.poses.size()<<endl;
    
    if(path_ekf.poses.size() == 0)
      return;

    path_ekf_new = path_ekf;
    _path_ekf_new.publish(path_ekf_new);    
}*/

int main (int argc, char** argv) {
        
      ros::init (argc, argv, "ConvertPlot");
      ros::NodeHandle n( "~" );
      
      _path_sub = 
            n.subscribe( "/pcd_trajectory_node/path", 1, pathCallBack );
      _traj_sub  = 
            n.subscribe( "/pcd_trajectory_node/trajectory_vis",1, trajCallBack);
      /*_path_ekf  = 
            n.subscribe("/ekf_node/path_ekf",1, pathekfCallBack);*/

      _path_new = 
            n.advertise<visualization_msgs::MarkerArray>( "path_new", 2 );     
      _traj_new = 
            n.advertise<visualization_msgs::Marker>( "traj_new", 2 );
      /*_path_ekf_new = 
            n.advertise<nav_msgs::Path>("ekf_path_new",10);
*/
      ros::spin ();
}
