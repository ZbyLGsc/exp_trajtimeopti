#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include <ros/ros.h>
#include <iostream>

int main(int argc, char** argv)
{
  ros::init(argc, argv, "waypoint_publisher");
  ros::NodeHandle node("~");

  ros::Publisher wp_pub = node.advertise<nav_msgs::Path>("/waypoint_generator/waypoints", 1);
  ros::Duration(1.0).sleep();

  while(ros::ok())
  {
    //--------------get coordinate--------------------------------------------
    std::cout << "Please input coordinate!" << std::endl;

    int x = 0, y = 0, z = 0;
    std::cin >> x >> y >> z;

    std::cout << "coordinate:(" << x << "," << y << "," << z << ")\n";

    //--------------publish message to trigger optimization-------------------
    nav_msgs::Path path;
    path.header.frame_id = "/map";
    path.header.stamp = ros::Time::now();

    geometry_msgs::PoseStamped pose;
    pose.header.frame_id = "/map";

    pose.pose.position.x = x;
    pose.pose.position.y = y;
    pose.pose.position.z = z;

    pose.pose.orientation.w = 1;
    pose.pose.orientation.x = 0;
    pose.pose.orientation.y = 0;
    pose.pose.orientation.z = 0;

    path.poses.push_back(pose);

    wp_pub.publish(path);

    ROS_INFO("Waypoint published!");
    ros::Duration(0.5).sleep();
  }

  return 0;
}