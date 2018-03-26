
// ROS
#include <ros/ros.h>

// PCL
#include <pcl/filters/voxel_grid.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

// STD
#include <iostream>

// Grad planner
#include "grad_traj_optimizer.h"

using namespace std;

// definition
ros::Publisher _sdf_pub;
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2& pointcloud_map)
{
    pcl::PointCloud<pcl::PointXYZ> CloudIn;
    pcl::fromROSMsg(pointcloud_map, CloudIn);

    // filter the map
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);

    pcl::VoxelGrid<pcl::PointXYZ> voxel_filter;
    voxel_filter.setLeafSize(0.2f, 0.2f, 0.2f);
    voxel_filter.setInputCloud(CloudIn.makeShared());
    voxel_filter.filter(*cloud_filtered);

    // ----------- build signed distance field for gradient planner
    // sdf collision map parameter
    const double resolution = 0.2;
    const double x_size = 30.0;
    const double z_size = 5.0;
    double y_size = 30.0;
    Eigen::Translation3d origin_translation(-15.0, -15.0, 0.0);
    Eigen::Quaterniond origin_rotation(1.0, 0.0, 0.0, 0.0);
    const Eigen::Isometry3d origin_transform = origin_translation * origin_rotation;
    const std ::string frame = "/map";

    // create map
    sdf_tools ::COLLISION_CELL cell;
    cell.occupancy = 0.0;
    cell.component = 0;
    const sdf_tools ::COLLISION_CELL oob_cell = cell;
    sdf_tools ::CollisionMapGrid collision_map(origin_transform, frame, resolution, x_size, y_size,
                                               z_size, oob_cell);

    // add point cloud to sdf map
    sdf_tools::COLLISION_CELL obstacle_cell(1.0);
    std::vector<pcl::PointXYZ, Eigen::aligned_allocator<pcl::PointXYZ>> pts =
        cloud_filtered->points;
    for(int i = 0; i < int(pts.size()); ++i)
    {
        if((pts[i].x < 14.9 && pts[i].x > -14.9) && (pts[i].y < 14.9 && pts[i].y > -14.9) &&
           (pts[i].z < 4.5 && pts[i].z > 0.2))
        {
            collision_map.Set(pts[i].x, pts[i].y, pts[i].z, obstacle_cell);
            std::cout << setprecision(2) << i << " th pt x:" << pts[i].x << " ," << pts[i].y << " ,"
                      << pts[i].z << std::endl;
        }
    }

    // visualize the collision map
    std_msgs::ColorRGBA collision_color;
    collision_color.r = 0.0;
    collision_color.g = 0.0;
    collision_color.b = 1.0;
    collision_color.a = 0.8;

    std_msgs::ColorRGBA free_color, unknown_color;
    unknown_color.a = free_color.a = 0.0;

    visualization_msgs::Marker collision_map_marker =
        collision_map.ExportForDisplay(collision_color, free_color, unknown_color);
    collision_map_marker.ns = "collision_map";
    collision_map_marker.id = 1;

    _sdf_pub.publish(collision_map_marker);
    cout << "Point cloud published!" << endl;

    // Build the signed distance field
    float oob_value = INFINITY;
    std::pair<sdf_tools::SignedDistanceField, std::pair<double, double>> sdf_with_extrema =
        collision_map.ExtractSignedDistanceField(oob_value);
    sdf_tools::SignedDistanceField sdf = sdf_with_extrema.first;

    // test sdf
    for(int i = -10; i < 10; ++i)
        for(int j = -10; j < 10; ++j)
            for(int k = 0; k < 5; ++k)
            {
                std::vector<double> location_gradient_query =
                    sdf.GetGradient(double(i), double(j), double(k), true);
                cout << "pt: " << i << ", " << j << ", " << k << setprecision(2)
                     << "grad: " << location_gradient_query[0] << ", " << location_gradient_query[1]
                     << ", " << location_gradient_query[2] << endl;
            }
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "test_any");
    ros::NodeHandle n("~");

    // subscriber to pcl point cloud
    ros::Subscriber pc_sub = n.subscribe("/random_map/AllMap", 1, rcvPointCloudCallBack);
    _sdf_pub = n.advertise<visualization_msgs::Marker>("sdf_visualization", 1, true);

    ros::Duration(1.0).sleep();

    ros::spin();

    return 0;
}