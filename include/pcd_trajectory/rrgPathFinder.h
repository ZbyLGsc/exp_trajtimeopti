#include <iostream>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/octree/octree_search.h>
#include <pcl/octree/octree.h>

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/registration/transforms.h>
#include <pcl/search/kdtree.h>
#include <pcl/io/pcd_io.h>
#include <pcl/search/impl/kdtree.hpp>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <fstream>

#include <Eigen/Eigen>
#include <math.h>

#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <random>

#include "pcd_trajectory/rrgDataType.h"

using namespace std;

class rrgPathFinder
{

	private:
			pcl::search::KdTree<pcl::PointXYZ> kdtreeForMap;
			//pcl::PointCloud<pcl::PointXYZ> CloudIn;

			vector<int>     pointIdxRadiusSearch;
			vector<float>   pointRadiusSquaredDistance;        
			vector<NodePtr> NodeList;
			vector<NodePtr> EndList;
			NodePtr best_node;
			NodeVec NodeMap[14][14];
			NodePtr rootPtr;

			NodePtr p_Astar_target;
			NodePtr p_Astar_start;

			Point start_pt;
			Point end_pt;                         

			double x_l, x_h, y_l, y_h, z_l, z_h, bias_l, bias_h;

			Eigen::MatrixXd Path;
			Eigen::VectorXd Radius;
			
			random_device rd;
			default_random_engine eng;

			uniform_real_distribution<double>  rand_x;
			uniform_real_distribution<double>  rand_y;
			uniform_real_distribution<double>  rand_z;
			uniform_real_distribution<double>  rand_bias;

	public:

			bool path_find_state;
			
			void reset();

			rrgPathFinder(double xl, double xh, double yl, double yh, double zl, double zh, double biasl, double biash);

	        ~rrgPathFinder();
			
			void setInput(pcl::PointCloud<pcl::PointXYZ> CloudIn);

			void setPt( Point startPt, Point endPt);
	        /*void setInit(double x_l, double x_h, double y_l, double y_h, double z_l, double z_h, double bias_l, double bias_h);*/

	        double getDis(const NodePtr node1, const NodePtr node2);

			double getDis(const NodePtr node1, const Point pt);

			double getDisL1(const NodePtr node1, const Point pt);

			double radius_search(Point search_Pt);//, pcl::search::KdTree<pcl::PointXYZ>::Ptr Tree);

			Point genSample();

			NodePtr findNearst( Point pt_sample );

			NodePtr findNearst( Point pt_sample, vector<NodePtr> NodeList );

			NodePtr genNewNode( Point pt_sample, NodePtr node_nearst_ptr );

			bool check_end( NodePtr ptr );

			bool check_end_AStar(NodePtr ptr);
			
			void Add_new_node(NodePtr ptr_node_new);

			void Clear_NodeMap();

			bool check_no_Empty(Point pt_sample);

			bool CheckConnect( double dis, NodePtr node_1, NodePtr node_2 );

			void AddtoGraph( NodePtr node_new_ptr, NodePtr node_nearst_ptr);

			Eigen::MatrixXd getPath();

			Eigen::VectorXd getRadius();

			void AStarPath(NodePtr p_start, NodePtr p_target);

			void RRGpathFind();
};