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
#include <sensor_msgs/LaserScan.h>
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
#include <kdtree/kdtree.h>
#include "pcd_trajectory/dataType.h"

using namespace std;

class rrtPathFinder
{
	private:
			pcl::search::KdTree<pcl::PointXYZ> kdtreeForMap;
			pcl::search::KdTree<pcl::PointXYZ> kdtreeAddMap;
			pcl::search::KdTree<pcl::PointXYZ> kdtreeDelMap;
			//pcl::PointCloud<pcl::PointXYZ> CloudIn;
			
			kdtree * kdTree_; // dynamic light-weight Kd-tree, for organizing the nodes in the exploring tree

			vector<int>     pointIdxRadiusSearch;
			vector<float>   pointRadiusSquaredDistance;        
			
			// All nodes,            for visualize,  all nodes reach the target
			vector<NodePtr> NodeList, visNodeList,   EndList;
			
			// all nodes in the current path, for easy usage in root re-decleration 
			vector<NodePtr> PathList, invalidSet;

			// all nodes lost its pre-node connection, for checking new connection in each callbck
			vector<NodePtr> reCheckList;
			
			// record the current best node which can reach the target and has the lowest path cost
			NodePtr best_end_ptr;

			// record the root of the rapidly exploring tree
			NodePtr rootNode;

			// start point,   target point,  centroid point of the ellipsoide sampling region
			Point start_pt,     end_pt,                  inform_centroid,                                   commit_root;

			// useful counter                  // ctrl the size of the trash nodes cach, once it's full, remove them and rebuild the NN
			int node_id,       invalid_cnt,    cach_size; 

			// used for basic of sampling-based path finding
			float x_l, x_h, y_l, y_h, z_l, z_h, bias_l, bias_h, inlier_ratio, goal_ratio;
			
			// used for range query and search strategy
			float safety_margin, max_radius, search_margin, sample_range;
			
			// used for the informed sampling strategy
			float min_distance, best_distance, elli_l, elli_s, ctheta, stheta;   
			
			// FLAGs
			bool inform_sample; // indicate whether the path refine has been finished

			Eigen::MatrixXd Path;
			Eigen::VectorXd Radius;
			
			random_device rd;
			default_random_engine eng;

			uniform_real_distribution<float>  rand_rho = uniform_real_distribution<float>(0.0, 1.0);  // random distribution for generating samples inside a unit circle
			uniform_real_distribution<float>  rand_phi = uniform_real_distribution<float>(0.0, 2 * M_PI);
			uniform_real_distribution<float>  rand_x,    rand_y,    rand_z,   rand_bias; // basic random distributions for generating samples, in all feasible regions
			uniform_real_distribution<float>  rand_x_in, rand_y_in, rand_z_in; // random distribution, especially for generating samples inside the local map's boundary

			vector<Point> sampleSet;

			int max_samples;
			int nr_nodes;

	public:
			bool path_find_state, refine_status, commitEnd;

			rrtPathFinder( );

	        ~rrtPathFinder();
			
			void setParam( float safety_margin_, float search_margin_, float max_radius_, float sample_range_ );

			void reset();

			void resetRoot(Point commitTarget);

			void costRecast(double costReduction, Point target);

			void setInput(pcl::PointCloud<pcl::PointXYZ> CloudIn);
			
			void rcvAddMap(pcl::PointCloud<pcl::PointXYZ> CloudAdd);

			void rcvDelMap(pcl::PointCloud<pcl::PointXYZ> CloudDel);

			void setPt( Point startPt, Point endPt, float xl, float xh, float yl, float yh, float zl, float zh, float biasl, float biash, 
						float localRange, int max_iter, float sample_portion, float goal_portion );

			void setStartPt( Point startPt, Point endPt);

			void clearBranchW(NodePtr node_delete); // weak clear:    clear branches while avoid deleting the nodes on the best path
 
			void clearBranchS(NodePtr node_delete); // strong clear:  clear all nodes on a branch

	        inline float getDis(const NodePtr node1, const NodePtr node2);

			inline float getDis(const NodePtr node1, const Point pt);

			inline float getDis(const Point p1, const Point p2);

			inline float getDisL1(const NodePtr node1, const Point pt);

			float radius_search(Point search_Pt);//, pcl::search::KdTree<pcl::PointXYZ>::Ptr Tree);
			
			float radius_search_add(Point search_Pt);

			Point genSample();
			
			NodePtr findNearstVertex( Point pt_sample );
			
			NodePtr findNearst( Point pt_sample );

			NodePtr genNewNode( Point pt_sample, NodePtr node_nearst_ptr );

			void UpdateHeuristicRegion( NodePtr update_end_node );

			bool checkTrajInvalid(Point traj_pt);

			bool check_end( NodePtr ptr );

			bool check_boundary_reach( NodePtr ptr );

			bool check_end_AStar(NodePtr ptr);
			
			void AddtoGraph(NodePtr new_node);			

			void TreeSparsify( NodePtr newPtr );

			inline bool CheckConnect( float dis, NodePtr node_1, NodePtr node_2 );

			inline int  CheckRelation( float dis, NodePtr node_1, NodePtr node_2 );

			void TreeRewire( NodePtr node_new_ptr, NodePtr node_nearst_ptr );
			
			void TreeExpansion( NodePtr node_new_ptr, NodePtr node_nearst_ptr );

			pair<Eigen::MatrixXd, Eigen::VectorXd> getPath();

			void AStarPath(NodePtr p_start, NodePtr p_target);

			void RRTpathFind( float path_time_limit);

			vector<NodePtr> getTree();

			vector<Point> getSamples(){return sampleSet;}

			void RRTpathRefine( float refine_limit, bool isPrint );

			void RRTpathReEvaluate( float revaluate_limit, bool isPrint);

			bool RefineStatus();

			void removeInvalid();
			
			void TreeDestruct();

			void Nodefree( NodePtr ptr ){ delete ptr; }

			void ValidateNodes();

			inline bool insideEllipsoid( float* pos);

			float reCheckRadius(Point pt);

			float reCheckAddRadius(Point pt);

			int ReCheckConnection(float new_radius, float old_radius);

			bool CheckValidEnd(NodePtr endPtr);

			void parentForgetChild(NodePtr ptrParent, NodePtr nodeptr);

			bool isSuccessor(NodePtr curPtr, NodePtr nearPtr);

			int Repair();

			void tracePath();

			void treeRepair(double repair_limit, vector< pair<Point, double> > evaFailList);
};