//#define PCL_NO_PRECOMPILE
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
#include <pcl/search/impl/kdtree.hpp>

#include "tf/tf.h"
#include "tf/transform_datatypes.h"
#include <tf/transform_broadcaster.h>
#include "visualization_msgs/MarkerArray.h"
#include "visualization_msgs/Marker.h"
#include "pcd_trajectory/mosek.h"
#include "pcd_trajectory/trajectory_generator_lite.h"
#include "pcd_trajectory/trajectory_generator_socp_lite.h"
#include "pcd_trajectory/path_finder.h"
#include "pcd_trajectory/backward.hpp"
#include "pcd_trajectory/dataType.h"
#include "pcd_trajectory/bezier_base.h"
#include "quadrotor_msgs/PositionCommand.h"
#include "quadrotor_msgs/PolynomialTrajectory.h"

namespace backward {
backward::SignalHandling sh;
}

using namespace std;
using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;

visualization_msgs::MarkerArray path_vis;
visualization_msgs::MarkerArray ctrlPt_vis;
visualization_msgs::MarkerArray locTarget_vis;

ros::Publisher _traj_vis_pub;
ros::Publisher _checkTraj_vis_pub;
ros::Publisher _path_vis_pub;
ros::Publisher _ctrlPt_vis_pub;
ros::Publisher _locTarget_vis_pub;
ros::Publisher _traj_pub;
ros::Publisher _sensing_range_pub;
ros::Publisher _vis_rrt_star_pub;
ros::Publisher _vis_samples_pub;

ros::Subscriber _odom_sub;
ros::Subscriber _dest_pts_sub;
ros::Subscriber _map_sub;
ros::Subscriber _add_map_sub;
ros::Subscriber _del_map_sub;
ros::Subscriber _cmd_sub;

nav_msgs::Odometry _odom;
nav_msgs::Path _waypoints;
bool _has_odom = false;
const size_t _odom_queue_size = 200;
deque<nav_msgs::Odometry> _odom_queue;
vector<double> _state;

vector<NodePtr> nodesAll;

quadrotor_msgs::PolynomialTrajectory _traj;
uint32_t _traj_id = 0;
float _vis_traj_width = 0.15; 

ros::Time _start_time = ros::TIME_MAX;
ros::Time init_time;
double _MAX_Vel,_MAX_Acc, _Vel, _Acc, _eps;

float _refine_portion ,time_limit, refine_limit, revaluate_limit, _path_find_limit, _sample_portion, _goal_portion;
float _safety_margin, _search_margin, _max_radius, _sensing_range, _sensing_rate;
int  _max_samples;
int  _minimize_order, _poly_order_min, _poly_order_max;
int  _segment_num;
vector<int> _poly_orderList;

double _time_commit;

bool traFinish     = false;
bool traHalt       = false;
bool targetReach   = false;
bool targetReceive = false;
bool missionFinish = false;

float _x_l, _x_h, _y_l, _y_h, _z_l, _z_h;  // For random map simulation : map boundary
float bias_l = 0.0, bias_h = 1.0; // sampling goal bias 

Eigen::MatrixXd _PolyCoeff;
Eigen::VectorXd _Time;
Eigen::VectorXd _Scale;
Eigen::VectorXd _Radius;
Eigen::MatrixXd _Path;

Point startPt  = { 0.0, 0.0, 1.0}; // For outdoor CYT model .
Point startVel = { 0.0, 0.0, 0.0};
Point startAcc = { 0.0, 0.0, 0.0};
Point endPt;      
Point commitTarget;

vector<MatrixXd> _MQMList, _FMList;
vector<VectorXd> _CList, _CvList, _CaList;

/*MatrixXd _MQM; // M' * Q * M in each block of the objective. No scale, only meta elements .
VectorXd _C;   // Position coefficients vector, used to record all the pre-compute 'n choose k' combinatorial for the bernstein coefficients .
VectorXd _C_v; // Velocity coefficients vector.
VectorXd _C_a; // Acceleration coefficients vector.
*/
Vector3d getStateFromBezier(const Eigen::MatrixXd & polyCoeff,  double t_now, int seg_now );

VectorXd getFullStateFromBezier(const Eigen::MatrixXd & polyCoeff,  double t_now, int seg_now );

Eigen::VectorXd TimeAllocate(Eigen::MatrixXd Path, Eigen::VectorXd Radius, Eigen::Vector3d iniVel );

quadrotor_msgs::PolynomialTrajectory getBezierTraj();

void visBezierTrajectory(Eigen::MatrixXd polyCoeff);

void visCtrlPoint(Eigen::MatrixXd polyCoeff);

void visPath(Eigen::MatrixXd path, Eigen::VectorXd radius);

void visSenseRange();

void visLocalTarget(Point p);

void visRrt();

void visSamples( vector<Point> SampleSet );

bool checkHalfWay();

int trajGenerator(Eigen::MatrixXd path, Eigen::VectorXd radius, double time_pre);

Point getCommitedTarget();

rrtPathFinder _pathPlaner;

//TrajectoryGenerator     _trajectoryGenerator;
TrajectoryGeneratorSOCP _trajectoryGeneratorSocp;

void rcvWaypointsCallback(const nav_msgs::Path & wp)
{     
      if(wp.poses[0].pose.position.z < 0.0)
        return;

      endPt.x = wp.poses[0].pose.position.x;
      endPt.y = wp.poses[0].pose.position.y;
      endPt.z = wp.poses[0].pose.position.z;

      targetReceive = true;
      
      targetReach   = false;
      traFinish     = false;
      missionFinish = false;

      _waypoints = wp;
}

int trajGenerator(Eigen::MatrixXd path, Eigen::VectorXd radius, double time_pre)
{           
      Eigen::MatrixXd pos(2,3);
      Eigen::MatrixXd vel(2,3);
      Eigen::MatrixXd acc(2,3);
  
      if(traFinish)
      {
          double time_est_opt = 0.02;      
          double t_s = ( (_odom.header.stamp - init_time).toSec()  + time_est_opt + time_pre) ;                
          int idx;
          for (idx = 0; idx < _segment_num; ++idx){
              if (t_s > _Time(idx) && idx + 1 < _segment_num)
                  t_s -= _Time(idx);
              else break;
          }

          t_s /= _Time(idx);
          VectorXd state = getFullStateFromBezier(_PolyCoeff, t_s, idx);

          for(int i = 0; i < 3; i++ ){
              pos(0, i) = state(i) * _Time(idx);
              vel(0, i) = state(i + 3); 
              acc(0, i) = state(i + 6) / _Time(idx); 
          }  
      }
      else
      {   
          pos.row(0) << startPt.x,  startPt.y,  startPt.z;
          if(traHalt){
              ROS_WARN("[Planning Node] Navigation interupted, re-start from current states");
              vel.row(0) << 0.0, 0.0, 0.0;
              acc.row(0) << 0.0, 0.0, 0.0;        
              traHalt = false;
          }
          else{   
              ROS_WARN("[Planning Node] Navigation re-plan, get a new target");
              vel.row(0) << startVel.x, startVel.y, startVel.z;
              acc.row(0) << startAcc.x, startAcc.y, startAcc.z;
          }
      }

      pos.row(1) << endPt.x, endPt.y, endPt.z;
      vel.row(1) << 0.0, 0.0, 0.0;
      acc.row(1) << 0.0, 0.0, 0.0;
      
      /*Eigen::MatrixXd path_tmp;
      path_tmp.resize(path.rows() - 2, 3);
      path_tmp = path.block(1, 0, path.rows() - 2, 3);

      int i;
      for(i = path_tmp.rows() - 1; i >= 0 ; i-- )
      {
          if( sqrt( (path_tmp(i,0) - pos(0,0)) * (path_tmp(i,0) - pos(0,0)) + (path_tmp(i,1) - pos(0,1)) * (path_tmp(i,1) - pos(0,1)) 
            + (path_tmp(i,2) - pos(0,2)) * (path_tmp(i,2) - pos(0,2)) ) <  radius(i) )
            break;
      }
      
      Eigen::MatrixXd path_(path.rows() - i, 3);
      Eigen::VectorXd radius_(radius.size() - i);
      path_.row(0) = pos.row(0);
      path_.block(1, 0, path_tmp.rows() - i, 3) = path_tmp.block(i, 0, path_tmp.rows() - i, 3);
      radius_.segment(0, radius.size() - i)     = radius.segment(i, radius.size() - i); 
      path_.row(path_.rows() - 1) = pos.row(1);*/

      Eigen::MatrixXd path_   = path;
      Eigen::VectorXd radius_ = radius;
      _Time    = TimeAllocate( path_, radius_, vel.row(0) );
      _Scale   = _Time;
      _segment_num = radius_.size();

      //vector<int> _poly_orderList;
      // Now we assign a order to each piece od the trajectory according to its size
      // The priciple is , the smallest is the minimum order, and the largest one is the highest order

      auto min_r = radius.minCoeff();
      auto max_r = radius.maxCoeff();

      //ROS_BREAK();
      _poly_orderList.clear();

      //double d_order = (_poly_order_max - _poly_order_min ) / (max_r - min_r);
      double d_order = (_poly_order_max - _poly_order_min ) / (_max_radius - min_r);
      
      for(int i =0; i < _segment_num; i ++ )
      {  
          if(i == 0)
            _poly_orderList.push_back( _poly_order_max ); // 8
          else
            _poly_orderList.push_back( int( _poly_order_min + (radius(i) - min_r) * d_order ) );
        //_poly_orderList.push_back( 6 );
      }

  /*    for(int i =0; i < _segment_num; i ++ )
      {  
          _poly_orderList.push_back( _poly_order_max );
      }*/

      ROS_WARN("check order in each segment");
      for(auto p:_poly_orderList)
        cout<<p<<endl;

/*      cout<<"check the order in each piece"<<endl;
      for(auto ptr:_poly_orderList)
        cout<<ptr<<endl;*/
      //ROS_BREAK();
      /*_PolyCoeff = _trajectoryGenerator.BezierPloyCoeffGeneration(  
                  path_, radius_, _Time, _poly_orderList, _MQMList, pos, vel, acc, _MAX_Vel, _MAX_Acc, _minimize_order );*/
      
      _PolyCoeff = _trajectoryGeneratorSocp.BezierPloyCoeffGenerationSOCP(  
                  path_, radius_, _Time, _poly_orderList, _FMList, pos, vel, acc, _MAX_Vel, _MAX_Acc, _minimize_order );
      
      //traFinish     = false;
      //targetReceive = false;
      //return -1; 
      if(_PolyCoeff.rows() == 3 && _PolyCoeff.cols() == 3){
          ROS_WARN("Cannot find a feasible and optimal solution, somthing wrong with the mosek solver ... ");
          _traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_WARN_IMPOSSIBLE;
          _traj_pub.publish(_traj);
          traFinish     = false;
          targetReceive = false;
          return -1;
      }
      else{
          ROS_WARN("Trajectory generated successed");
          _traj = getBezierTraj();
          _traj_pub.publish(_traj);
          traFinish = true;
      }

      _traj_id ++;
      return 1;
}

inline double getDis(Point p1, Point p2){
      return sqrt( pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2) );
}

inline double getDis(Vector3d p1, Point p2){
      return sqrt( pow(p1(0) - p2.x, 2) + pow( p1(1) - p2.y, 2) + pow( p1(2) - p2.z, 2) );
}

void rcvPosCmdCallBack(const quadrotor_msgs::PositionCommand cmd)
{   
/*      _has_odom = true;
      _odom.pose.pose.position.x = cmd.position.x;
      _odom.pose.pose.position.y = cmd.position.y;
      _odom.pose.pose.position.z = cmd.position.z;

      _odom.twist.twist.linear.x = cmd.velocity.x;
      _odom.twist.twist.linear.y = cmd.velocity.y;
      _odom.twist.twist.linear.z = cmd.velocity.z;
      _odom.header = cmd.header;

      _odom.pose.pose.orientation.x = 0.0;
      _odom.pose.pose.orientation.y = 0.0;
      _odom.pose.pose.orientation.z = 0.0;
      _odom.pose.pose.orientation.w = 1.0;

      _state = {
          _odom.pose.pose.position.x, 
          _odom.pose.pose.position.y, 
          _odom.pose.pose.position.z, 
          _odom.twist.twist.linear.x,
          _odom.twist.twist.linear.y,
          _odom.twist.twist.linear.z,
          0.0, 0.0, 0.0
      };

      startPt.x  = _odom.pose.pose.position.x;
      startPt.y  = _odom.pose.pose.position.y;
      startPt.z  = _odom.pose.pose.position.z;

      startVel.x = _odom.twist.twist.linear.x;
      startVel.y = _odom.twist.twist.linear.y;
      startVel.z = _odom.twist.twist.linear.z;*/

      startAcc.x = cmd.acceleration.x;
      startAcc.y = cmd.acceleration.y;
      startAcc.z = cmd.acceleration.z;

/*      _pathPlaner.setStartPt(startPt, endPt);
      
      _odom_queue.push_back(_odom);
      while (_odom_queue.size() > _odom_queue_size) _odom_queue.pop_front();

      visSenseRange();
      visLocalTarget(commitTarget);

      if( !traFinish ) return;

      if(getDis(startPt, commitTarget) < _eps)
          targetReach = true;

      if(getDis(startPt, endPt) < _eps){
          missionFinish = true;
      }*/
}

void rcvOdometryCallbck(const nav_msgs::Odometry odom)
{
      //ROS_WARN("odometry callbck start");
      if (odom.child_frame_id == "X" || odom.child_frame_id == "O") return ;
      _odom = odom;
      _has_odom = true;

      _state = {
          _odom.pose.pose.position.x, 
          _odom.pose.pose.position.y, 
          _odom.pose.pose.position.z, 
          _odom.twist.twist.linear.x,
          _odom.twist.twist.linear.y,
          _odom.twist.twist.linear.z,
          0.0, 0.0, 0.0
      };

      startPt.x  = _odom.pose.pose.position.x;
      startPt.y  = _odom.pose.pose.position.y;
      startPt.z  = _odom.pose.pose.position.z;

      startVel.x = _odom.twist.twist.linear.x;
      startVel.y = _odom.twist.twist.linear.y;
      startVel.z = _odom.twist.twist.linear.z;

      _odom_queue.push_back(odom);
      while (_odom_queue.size() > _odom_queue_size) _odom_queue.pop_front();
      
      static tf::TransformBroadcaster br;
      tf::Transform transform;
      transform.setOrigin( tf::Vector3(_odom.pose.pose.position.x, _odom.pose.pose.position.y, _odom.pose.pose.position.z) );
      transform.setRotation(tf::Quaternion(0, 0, 0, 1.0));
      br.sendTransform(tf::StampedTransform(transform, _odom.header.stamp, "map", "quadrotor"));

      _pathPlaner.setStartPt(startPt, endPt);
      
      _odom_queue.push_back(_odom);
      while (_odom_queue.size() > _odom_queue_size) _odom_queue.pop_front();

      visSenseRange();
      visLocalTarget(commitTarget);

      if( !traFinish ) return;

      if(getDis(startPt, commitTarget) < _eps)
          targetReach = true;

      if(getDis(startPt, endPt) < _eps){
          missionFinish = true;
      }
}

bool checkEndOfCommitedTraj()
{     
      if(targetReach){     
        targetReach = false;
        return true;
      }
      else 
        return false;
}

double _r_eps = 0.2;
Point getCommitedTarget()
{
      Point commit_target;
      double t_s = _time_commit;      
  
      int idx;
      for (idx = 0; idx < _segment_num; ++idx){
          if (t_s > _Time(idx) && idx + 1 < _segment_num)
              t_s -= _Time(idx);
          else break;
      }          
      t_s /= _Time(idx);

      Vector3d ball_center = _Path.row(idx + 1);
      double   radius = _Radius(idx);
      /*cout<<"call center: \n"<<_Path.row(idx + 1)<<endl;
      cout<<"R: "<<_Radius(idx)<<endl;*/
      
      double d_t = 0.0;
      double t_start = t_s;
      double dir = 1.0;

      while(true)
      {      
        t_start += d_t; 

        Vector3d tentative_target = _Time(idx) * getStateFromBezier(_PolyCoeff, t_start, idx);
        cout<<"tentative_target: \n"<<tentative_target<<endl;

        double dis = (tentative_target - ball_center).norm();
        cout<<"distance to ball center is : "<<dis<<endl;
        if( dis < (radius - _r_eps) )
        {
            commit_target.x = tentative_target(0); 
            commit_target.y = tentative_target(1);
            commit_target.z = tentative_target(2);
            break; 
        }

        if(dis >= radius)
        {
          t_start = t_s;
          d_t = 0.0;
          dir = -1.0;
        }
        d_t += 0.01 * dir;
      }

      return commit_target;
}

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 & pointcloud_map )
{
      pcl::PointCloud<pcl::PointXYZ> CloudIn;      
      pcl::fromROSMsg(pointcloud_map, CloudIn);
      _pathPlaner.setInput(CloudIn);
}

void rcvAddedPointCloudCallBack(const sensor_msgs::PointCloud2 & add_pointcloud_map )
{     
      pcl::PointCloud<pcl::PointXYZ> CloudAdd;        
      pcl::fromROSMsg( add_pointcloud_map, CloudAdd );
      if(CloudAdd.points.size() == 0) return;
      _pathPlaner.rcvAddMap( CloudAdd );
}

void rcvDeletedPointCloudCallBack(const sensor_msgs::PointCloud2 & delete_pointcloud_map )
{     
      pcl::PointCloud<pcl::PointXYZ> CloudDel;      
      pcl::fromROSMsg( delete_pointcloud_map, CloudDel );
      if(CloudDel.points.size() == 0) return;
      _pathPlaner.rcvDelMap( CloudDel );
}

void InitPlanning()
{
      _pathPlaner.reset();
      _pathPlaner.setPt( startPt, endPt, _x_l, _x_h, _y_l, _y_h, _z_l, _z_h, bias_l, bias_h, _sensing_range, _max_samples, _sample_portion, _goal_portion );

      ros::Time timeBef = ros::Time::now();
      _pathPlaner.RRTpathFind(_path_find_limit);
      ros::Time timeAft = ros::Time::now();
      double _path_time = (timeAft-timeBef).toSec();
  
      pair<Eigen::MatrixXd, Eigen::VectorXd> res = _pathPlaner.getPath();
      if(_pathPlaner.path_find_state == false)
      {
          ROS_WARN("[Planning Node] Can't find a path, mission stall, please reset the target");
          _traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_WARN_IMPOSSIBLE;
          _traj_pub.publish(_traj);
          traFinish = false;
          targetReceive = false;
      }
      else  // Path finding succeed .
      {
          _Path   = res.first;
          _Radius = res.second;
          
          if(trajGenerator( _Path, _Radius, _path_time ) == 1) 
          {   
              ROS_WARN("[Planning Node] trajectory generation succeed .");
              commitTarget = getCommitedTarget();
              _pathPlaner.resetRoot(commitTarget);

              visBezierTrajectory(_PolyCoeff);
              visCtrlPoint( _PolyCoeff );
              visPath(_Path, _Radius);
          }
      }

      vector<Point> sampleSet = _pathPlaner.getSamples();
      nodesAll = _pathPlaner.getTree();
      visRrt();
      visSamples(sampleSet);
}

void IncrePlanning()
{     
      if(_pathPlaner.commitEnd == true){
        nodesAll = _pathPlaner.getTree();
        visRrt();
        return;
      }

      if( checkEndOfCommitedTraj() ) // almost reach the end of the commited trajectory.
      { 
          ROS_WARN("[Planning Node] almost reach the commited target, start generating the next trajectory"); 
          //pair<Eigen::MatrixXd, Eigen::VectorXd> res = _pathPlaner.getPath(); // find the path in tht tree
          if(_pathPlaner.path_find_state == false) // no feasible path exists
          {
              ROS_WARN("[Planning Node] reach commited target but no path exists, waiting for a path");
              
              _traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_WARN_IMPOSSIBLE;
              _traj_pub.publish(_traj);
              traFinish = false;  // continue to fly
              traHalt   = true;  // wait to re-generate from 0 states
              //targetReceive = false; // try to update this logic, when can't find a path, continue expand the tree, stop and wait
              return;
          }
          else
          {   
              //targetReach = false;
              visPath(_Path, _Radius);
              _Radius(0) += _search_margin; // TEST: relax the first sphere's constraint, hope to improve robustness
              if(trajGenerator( _Path, _Radius, 0.15 ) == 1)      // Generate a new trajectory.
              {
                commitTarget = getCommitedTarget();   
                _pathPlaner.resetRoot(commitTarget);  
                visBezierTrajectory(_PolyCoeff);
                visCtrlPoint( _PolyCoeff );
              }

          }
      }
      else
      { // continue to refine the uncommited trajectory
          //if(_pathPlaner.refine_status == true)
          _pathPlaner.RRTpathRefine( refine_limit, false );    // add samples to the tree
          ros::Time time_bef_evaluate = ros::Time::now();
          _pathPlaner.RRTpathReEvaluate(revaluate_limit, true ); // ensure that the path is collision-free
          ros::Time time_aft_evaluate = ros::Time::now();
          ROS_WARN("[Planning Node] Time in re-evaluate the path is %f", (time_aft_evaluate - time_bef_evaluate).toSec());

          if(_pathPlaner.path_find_state == true)
          {
              pair<Eigen::MatrixXd, Eigen::VectorXd> res = _pathPlaner.getPath();
              _Path   = res.first;
              _Radius = res.second;
              visPath(_Path, _Radius);
          }

          nodesAll = _pathPlaner.getTree();
          visRrt(); 
      }

      /*  Debug use only  */
/*      _pathPlaner.RRTpathRefine( refine_limit, false );
      nodesAll = _pathPlaner.getTree();
      visRrt(); 
      pair<Eigen::MatrixXd, Eigen::VectorXd> res = _pathPlaner.getPath();
      _Path   = res.first;
      _Radius = res.second;
      visPath(_Path, _Radius);*/
}

void Task()
{
    if(targetReceive == false)
      return;
    
    if(traFinish == false)
      InitPlanning();

      IncrePlanning();  
}

int main (int argc, char** argv) 
{        
      ros::init (argc, argv, "pcd_RRT");
      ros::NodeHandle n( "~" );

      n.param("sensing_rate",    _sensing_rate,  10.0f);      
      n.param("PlanParam/safety_margin",    _safety_margin,  0.65f);
      n.param("PlanParam/search_margin",    _search_margin,  0.35f);
      n.param("PlanParam/max_radius",       _max_radius,     10.0f);
      n.param("PlanParam/sensing_range",    _sensing_range,  10.0f);     
      n.param("PlanParam/refine_portion",   _refine_portion,  0.8f);     
      n.param("PlanParam/sample_portion",   _sample_portion,  0.25f);     // the ratio to generate samples inside the map range
      n.param("PlanParam/goal_portion",     _goal_portion,    0.05f);     // the ratio to generate samples on the goal
      n.param("PlanParam/path_find_limit",  _path_find_limit, 0.05f);     
      n.param("PlanParam/max_samples",      _max_samples,     3000);     
      n.param("dynamic/vec",   _Vel,  2.0);
      n.param("dynamic/acc",   _Acc,  1.0);
      n.param("dynamic/max_vec",   _MAX_Vel,  3.0);
      n.param("dynamic/max_acc",   _MAX_Acc,  1.5);
      n.param("mapBoundary/lower_x", _x_l,  -50.0f);
      n.param("mapBoundary/upper_x", _x_h,   50.0f);
      n.param("mapBoundary/lower_y", _y_l,  -50.0f);
      n.param("mapBoundary/upper_y", _y_h,   50.0f);
      n.param("mapBoundary/lower_z", _z_l,    0.0f);
      n.param("mapBoundary/upper_z", _z_h,    3.0f);
      n.param("optimization/poly_order_min", _poly_order_min,  5);
      n.param("optimization/poly_order_max", _poly_order_max,  10);
      n.param("optimization/minimize_order", _minimize_order,  3);
      n.param("commitTime", _time_commit,  1.0);

      Bernstein _bernstein;
      if(_bernstein.setParam(_poly_order_min, _poly_order_max, _minimize_order) == -1)
      {
          ROS_ERROR(" The trajectory order is set beyond the library's scope, please re-set ");
      }

      _MQMList = _bernstein.getMQM();
      _FMList   = _bernstein.getFM();
      _CList   = _bernstein.getC();
      _CvList  = _bernstein.getC_v();
      _CaList  = _bernstein.getC_a();

      _eps = 0.25; 

      _pathPlaner.setParam(_safety_margin, _search_margin, _max_radius, _sensing_range);

      time_limit      = 1.0f / _sensing_rate;
      refine_limit    =  _refine_portion * time_limit;
      revaluate_limit = (1 - _refine_portion) * time_limit; 
      cout<<"refine_limit: "<<refine_limit<<endl;
            
      // subcribed msgs
      _dest_pts_sub = n.subscribe( "waypoints", 1, rcvWaypointsCallback );
      _map_sub     = n.subscribe( "PointCloud", 1, rcvPointCloudCallBack);
      _add_map_sub = n.subscribe( "DeltaPointCloud", 1, rcvAddedPointCloudCallBack);
      //_del_map_sub = n.subscribe( "DeltaPointCloud", 1, rcvDeletedPointCloudCallBack);
      _cmd_sub    = n.subscribe( "position_cmd",50, rcvPosCmdCallBack);
      _odom_sub   = n.subscribe( "odometry", 50, rcvOdometryCallbck );

      // publish visualize related msgs
      _ctrlPt_vis_pub       = n.advertise<visualization_msgs::MarkerArray>("ctrlPt", 1);                      
      _locTarget_vis_pub    = n.advertise<visualization_msgs::MarkerArray>("commitedTarget", 1);                      
      _traj_vis_pub         = n.advertise<visualization_msgs::Marker>("trajectory_vis", 1);
      _checkTraj_vis_pub    = n.advertise<visualization_msgs::Marker>("check_trajectory_vis", 1);
      _path_vis_pub         = n.advertise<visualization_msgs::MarkerArray>("path", 1);
      _traj_pub             = n.advertise<quadrotor_msgs::PolynomialTrajectory>("trajectory", 10);
      _sensing_range_pub    = n.advertise<visualization_msgs::Marker>("sensing_range", 1);
      _vis_rrt_star_pub     = n.advertise<visualization_msgs::Marker>("tree", 1);
      _vis_samples_pub      = n.advertise<visualization_msgs::Marker>("samples", 1);
     
      ros::Rate rate(10);
      bool status = ros::ok();
      while (status) {
          ros::spinOnce();

          Task();

          status = ros::ok();
          rate.sleep();
      }
}

quadrotor_msgs::PolynomialTrajectory getBezierTraj()
{
      quadrotor_msgs::PolynomialTrajectory traj;
      traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;
      traj.num_segment = _Time.size();

//      ROS_WARN("Re-stack the coefficients");

      int polyTotalNum = 0;
      for(auto order:_poly_orderList)
          polyTotalNum += (order + 1);

      traj.coef_x.resize(polyTotalNum);
      traj.coef_y.resize(polyTotalNum);
      traj.coef_z.resize(polyTotalNum);

/*      cout<<"polyTotalNum: "<<polyTotalNum<<endl;
*/
      /*cout<<"_poly_orderList"<<endl;
      for(auto it:_poly_orderList)
        cout<<it<<endl;*/

      //cout<<"_PolyCoeff: \n"<<_PolyCoeff<<endl;
      int idx = 0;
      for(int i = 0; i < _segment_num; i++ )
      {    
          int order = _poly_orderList[i];
          int poly_num1d = order + 1;
          
          for(int j =0; j < poly_num1d; j++)
          { 
        /*      cout<<"idx:"<< idx<<endl;
              cout<<"x coeff: "<<_PolyCoeff(i,                  j)<<endl;
              cout<<"y coeff: "<<_PolyCoeff(i,     poly_num1d + j)<<endl;
              cout<<"z coeff: "<<_PolyCoeff(i, 2 * poly_num1d + j)<<endl;*/

              traj.coef_x[idx] = _PolyCoeff(i,                  j);
              traj.coef_y[idx] = _PolyCoeff(i,     poly_num1d + j);
              traj.coef_z[idx] = _PolyCoeff(i, 2 * poly_num1d + j);
              idx++;
          }
      }


//      ROS_WARN("Fill the coefficients into the msgs");
      traj.header.frame_id = "/bernstein";
      traj.header.stamp = _odom.header.stamp; //ros::Time(_odom.header.stamp.toSec()); 
      init_time = traj.header.stamp;

      traj.time.resize(_Time.size());
      traj.order.resize(_Time.size());

      _start_time = traj.header.stamp; 
      traj.mag_coeff = 1.0;

      for (int idx = 0; idx < _Time.size(); ++idx){
          traj.time[idx] = _Time(idx);
          traj.order[idx] = _poly_orderList[idx];
      }
      
      traj.start_yaw = 0.0;//tf::getYaw(_odom.pose.pose.orientation);
      traj.final_yaw = 0.0;//tf::getYaw(_waypoints.poses.back().pose.orientation);

      if (abs( (traj.start_yaw + 2 * _PI) - traj.final_yaw) < _PI)
          traj.start_yaw += 2 * _PI;
      if (abs( (traj.start_yaw - 2 * _PI) - traj.final_yaw) < _PI)
          traj.start_yaw -= 2 * _PI;

      traj.trajectory_id = _traj_id;
      traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;

      //ROS_BREAK();
      return traj;
}

Eigen::VectorXd TimeAllocate( Eigen::MatrixXd Path, Eigen::VectorXd Radius, Eigen::Vector3d iniVel )
{
      int ball_num = Path.rows() - 2;

      Eigen::MatrixXd check_pt( ball_num - 1 , 3 );
      Eigen::MatrixXd Path_ball = Path.block( 1, 0, ball_num, 3 );

      const static auto getdist = []( const Eigen::Vector3d u, const Eigen::Vector3d v ){
            const static auto sq = [] (double && val){ 
                  return val * val;
              };

            return sqrt( sq(v[0] - u[0]) + sq(v[1] - u[1]) + sq(v[2] - u[2])  );
      };

      Eigen::Vector3d center, last_center;
      double radius, last_radius;

      int index = 0 ;
      last_center <<  Path_ball(index, 0), Path_ball(index, 1), Path_ball(index, 2); 
      last_radius  =  Radius[index];

      for(index = 1; index < ball_num; index ++ ){   
          center <<  Path_ball(index, 0), Path_ball(index, 1), Path_ball(index, 2); 
          radius  =  Radius[index];

          double dist = getdist(last_center, center);  
          
          Eigen::Vector3d delta_Vec = center - last_center;
          Eigen::Vector3d joint_pt;
          joint_pt = delta_Vec * ( dist + last_radius - radius) / ( 2.0 * dist ) + last_center; 

          check_pt.block( index - 1, 0, 1, 3 ) = joint_pt.transpose();

          last_center = center;
          last_radius = radius;
      }

      Eigen::MatrixXd all_points( ball_num + 1 , 3 );
      all_points.row( 0 )                     = Path.row( 0 );
      all_points.block( 1, 0, ball_num-1, 3 ) = check_pt;
      all_points.row( ball_num )              = Path.row( ball_num + 1 );

      Eigen::VectorXd time_allocate(all_points.rows() - 1);

      Eigen::Vector3d initv = iniVel;
      for (int k = 0; k < all_points.rows() - 1; k++){
          // Position time
          double dtxyz;

          Eigen::Vector3d p0   = all_points.row(k);           // The start point of this segment
          Eigen::Vector3d p1   = all_points.row(k + 1);       // The end point of this segment
          Eigen::Vector3d d    = p1 - p0;                     // The position difference
          Eigen::Vector3d v0(0.0, 0.0, 0.0);                  // The init velocity
          
          if( k == 0) v0 = initv;

          double D    = d.norm();                             // The norm of position difference 
          double V0   = v0.dot(d / D);                        // Map velocity to the vector (p1-p0)
          double aV0  = fabs(V0);                             // The absolote mapped velocity

          double acct = (_Vel - V0) / _Acc * ((_Vel > V0)?1:-1);                      // The time to speed up to to the max veloctiy
          double accd = V0 * acct + (_Acc * acct * acct / 2) * ((_Vel > V0)?1:-1);    // The distance to speed up
          double dcct = _Vel / _Acc;                                                  // The time to slow down.
          double dccd = _Acc * dcct * dcct / 2;                                       // The distance to slow down.

          if (D < aV0 * aV0 / (2 * _Acc)){                 // if not enough distance to slow down
            double t1 = (V0 < 0)?2.0 * aV0 / _Acc:0.0;
            double t2 = aV0 / _Acc;
            dtxyz     = t1 + t2;                 
          }
          else if (D < accd + dccd){                      // if not enough distance to speed up and slow dwon 
            double t1 = (V0 < 0)?2.0 * aV0 / _Acc:0.0;
            double t2 = (-aV0 + sqrt(aV0 * aV0 + _Acc * D - aV0 * aV0 / 2)) / _Acc;
            double t3 = (aV0 + _Acc * t2) / _Acc;
            dtxyz     = t1 + t2 + t3;    
          }
          else{                                          // can reach max veloctiy and keep for a while.
            double t1 = acct;                              
            double t2 = (D - accd - dccd) / _Vel;
            double t3 = dcct;
            dtxyz     = t1 + t2 + t3;                                                                  
          }
          time_allocate(k) = dtxyz;

          /*if(k == 0)
            time_allocate(k) *= 2.0;*/
          
          //time_allocate(k) = D/_Vel;
          //time_allocate(k) = 1.0;
          //cout<<"ratio: "<<D/dtxyz<<endl;
      }

      return time_allocate;
}

VectorXd getFullStateFromBezier(const Eigen::MatrixXd & polyCoeff,  double t_now, int seg_now )
{
      VectorXd ret = VectorXd::Zero(9);
      Eigen::VectorXd ctrl_now = polyCoeff.row(seg_now);

      int order = _poly_orderList[seg_now];
      int ctrl_num1D = order + 1;

      for(int i = 0; i < 3; i++)
      {   
          for(int j = 0; j < ctrl_num1D; j++){
              ret(i) += _CList[order](j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (order - j) ); 
              
              if(j < ctrl_num1D - 1 )
                ret(i+3) += _CvList[order](j) * order 
                          * ( ctrl_now(i * ctrl_num1D + j + 1) - ctrl_now(i * ctrl_num1D + j))
                          * pow(t_now, j) * pow((1 - t_now), (order - j - 1) ); 
              
              if(j < ctrl_num1D - 2 )
                ret(i+6) += _CaList[order](j) * order * (order - 1) 
                          * ( ctrl_now(i * ctrl_num1D + j + 2) - 2 * ctrl_now(i * ctrl_num1D + j + 1) + ctrl_now(i * ctrl_num1D + j))
                          * pow(t_now, j) * pow((1 - t_now), (order - j - 2) );                         
          }

      }
      return ret;  
}

Vector3d getStateFromBezier(const Eigen::MatrixXd & polyCoeff,  double t_now, int seg_now )
{
      Vector3d ret = VectorXd::Zero(3);
      Eigen::VectorXd ctrl_now = polyCoeff.row(seg_now);

      int order = _poly_orderList[seg_now];
      int ctrl_num1D = order + 1;
      for(int i = 0; i < 3; i++)
      {   
          for(int j = 0; j < ctrl_num1D; j++){
              ret(i) += _CList[order](j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (order - j) ); 
          }
      }

      //ROS_WARN("finish get state from Bezier curve ...");
      return ret;  
}

bool checkHalfWay()
{   
    if(!traFinish) return false;

    Point check_traj_pt;
    VectorXd state;

    visualization_msgs::Marker _traj_vis;

    geometry_msgs::Point pt;
    _traj_vis.header.stamp       = ros::Time::now();
    _traj_vis.header.frame_id    = "map";

    _traj_vis.ns = "trajectory/trajectory";
    _traj_vis.id = 0;
    _traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    _traj_vis.action = visualization_msgs::Marker::ADD;
    _traj_vis.scale.x = 1.2 * _vis_traj_width;
    _traj_vis.scale.y = 1.2 * _vis_traj_width;
    _traj_vis.scale.z = 1.2 * _vis_traj_width;
    _traj_vis.pose.orientation.x = 0.0;
    _traj_vis.pose.orientation.y = 0.0;
    _traj_vis.pose.orientation.z = 0.0;
    _traj_vis.pose.orientation.w = 1.0;
    _traj_vis.color.r = 1.0;
    _traj_vis.color.g = 0.0;
    _traj_vis.color.b = 0.0;
    _traj_vis.color.a = 1.0;

    double t_s = max(0.0, (_odom.header.stamp - _start_time).toSec());      
    int idx;
    for (idx = 0; idx < _segment_num; ++idx){
      if (t_s > _Time(idx) && idx + 1 < _segment_num)
          t_s -= _Time(idx);
      else break;
    }
    
    double t_ss;
    for(int i = idx; i < _segment_num; i++ ){
      t_ss = (i == idx) ? t_s : 0.0;
      for (double t = t_ss; t < _Time(i); t += 0.01){
        state = getFullStateFromBezier(_PolyCoeff, t/_Time(i), i);
        check_traj_pt.x = _Time(i) * state(0); 
        check_traj_pt.y = _Time(i) * state(1);
        check_traj_pt.z = _Time(i) * state(2);
        
        pt.x = _Time(i) * state(0);
        pt.y = _Time(i) * state(1);
        pt.z = _Time(i) * state(2);
        
        _traj_vis.points.push_back(pt);
        if(_pathPlaner.checkTrajInvalid(check_traj_pt)){
          ROS_WARN("[Planning Node] Collision Occur, Need Re-Plan");
          _checkTraj_vis_pub.publish(_traj_vis);
          return true;
        }
      }
    }
    _checkTraj_vis_pub.publish(_traj_vis); 

    return false;
}

void visSenseRange()
{    
      visualization_msgs::Marker mk;
      mk.header.frame_id = "map";
      mk.header.stamp = ros::Time::now();
      mk.ns = "Pcdplanner/sensing_range";
      mk.type = visualization_msgs::Marker::CYLINDER;
      mk.action = visualization_msgs::Marker::ADD;
      mk.pose.orientation.x = 0.0;
      mk.pose.orientation.y = 0.0;
      mk.pose.orientation.z = 0.0;
      mk.pose.orientation.w = 1.0;
      mk.color.a = 0.4;
      mk.color.r = 0.0;
      mk.color.g = 1.0;
      mk.color.b = 0.0;
      mk.scale.x = 2 * _sensing_range;
      mk.scale.y = 2 * _sensing_range;
      mk.scale.z = 2 * _sensing_range * sin(15.0 / 180.0 * M_PI);

      mk.id = 0;
      mk.pose.position.x = startPt.x; 
      mk.pose.position.y = startPt.y; 
      mk.pose.position.z = startPt.z; 
          
      _sensing_range_pub.publish(mk);
}

void visRrt()
{     
    //ROS_WARN("Prepare to vis rrt result");
    visualization_msgs::Marker points, line_list, root;

    root.header.frame_id    = points.header.frame_id    = line_list.header.frame_id    = "/map";
    root.header.stamp       = points.header.stamp       = line_list.header.stamp       = ros::Time::now();
    root.ns                 = "/root";
    points.ns               = "/vertex";
    line_list.ns            = "/edges";
    root.action             = points.action             = line_list.action             = visualization_msgs::Marker::ADD;
    root.pose.orientation.w = points.pose.orientation.w = line_list.pose.orientation.w = 1.0;
    root.pose.orientation.x = points.pose.orientation.x = line_list.pose.orientation.x = 0.0;
    root.pose.orientation.y = points.pose.orientation.y = line_list.pose.orientation.y = 0.0;
    root.pose.orientation.z = points.pose.orientation.z = line_list.pose.orientation.z = 0.0;
    root.id                 = points.id                 = line_list.id = 0;
    root.type               = points.type               = visualization_msgs::Marker::SPHERE_LIST;
    line_list.type = visualization_msgs::Marker::LINE_LIST;
    
    root.scale.x      = root.scale.y      = root.scale.z = 0.25;
    points.scale.x    = points.scale.y    = points.scale.z = 0.15;
    line_list.scale.x = line_list.scale.y = line_list.scale.z = 0.025;

    root.color.a = 1.0;
    root.color.r = 0.0;
    root.color.g = 1.0;
    root.color.b = 1.0;
    root.points.clear();

    points.color.a = 1.0;
    points.color.r = 1.0;
    points.color.g = 0.0;
    points.color.b = 0.0;
    points.points.clear();

    line_list.color.a = 0.7;
    line_list.color.r = 0.0;
    line_list.color.g = 1.0;
    line_list.color.b = 0.0;
    line_list.points.clear();

    NodePtr nodeptr;
    for(int i = 0; i < int(nodesAll.size()); i++){
        nodeptr = nodesAll[i];
        if (nodeptr->preNode_ptr == NULL){ 
            geometry_msgs::Point r;
            r.x = nodeptr->coord.x;
            r.y = nodeptr->coord.y; 
            r.z = nodeptr->coord.z; 
            root.points.push_back(r);
            continue;
        }

        geometry_msgs::Point p;
        p.x = nodeptr->coord.x;
        p.y = nodeptr->coord.y; 
        p.z = nodeptr->coord.z; 
        points.points.push_back(p);

        geometry_msgs::Point p_line;
        p_line = p;
        line_list.points.push_back(p_line);
        
        p_line.x = nodeptr->preNode_ptr->coord.x;
        p_line.y = nodeptr->preNode_ptr->coord.y; 
        p_line.z = nodeptr->preNode_ptr->coord.z; 
        line_list.points.push_back(p_line);
    }

    _vis_rrt_star_pub.publish(root);
    _vis_rrt_star_pub.publish(points);
    _vis_rrt_star_pub.publish(line_list);
}

void visSamples( vector<Point> SampleSet )
{     
    visualization_msgs::Marker points;

    points.header.frame_id    = "/map";
    points.header.stamp       = ros::Time::now();
    points.ns                 = "/vertex";
    points.action             = visualization_msgs::Marker::ADD;
    points.pose.orientation.w = 1.0;
    points.pose.orientation.x = 0.0;
    points.pose.orientation.y = 0.0;
    points.pose.orientation.z = 0.0;
    points.id                 = 0;
    points.type    = visualization_msgs::Marker::SPHERE_LIST;
    
    points.scale.x = 0.5;
    points.scale.y = 0.5;
    points.scale.z = 0.5;
    points.color.a = 1.0;
    points.color.r = 0.0;
    points.color.g = 1.0;
    points.color.b = 0.0;

    geometry_msgs::Point p;
    for(auto ptr : SampleSet){
        p.x = ptr.x;
        p.y = ptr.y; 
        p.z = ptr.z; 
        points.points.push_back(p);
    }

    _vis_samples_pub.publish(points);
}

void visBezierTrajectory(Eigen::MatrixXd polyCoeff)
{
    visualization_msgs::Marker _traj_vis;

    _traj_vis.header.stamp       = ros::Time::now();
    _traj_vis.header.frame_id    = "map";

    _traj_vis.ns = "trajectory/trajectory";
    _traj_vis.id = 0;
    _traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    _traj_vis.action = visualization_msgs::Marker::ADD;
    _traj_vis.scale.x = _vis_traj_width;
    _traj_vis.scale.y = _vis_traj_width;
    _traj_vis.scale.z = _vis_traj_width;
    _traj_vis.pose.orientation.x = 0.0;
    _traj_vis.pose.orientation.y = 0.0;
    _traj_vis.pose.orientation.z = 0.0;
    _traj_vis.pose.orientation.w = 1.0;
    _traj_vis.color.r = 0.0;
    _traj_vis.color.g = 1.0;
    _traj_vis.color.b = 0.0;
    _traj_vis.color.a = 0.3;

    double traj_len = 0.0;
    int count = 0;
    Eigen::Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    _traj_vis.points.clear();

    Vector3d coord;
    geometry_msgs::Point pt;

    for(int i = 0; i < _segment_num; i++ ){
        for (double t = 0.0; t < 1.0; t += 0.01 / _Scale(i), count += 1){
            coord = getStateFromBezier( polyCoeff, t, i );
            cur(0) = pt.x = _Scale(i) * coord(0);
            cur(1) = pt.y = _Scale(i) * coord(1);
            cur(2) = pt.z = _Scale(i) * coord(2);
            _traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);

    _traj_vis_pub.publish(_traj_vis);
}

void visPath(Eigen::MatrixXd path, Eigen::VectorXd radius)
{     
      
      Eigen::MatrixXd path_new  = path.block(1, 0, path.rows() - 2, path.cols() );

      for (auto & mk: path_vis.markers) 
        mk.action = visualization_msgs::Marker::DELETE;

      _path_vis_pub.publish(path_vis);
      path_vis.markers.clear();

      visualization_msgs::Marker mk;
      mk.header.frame_id = "map";
      mk.header.stamp = ros::Time::now();
      mk.ns = "pcd_RRT/path";
      mk.type = visualization_msgs::Marker::SPHERE;
      mk.action = visualization_msgs::Marker::ADD;
      mk.pose.orientation.x = 0.0;
      mk.pose.orientation.y = 0.0;
      mk.pose.orientation.z = 0.0;
      mk.pose.orientation.w = 1.0;
      mk.color.a = 0.4;
      mk.color.r = 1.0;
      mk.color.g = 1.0;
      mk.color.b = 1.0;

      for(int i = 0; i < int(path_new.rows()); i++){
          mk.id = i;
          mk.pose.position.x = path_new(i, 0); 
          mk.pose.position.y = path_new(i, 1); 
          mk.pose.position.z = path_new(i, 2); 
          mk.scale.x = 2 * radius(i);
          mk.scale.y = 2 * radius(i);
          mk.scale.z = 2 * radius(i);
          
          path_vis.markers.push_back(mk);
    }

    _path_vis_pub.publish(path_vis);
}

void visCtrlPoint(Eigen::MatrixXd polyCoeff)
{
      for (auto & mk: ctrlPt_vis.markers) 
          mk.action = visualization_msgs::Marker::DELETE;
      
      _ctrlPt_vis_pub.publish(ctrlPt_vis);

      ctrlPt_vis.markers.clear();
      visualization_msgs::Marker mk;
      mk.header.frame_id = "map";
      mk.header.stamp = ros::Time::now();
      mk.ns = "pcd_RRT/ctrlPt";
      mk.type = visualization_msgs::Marker::SPHERE;
      mk.action = visualization_msgs::Marker::ADD;
      mk.pose.orientation.x = 0.0;
      mk.pose.orientation.y = 0.0;
      mk.pose.orientation.z = 0.0;
      mk.pose.orientation.w = 1.0;
      mk.color.a = 1.0;
      mk.color.r = 1.0;
      mk.color.g = 0.0;
      mk.color.b = 0.0;

      int idx = 0;
      for(int i = 0; i < _segment_num; i++)
      {   
          int order = _poly_orderList[i];
          int ctrl_num = order + 1;
          
          for(int j = 0; j < ctrl_num; j++)
          {
              mk.id = idx;
              mk.pose.position.x = _Scale(i) * polyCoeff(i, j);
              mk.pose.position.y = _Scale(i) * polyCoeff(i, ctrl_num + j);
              mk.pose.position.z = _Scale(i) * polyCoeff(i, 2 * ctrl_num + j);
              mk.scale.x = 0.25;
              mk.scale.y = 0.25;
              mk.scale.z = 0.25;
              ctrlPt_vis.markers.push_back(mk);
              idx ++;
          }
    }

    _ctrlPt_vis_pub.publish(ctrlPt_vis);
    //ROS_WARN("vis control points OK");
}

void visLocalTarget(Point p)
{
      if(!traFinish)
        return;

      for (auto & mk: locTarget_vis.markers) 
          mk.action = visualization_msgs::Marker::DELETE;
      
      locTarget_vis.markers.clear();
      _locTarget_vis_pub.publish(locTarget_vis);

      visualization_msgs::Marker mk;
      mk.header.frame_id = "map";
      mk.header.stamp = ros::Time::now();
      mk.ns = "pcd_RRT/lcoalTarget";
      mk.type = visualization_msgs::Marker::SPHERE;
      mk.action = visualization_msgs::Marker::ADD;
      mk.pose.orientation.x = 0.0;
      mk.pose.orientation.y = 0.0;
      mk.pose.orientation.z = 0.0;
      mk.pose.orientation.w = 1.0;
      mk.color.a = 1.0;
      mk.color.r = 0.0;
      mk.color.g = 0.0;
      mk.color.b = 1.0; 
      mk.scale.x = 0.25;
      mk.scale.y = 0.25;
      mk.scale.z = 0.25;

      mk.pose.position.x  = p.x;
      mk.pose.position.y  = p.y;
      mk.pose.position.z  = p.z;
      locTarget_vis.markers.push_back(mk);
      _locTarget_vis_pub.publish(locTarget_vis);

      //ROS_WARN("[Planning Node] Finish vis the local target");
  }