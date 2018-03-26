//#define PCL_NO_PRECOMPILE
#include <iostream>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>

#include <pcl/kdtree/kdtree_flann.h>
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

#include "pcd_trajectory/time_optimizer.h"
#include "pcd_trajectory/trajectory_generator.h"
#include "pcd_trajectory/trajectory_generator_discrete.h"
#include "pcd_trajectory/trajectory_generator_lite.h"
#include "pcd_trajectory/trajectory_generator_socp_lite.h"
#include "pcd_trajectory/rrgPathFinder.h"
#include "pcd_trajectory/backward.hpp"
#include "pcd_trajectory/bezier_base.h"
#include "quadrotor_msgs/PositionCommand.h"
#include "quadrotor_msgs/PolynomialTrajectory.h"

namespace backward {
backward::SignalHandling sh;
}

using namespace std;
using namespace Eigen;
using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;

visualization_msgs::MarkerArray path_vis;

ros::Publisher _traj_vis_pub;
ros::Publisher _traj_bezier_vis_pub;
ros::Publisher _checkPt_vis_pub;
ros::Publisher _poly_traj_vis_pub;
ros::Publisher _dis_traj_vis_pub;
ros::Publisher _checkTraj_vis_pub;
ros::Publisher _path_vis_pub;
ros::Publisher _ctrlPt_vis_pub;
ros::Publisher _disPath_vis_pub;
ros::Publisher _traj_pub;

ros::Subscriber _dest_pts_sub;
ros::Subscriber _map_sub;

nav_msgs::Path _waypoints;

float _vis_traj_width = 0.15; 

double _MAX_Vel,_MAX_Acc, _Vel, _Acc, _eps;
double _l;

int  _minimize_order, _poly_order_min, _poly_order_max;
int  _segment_num;
vector<int> _poly_orderList;

MatrixXd _PolyCoeff_p;
MatrixXd _PolyCoeff_b;
MatrixXd _PolyCoeff_b_socp;

VectorXd _Time;
VectorXd _Radius;
MatrixXd _Path;

Point startPt  = { 0.0, 0.0, 0.0}; // For outdoor CYT model .
Point startVel = { 0.0, 0.0, 0.0};
Point startAcc = { 0.0, 0.0, 0.0};
Point endPt;      

vector<MatrixXd> _MQMList, _FMList;
vector<VectorXd> _CList, _CvList, _CaList;

vector<double> getStateFromBezier(const MatrixXd & polyCoeff,  double t_now, int seg_now );
vector<double> getFullStateFromBezier(const MatrixXd & polyCoeff,  double t_now, int seg_now );

vector<double> getStateFromPolynomial( const MatrixXd & PolyCoeff, const VectorXd & Time, double t_now, int seg_now );
vector<double> getFullStateFromPolynomial( const MatrixXd & PolyCoeff, const VectorXd & Time, double t_now, int seg_now );
VectorXd TimeAllocate(MatrixXd Path, VectorXd Radius, Vector3d iniVel );
void trajGenerationTest(MatrixXd path, VectorXd radius);

void visBezierTrajectory(MatrixXd polyCoeff, int rgb);
void visCtrlPoint(MatrixXd polyCoeff, int rgb);
void visPath(MatrixXd path, VectorXd radius);
void visPolyTraj(vector<MatrixXd> polyCoeffList);
void visCheckPt(vector<MatrixXd> polyCoeffList);
void visDiscretePath(vector<Vector3d> path_dis);
void visDisTraj(MatrixXd disAcc, double h);

double x_l = -40.0, x_h = 40.0, y_l = -40.0, y_h = 40.0, z_l = 0.0, z_h = 5.0, bias_l = 0.0, bias_h = 1.0; 
rrgPathFinder _pathPlaner(x_l, x_h, y_l, y_h, z_l, z_h, bias_l, bias_h);

TrajectoryGeneratorDiscrete  _trajectoryGeneratorDiscrete;
TrajectoryGeneratorPoly      _trajectoryGeneratorPoly;
TrajectoryGeneratorBezier    _trajectoryGeneratorBezier;
TrajectoryGeneratorSOCP      _trajectoryGeneratorSocp;

void rcvWaypointsCallback(const nav_msgs::Path & wp)
{     
      if(wp.poses[0].pose.position.z < 0.0)
        return;

      endPt.x = wp.poses[0].pose.position.x;
      endPt.y = wp.poses[0].pose.position.y;
      endPt.z = wp.poses[0].pose.position.z;

      _pathPlaner.reset();
      _pathPlaner.setPt(startPt, endPt);
      _pathPlaner.RRGpathFind();

      if(!_pathPlaner.path_find_state){
        ROS_WARN("In planning node, can't find a path");
        return;
      }

      _Path   = _pathPlaner.getPath();
      _Radius = _pathPlaner.getRadius();

      /*cout<<"check path : \n"<<_Path<<endl;
      cout<<"check radius:\n"<<_Radius<<endl;*/

      visPath(_Path, _Radius);
      trajGenerationTest(_Path, _Radius);      
}

pair< vector<Vector3d>, vector<double> > PathDiscretize(MatrixXd path, double h)
{
      cout<<"check path: \n"<<path<<endl;

      vector<Vector3d> dis_path;
      vector<double>   dis_time;

      double t = 0.0;
      for(int i = 0; i < path.rows() - 1; i ++)
      {
          Vector3d point0 = path.row(i);
          Vector3d pointf = path.row(i + 1);
          Vector3d vec    = pointf - point0; 

          double dis = (pointf - point0).norm();
          Vector3d vec_dir = vec / dis; 

          int step = (int) (dis / _l);

          for(int m = 0; m < step; m ++)
          {
              Vector3d point  = point0 + vec_dir * m * _l;
              dis_path.push_back(point);
              dis_time.push_back(t);
              t += h;
          }

          dis_path.push_back(pointf);
          dis_time.push_back(t);
      }

      return make_pair(dis_path, dis_time);
}

double vecSum(VectorXd vec)
{   
    double sum = 0.0;
    for(int i = 0; i < vec.size(); i++)
      sum += vec(i);

    return sum;
}

double _scale;
void trajGenerationTest(MatrixXd path, VectorXd radius)
{           
      MatrixXd pos(2,3);
      MatrixXd vel(2,3);
      MatrixXd acc(2,3);

      pos <<  startPt.x,  startPt.y,  startPt.z,
              endPt.x,    endPt.y,    endPt.z;
      
      vel <<  startVel.x, startVel.y, startVel.z,
      //vel <<  -1.0,     -1.0,       0.5,
              0.0,       0.0,       0.0;

      acc <<  startAcc.x, startAcc.y, startAcc.z,
              0.0,        0.0,        0.0;
      
      MatrixXd path_   = path;
      VectorXd radius_ = radius;

      for(int i = 0; i < radius_.size(); i++)
          radius_(i) -= 0.1;

      _Time        = TimeAllocate( path_, radius_, vel.row(0) );
      _segment_num = radius_.size();

      double h;
      double max_v, max_a;
      
      max_a = _MAX_Acc;
      max_v = sqrt(_l * max_a);
      h = sqrt(4 * _l / max_a);

      pair< vector<Vector3d>, vector<double> > path_time = PathDiscretize(path, h);
      vector<Vector3d> path_dis = path_time.first;
      vector<double>   time_dis = path_time.second;

      visDiscretePath(path_dis);
      //ROS_WARN("[Benchmark] discrete path size is %d, discrete time index size is %d", (int)path_dis.size(), (int)time_dis.size());
/*      ROS_WARN("[Benchmark] check discrete path");
      for(auto ptr:path_dis)
        cout<<ptr<<endl;*/

/*      ROS_WARN("[Benchmark] check discrete time index");
      for(auto ptr:time_dis)
        cout<<ptr<<endl;*/

      MatrixXd acc_discrete;
      acc_discrete = _trajectoryGeneratorDiscrete.AccDiscreteGeneration(path_dis, vel, acc, max_v, max_a, h, _l); 
      ROS_WARN("[Benchmark] total time allocated in the discrete acc trajectory planning is %f", time_dis.back());
      ROS_WARN("[Benchmark] total time allocated in the Bezier and monomial trajectory planning is %f", vecSum(_Time));
      visDisTraj(acc_discrete, h);

      //cout<<"check time allocation:\n"<<_Time<<endl;
      
      _poly_orderList.clear();
      for(int i = 0; i < _segment_num; i ++ )
          _poly_orderList.push_back( 10 );

      _PolyCoeff_b_socp = _trajectoryGeneratorSocp.BezierPloyCoeffGenerationSOCP(  
                          path_, radius_, _Time, _poly_orderList, _FMList, pos, vel, acc, _MAX_Vel, _MAX_Acc, _minimize_order );

      if(_PolyCoeff_b_socp.rows() == 3 && _PolyCoeff_b_socp.cols() == 3)
          ROS_WARN("Bernstein solver failed");
      else
      {
          visBezierTrajectory(_PolyCoeff_b_socp, 0);
          visCtrlPoint(_PolyCoeff_b_socp, 0);
      }

      /*_PolyCoeff_b =  _trajectoryGeneratorBezier.BezierPloyCoeffGeneration(  
                      path_, radius_, _Time, _poly_orderList, _MQMList, pos, vel, acc, _MAX_Vel, _MAX_Acc, _minimize_order );
      visBezierTrajectory(_PolyCoeff_b, 1);
      visCtrlPoint(_PolyCoeff_b_socp, 1);*/

      ROS_WARN("finish Bezier trajectory solver");
      auto max_t = _Time.maxCoeff();
      _scale = (double)max_t;

      path_  /= _scale;
      radius_/= _scale;
      _Time  /= _scale;
      
      vector<MatrixXd> polyCoeff_pList;
      double weight_arc  = 0.0;
      double weight_jerk = 1.0;

      polyCoeff_pList =  _trajectoryGeneratorPoly.PloyCoeffGeneration(path_, radius_, _Time, vel, acc, _MAX_Vel, _MAX_Acc, 8, _scale, weight_arc, weight_jerk);

      if(polyCoeff_pList.back().rows() == 3 && polyCoeff_pList.back().cols() == 3)
          ROS_WARN("Monomial solver failed");
      else
      {
          visPolyTraj(polyCoeff_pList);
          visCheckPt(polyCoeff_pList);
      }

/*    ROS_WARN("[Benchmark] The global scale for monomial solver is %f", _scale);
      ROS_WARN("[Benchmark] Trajectory Cost in Berstein solver is %f", _trajectoryGeneratorSocp.getObjective());
      ROS_WARN("[Benchmark] Trajectory Cost in Monomial solver is %f", _trajectoryGeneratorPoly.getObjective().back());*/
}

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 & pointcloud_map )
{
      pcl::PointCloud<pcl::PointXYZ> CloudIn;      
      pcl::fromROSMsg(pointcloud_map, CloudIn);
      _pathPlaner.setInput(CloudIn);
}

int main (int argc, char** argv) 
{        
      ros::init (argc, argv, "benchmark_node");
      ros::NodeHandle n( "~" );

      n.param("discreteTraj/l",   _l,  0.1);

      n.param("dynamic/vec",   _Vel,  2.0);
      n.param("dynamic/acc",   _Acc,  1.0);
      n.param("dynamic/max_vec",   _MAX_Vel,  3.0);
      n.param("dynamic/max_acc",   _MAX_Acc,  1.5);
      n.param("optimization/poly_order_min", _poly_order_min,  5);
      n.param("optimization/poly_order_max", _poly_order_max,  10);
      n.param("optimization/minimize_order", _minimize_order,  3);

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

      // subcribed msgs
      _dest_pts_sub = n.subscribe( "waypoints", 1, rcvWaypointsCallback );
      _map_sub     = n.subscribe( "PointCloud", 1, rcvPointCloudCallBack);

      // publish visualize related msgs
      _disPath_vis_pub      = n.advertise<visualization_msgs::Marker>("disPath", 1);                      

      _ctrlPt_vis_pub       = n.advertise<visualization_msgs::Marker>("ctrlPt", 1);                      
      _traj_vis_pub         = n.advertise<visualization_msgs::Marker>("trajectory_vis", 1);
      _traj_bezier_vis_pub  = n.advertise<visualization_msgs::Marker>("bezier_trajectory_vis", 1);
      _checkPt_vis_pub      = n.advertise<visualization_msgs::MarkerArray>("checkpoint_vis", 1);

      _poly_traj_vis_pub    = n.advertise<visualization_msgs::Marker>("monomial_trajectory_vis", 1);
      _dis_traj_vis_pub     = n.advertise<visualization_msgs::Marker>("discrete_trajectory_vis", 1);

      _checkTraj_vis_pub    = n.advertise<visualization_msgs::Marker>("check_trajectory_vis", 1);
      _path_vis_pub         = n.advertise<visualization_msgs::MarkerArray>("path", 1);
      _traj_pub             = n.advertise<quadrotor_msgs::PolynomialTrajectory>("trajectory", 10);
     
      ros::Rate rate(100);
      bool status = ros::ok();
      while(status) 
      {
          ros::spinOnce();           
          status = ros::ok();
          rate.sleep();
      }

      return 0;
}

VectorXd TimeAllocate( MatrixXd Path, VectorXd Radius, Vector3d iniVel )
{
      int ball_num = Path.rows() - 2;

      //cout<<"ball num: "<<ball_num<<endl;

      MatrixXd check_pt( ball_num - 1 , 3 );
      MatrixXd Path_ball = Path.block( 1, 0, ball_num, 3 );

      //cout<<"Path_ball:\n"<<Path_ball<<endl;

      const static auto getdist = []( const Vector3d u, const Vector3d v ){
            const static auto sq = [] (double && val){ 
                  return val * val;
              };

            return sqrt( sq(v[0] - u[0]) + sq(v[1] - u[1]) + sq(v[2] - u[2])  );
      };

      Vector3d center, last_center;
      double radius, last_radius;

      int index = 0 ;
      last_center <<  Path_ball(index, 0), Path_ball(index, 1), Path_ball(index, 2); 
      last_radius  =  Radius[index];

      for(index = 1; index < ball_num; index ++ ){   
          center <<  Path_ball(index, 0), Path_ball(index, 1), Path_ball(index, 2); 
          radius  =  Radius[index];

          double dist = getdist(last_center, center);  
          
          Vector3d delta_Vec = center - last_center;
          Vector3d joint_pt;
          joint_pt = delta_Vec * ( dist + last_radius - radius) / ( 2.0 * dist ) + last_center; 

          check_pt.block( index - 1, 0, 1, 3 ) = joint_pt.transpose();

          last_center = center;
          last_radius = radius;
      }

      MatrixXd all_points( ball_num + 1 , 3 );
      all_points.row( 0 )                     = Path.row( 0 );
      all_points.block( 1, 0, ball_num-1, 3 ) = check_pt;
      all_points.row( ball_num )              = Path.row( ball_num + 1 );

      VectorXd time_allocate(all_points.rows() - 1);

      Vector3d initv = iniVel;
      for (int k = 0; k < all_points.rows() - 1; k++){
          // Position time
          double dtxyz;

          Vector3d p0   = all_points.row(k);           // The start point of this segment
          Vector3d p1   = all_points.row(k + 1);       // The end point of this segment
          Vector3d d    = p1 - p0;                     // The position difference
          Vector3d v0(0.0, 0.0, 0.0);                  // The init velocity
          
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

vector<double> getFullStateFromBezier(const MatrixXd & polyCoeff,  double t_now, int seg_now )
{
      vector<double > ret(9, 0);
      VectorXd ctrl_now = polyCoeff.row(seg_now);

      int order = _poly_orderList[seg_now];
      int ctrl_num1D = order + 1;

      for(int i = 0; i < 3; i++)
      {   
          for(int j = 0; j < ctrl_num1D; j++){
              ret[i] += _CList[order](j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (order - j) ); 
              
              if(j < ctrl_num1D - 1 )
                ret[i+3] += _CvList[order](j) * order 
                          * ( ctrl_now(i * ctrl_num1D + j + 1) - ctrl_now(i * ctrl_num1D + j))
                          * pow(t_now, j) * pow((1 - t_now), (order - j - 1) ); 
              
              if(j < ctrl_num1D - 2 )
                ret[i+6] += _CaList[order](j) * order * (order - 1) 
                          * ( ctrl_now(i * ctrl_num1D + j + 2) - 2 * ctrl_now(i * ctrl_num1D + j + 1) + ctrl_now(i * ctrl_num1D + j))
                          * pow(t_now, j) * pow((1 - t_now), (order - j - 2) );                         
          }

      }
      return ret;  
}

vector<double> getStateFromBezier(const MatrixXd & polyCoeff,  double t_now, int seg_now )
{
      vector<double > ret(3, 0);
      VectorXd ctrl_now = polyCoeff.row(seg_now);

      int order = _poly_orderList[seg_now];
      int ctrl_num1D = order + 1;
      for(int i = 0; i < 3; i++)
      {   
          for(int j = 0; j < ctrl_num1D; j++){
              ret[i] += _CList[order](j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (order - j) ); 
          }
      }

      //ROS_WARN("finish get state from Bezier curve ...");
      return ret;  
}

vector<double> getStateFromPolynomial( const MatrixXd & PolyCoeff, const VectorXd & Time, double t_now, int seg_now )
{
        vector<double > ret(3 * 3, 0);

        t_now = min(max(0.0, t_now), Time(seg_now) );
        int poly_num1D = PolyCoeff.cols() / 3;

        for ( int dim = 0; dim < 3; dim++ ){
            VectorXd coeff = (PolyCoeff.row(seg_now)).segment( dim * poly_num1D, poly_num1D );

            MatrixXd t = MatrixXd::Zero( 3, poly_num1D );
            
            for(int j = 0; j < poly_num1D; j ++){
                t(0 ,j) = pow(t_now, j);
                t(1 ,j) = j * pow(t_now, j-1);
                t(2 ,j) = j * (j - 1) * pow(t_now, j-2);
            }

            for (int i = 0; i < 3; i++)
                ret[dim + i * 3] = coeff.dot(t.row(i));
        }

      return ret;
}

vector<double> getFullStateFromPolynomial( const MatrixXd & PolyCoeff, const VectorXd & Time, double t_now, int seg_now )
{
        vector<double > ret(3 * 3, 0);

        t_now = min(max(0.0, t_now), Time(seg_now) );
        int poly_num1D = PolyCoeff.cols() / 3;

        for ( int dim = 0; dim < 3; dim++ ){
            VectorXd coeff = (PolyCoeff.row(seg_now)).segment( dim * poly_num1D, poly_num1D );

            MatrixXd t = MatrixXd::Zero( 3, poly_num1D );
            
            for(int j = 0; j < poly_num1D; j ++){
                t(0 ,j) = pow(t_now, j);
                t(1 ,j) = j * pow(t_now, j-1);
                t(2 ,j) = j * (j - 1) * pow(t_now, j-2);
            }

            for (int i = 0; i < 3; i++)
                ret[dim + i * 3] = coeff.dot(t.row(i));
        }

      return ret;
}

void visPolyTraj(vector<MatrixXd> polyCoeffList)
{        
    for(int k = 0; k < (int)polyCoeffList.size(); k++)
    {   
        visualization_msgs::Marker _traj_vis;

        _traj_vis.header.stamp       = ros::Time::now();
        _traj_vis.header.frame_id    = "map";

        string ns = "benchmark/trajectory" + to_string(k);
        _traj_vis.ns = ns;
        _traj_vis.id = k;
        _traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
        _traj_vis.action = visualization_msgs::Marker::ADD;
        _traj_vis.scale.x = _vis_traj_width;
        _traj_vis.scale.y = _vis_traj_width;
        _traj_vis.scale.z = _vis_traj_width;
        _traj_vis.pose.orientation.x = 0.0;
        _traj_vis.pose.orientation.y = 0.0;
        _traj_vis.pose.orientation.z = 0.0;
        _traj_vis.pose.orientation.w = 1.0;

        _traj_vis.color.r = 1.0 - (double)k / (int)polyCoeffList.size();
        _traj_vis.color.g = 0.0;
        _traj_vis.color.b = 0.0;
        _traj_vis.color.a = 1.0;

        double traj_len = 0.0;
        int count = 0;
        Vector3d cur, pre;
        cur.setZero();
        pre.setZero();

        _traj_vis.points.clear();
        vector<double> state;
        geometry_msgs::Point pt;
        Vector3d max_v(0.0, 0.0, 0.0);
        Vector3d max_a(0.0, 0.0, 0.0);

        for(int i = 0; i < _Time.size(); i++ )
        {   
            for (double t = 0.0; t < _Time(i); t += 0.02/_scale, count += 1)
            {
                state = getFullStateFromPolynomial(polyCoeffList[k], _Time, t, i);
                cur(0) = pt.x = _scale * state[0];
                cur(1) = pt.y = _scale * state[1];
                cur(2) = pt.z = _scale * state[2];
                _traj_vis.points.push_back(pt);

                if(t > 0.0)
                  for(int p = 0; p < 3; p++)
                  {   
                      double v_p = fabs(state[3 + p]);
                      if( v_p > max_v(p))
                        max_v(p) = v_p;

                      double a_p = fabs(state[6 + p]) / _scale;
                      if( a_p > max_a(p))
                        max_a(p) =  a_p;
                  }
                
                if (count) traj_len += (pre - cur).norm();
                pre = cur;
            }
        }

        ROS_WARN(" max velocity of monomial curve in x, y, z axis are %f, %f, %f ",     max_v(0), max_v(1), max_v(2));
        ROS_WARN(" max acceleration of monomial curve in x, y, z axis are %f, %f, %f ", max_a(0), max_a(1), max_a(2));
        ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);
        _poly_traj_vis_pub.publish(_traj_vis);

    }
}

void visDisTraj(MatrixXd disAcc, double h)
{        
    visualization_msgs::Marker _traj_vis;

    _traj_vis.header.stamp       = ros::Time::now();
    _traj_vis.header.frame_id    = "map";

    _traj_vis.ns = "benchmark/discrete_trajectory";
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
    _traj_vis.color.a = 1.0;

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();

    _traj_vis.points.clear();
    vector<double> state;
    geometry_msgs::Point pt;

    Vector3d pos, vel;
    pos << startPt.x,  startPt.y,  startPt.z;
    vel << startVel.x, startVel.y, startVel.z;
    
    Vector3d max_v(0.0, 0.0, 0.0);
    Vector3d max_a(0.0, 0.0, 0.0);
    
    for(int i = 0; i < disAcc.cols(); i++)
    {   
        Vector3d acc_t = disAcc.col(i);
        Vector3d pos_t, vel_t;

        for(double t = 0.0; t < h; t += 0.01, count ++)
        {
            pos_t = pos + vel * t + acc_t * (t*t/2.0);
            vel_t = vel + acc_t * t;

            cur(0) = pt.x = pos_t(0);
            cur(1) = pt.y = pos_t(1);
            cur(2) = pt.z = pos_t(2);
            _traj_vis.points.push_back(pt);

            for(int p = 0; p < 3; p++)
            {   
                double v_p = fabs(vel_t(p));
                if( v_p > max_v(p))
                  max_v(p) = v_p;

                double a_p = fabs(acc_t(p));
                if( a_p > max_a(p))
                  max_a(p) =  a_p;
            }
            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
        pos = pos_t;
        vel = vel_t;
    }

    ROS_WARN(" max velocity of discrete curve in x, y, z axis are %f, %f, %f ",     max_v(0), max_v(1), max_v(2));
    ROS_WARN(" max acceleration of discrete curve in x, y, z axis are %f, %f, %f ", max_a(0), max_a(1), max_a(2));
    ROS_INFO(" The length of the discrete trajectory; %.3lfm.", traj_len);
    _dis_traj_vis_pub.publish(_traj_vis);
}

void visBezierTrajectory(MatrixXd polyCoeff, int rgb)
{
    visualization_msgs::Marker _traj_vis;

    _traj_vis.header.stamp       = ros::Time::now();
    _traj_vis.header.frame_id    = "map";

    _traj_vis.ns = "benchmark/bezier_trajectory" + to_string(rgb);
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

    if(rgb == 0)
    {
      _traj_vis.color.r = 0.0;
      _traj_vis.color.g = 1.0;
      _traj_vis.color.b = 0.0;
      _traj_vis.color.a = 0.3;
    }
    else
    {
      _traj_vis.color.r = 0.0;
      _traj_vis.color.g = 0.0;
      _traj_vis.color.b = 1.0;
      _traj_vis.color.a = 0.3; 
    }

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    _traj_vis.points.clear();

    vector<double> state;
    geometry_msgs::Point pt;
    Vector3d max_v(0.0, 0.0, 0.0);
    Vector3d max_a(0.0, 0.0, 0.0);

    for(int i = 0; i < _segment_num; i++ ){
        for (double t = 0.0; t < 1.0; t += 0.01 / _Time(i), count += 1){
            state = getFullStateFromBezier( polyCoeff, t, i );
            cur(0) = pt.x = _Time(i) * state[0];
            cur(1) = pt.y = _Time(i) * state[1];
            cur(2) = pt.z = _Time(i) * state[2];
            _traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;

            for(int p = 0; p < 3; p++)
            {   
                double v_ = fabs(state[3 + p]);
                double a_ = fabs(state[6 + p]) / _Time(i);
                
                if( v_ > max_v(p) )
                    max_v(p) = v_;

                if( a_ > max_a(p) )
                    max_a(p) = a_;
            }
        }
    }

    ROS_WARN(" max velocity of the Bezier curve in x, y, z axis are %f, %f, %f ",     max_v(0), max_v(1), max_v(2));
    ROS_WARN(" max acceleration of the Bezier curve in x, y, z axis are %f, %f, %f ", max_a(0), max_a(1), max_a(2));
    ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);

    _traj_bezier_vis_pub.publish(_traj_vis);
}

void visPath(MatrixXd path, VectorXd radius)
{     
    MatrixXd path_new  = path.block(1, 0, path.rows() - 2, path.cols() );

    for (auto & mk: path_vis.markers) 
      mk.action = visualization_msgs::Marker::DELETE;

    _path_vis_pub.publish(path_vis);
    path_vis.markers.clear();

    visualization_msgs::Marker mk;
    mk.header.frame_id = "map";
    mk.header.stamp = ros::Time::now();
    mk.ns = "benchmark/path";
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

void visDiscretePath(vector<Vector3d> path_dis)
{
    visualization_msgs::Marker disPath_vis;

    disPath_vis.header.frame_id = "map";
    disPath_vis.header.stamp = ros::Time::now();
    disPath_vis.ns = "benchmark/disPath";
    disPath_vis.type = visualization_msgs::Marker::CUBE_LIST;
    disPath_vis.action = visualization_msgs::Marker::ADD;
    disPath_vis.pose.orientation.x = 0.0;
    disPath_vis.pose.orientation.y = 0.0;
    disPath_vis.pose.orientation.z = 0.0;
    disPath_vis.pose.orientation.w = 1.0;
    disPath_vis.id = 0;
    disPath_vis.scale.x = 2.0 * _l;
    disPath_vis.scale.y = 2.0 * _l;
    disPath_vis.scale.z = 2.0 * _l;

    disPath_vis.color.r = 0.0;
    disPath_vis.color.g = 1.0;
    disPath_vis.color.b = 0.0;
    disPath_vis.color.a = 0.3;
    
    disPath_vis.points.clear();
    geometry_msgs::Point pt;
    for(int i = 0; i < (int)path_dis.size(); i++)
    {   
        pt.x = path_dis[i](0);
        pt.y = path_dis[i](1);
        pt.z = path_dis[i](2);
        disPath_vis.points.push_back(pt);
    }

    _disPath_vis_pub.publish(disPath_vis);
}

void visCtrlPoint(MatrixXd polyCoeff, int rgb)
{
    visualization_msgs::Marker ctrlPt_vis;

    ctrlPt_vis.header.frame_id = "map";
    ctrlPt_vis.header.stamp = ros::Time::now();
    ctrlPt_vis.ns = "benchmark/ctrlPt" + to_string(rgb);
    ctrlPt_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    ctrlPt_vis.action = visualization_msgs::Marker::ADD;
    ctrlPt_vis.pose.orientation.x = 0.0;
    ctrlPt_vis.pose.orientation.y = 0.0;
    ctrlPt_vis.pose.orientation.z = 0.0;
    ctrlPt_vis.pose.orientation.w = 1.0;
    ctrlPt_vis.id = 0;
    ctrlPt_vis.scale.x = 2.0 * _vis_traj_width;
    ctrlPt_vis.scale.y = 2.0 * _vis_traj_width;
    ctrlPt_vis.scale.z = 2.0 * _vis_traj_width;

    if(rgb == 0)
    {
      ctrlPt_vis.color.r = 0.0;
      ctrlPt_vis.color.g = 1.0;
      ctrlPt_vis.color.b = 0.0;
      ctrlPt_vis.color.a = 0.3;
    }
    else
    {
      ctrlPt_vis.color.r = 0.0;
      ctrlPt_vis.color.g = 0.0;
      ctrlPt_vis.color.b = 1.0;
      ctrlPt_vis.color.a = 0.3; 
    }

/*      mk.color.a = 1.0;
    mk.color.r = 1.0;
    mk.color.g = 0.0;
    mk.color.b = 0.0;
*/    
    ctrlPt_vis.points.clear();
    geometry_msgs::Point pt;
    for(int i = 0; i < _segment_num; i++)
    {   
        int order = _poly_orderList[i];
        int ctrl_num = order + 1;
        
        for(int j = 0; j < ctrl_num; j++)
        {
            pt.x = _Time(i) * polyCoeff(i, j);
            pt.y = _Time(i) * polyCoeff(i, ctrl_num + j);
            pt.z = _Time(i) * polyCoeff(i, 2 * ctrl_num + j);
            ctrlPt_vis.points.push_back(pt);
        }
    }

    _ctrlPt_vis_pub.publish(ctrlPt_vis);
    //ROS_WARN("vis control points OK");
}

visualization_msgs::MarkerArray checkPt_vis;
void visCheckPt(vector<MatrixXd> polyCoeffList)
{     
    MatrixXd poly = polyCoeffList.back();

    for (auto & mk: checkPt_vis.markers) 
        mk.action = visualization_msgs::Marker::DELETE;
    
    _checkPt_vis_pub.publish(checkPt_vis);
    checkPt_vis.markers.clear();

    visualization_msgs::Marker mk;
    mk.header.frame_id = "map";
    mk.header.stamp = ros::Time::now();
    mk.ns = "benchamrk_test/checkPt";
    mk.type = visualization_msgs::Marker::SPHERE;
    mk.action = visualization_msgs::Marker::ADD;
    mk.pose.orientation.x = 0.0;
    mk.pose.orientation.y = 0.0;
    mk.pose.orientation.z = 0.0;
    mk.pose.orientation.w = 1.0;
    mk.color.a = 1.0;
    mk.color.r = 0.0;
    mk.color.g = 0.0;
    mk.color.b = 0.0;

    vector<double> state;
    
    for(int i = 0; i < _Time.size(); i++){
        mk.id = i;
        
        state = getStateFromPolynomial(poly, _Time, _Time(i), i);
        mk.pose.position.x = _scale * state[0];
        mk.pose.position.y = _scale * state[1];
        mk.pose.position.z = _scale * state[2];
        mk.scale.x = 0.25;
        mk.scale.y = 0.25;
        mk.scale.z = 0.25;
        
        checkPt_vis.markers.push_back(mk);
    }

    _checkPt_vis_pub.publish(checkPt_vis);
    ROS_WARN("vis check_pt OK");
}