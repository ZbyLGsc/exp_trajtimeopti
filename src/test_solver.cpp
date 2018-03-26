#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/Imu.h>
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
#include "pcd_trajectory/mosek.h"
//#include "pcd_trajectory/trajectory_generator.h"
#include "pcd_trajectory/bezier_base.h"
#include "pcd_trajectory/trajectory_generator_socp.h"

using namespace std;

void visBezierTrajectory(const ros::TimerEvent& evt);
void visualize_path(const ros::TimerEvent& evt);

vector<double> getDesiredState( const Eigen::MatrixXd & PolyCoeff, const Eigen::VectorXd & Time, double t_now, int seg_now );

vector<double> getStateFromBezier(const Eigen::MatrixXd & polyCoeff,  double t_now, int seg_now );

visualization_msgs::Marker _traj_vis;
visualization_msgs::MarkerArray path_vis;
visualization_msgs::MarkerArray checkPt_vis;

ros::Publisher _traj_vis_pub;
ros::Publisher path_vis_pub;
ros::Publisher checkPt_vis_pub;

ros::Timer vis_timer;
ros::Timer vis_timer2;

double _vis_traj_width = 0.25; 
double _MAX_Acc = 1.0;
double _MAX_Vel = 2.0;

Eigen::MatrixXd PolyCoeff;
Eigen::VectorXd _Time;
Eigen::VectorXd Radius;
Eigen::MatrixXd Path;
int _segment_num = 14; 

int _traj_order = 7;
int min_order   = 3;
VectorXd _C, _Cv, _Ca, _Cj;

int main (int argc, char** argv) {
        
      ros::init (argc, argv, "SOCP_SOLVER_TEST");
      ros::NodeHandle n( "~" );
      
      _traj_vis_pub =
            n.advertise<visualization_msgs::Marker>("trajectory_vis", 2);

      path_vis_pub =
            n.advertise<visualization_msgs::MarkerArray>("path", 2);

      checkPt_vis_pub =
            n.advertise<visualization_msgs::MarkerArray>("checkPt", 2);

      vis_timer = 
            n.createTimer(ros::Duration(1.0), &visBezierTrajectory);
      
      vis_timer2 = 
            n.createTimer(ros::Duration(1.0), &visualize_path);

      TrajectoryGenerator _trajGen;

      Eigen::MatrixXd Pos(2,3);
      Eigen::MatrixXd Vel(2,3);
      Eigen::MatrixXd Acc(2,3);
      
      double maxVel = _MAX_Vel;
      double maxAcc = _MAX_Acc;

      double coeff = 1.0;

      Radius.resize(_segment_num); 
      //Radius << 2.5, 2.8, 3.0;
      Radius <<1.14099, 0.719942, 0.964115, 2.10253, 1.88788, 1.78146, 1.15229, 0.881995, 0.948956, 1.43287, 1.80342, 6.9662, 9.65, 9.65;

      Radius /= coeff;
      
      _Time.resize(_segment_num);   
      //_Time << 1.0, 1.0, 1.0;
      _Time << 1.01999, 1.13218, 1.35922, 2.05877, 1.96878, 2.00296, 1.53475, 1.2645, 1.34313, 1.77562, 1.9258, 3.39531, 3.70139, 4.17193;

      _Time = _Time / coeff;

      cout<<_Time<<endl;
      Path.resize( _segment_num + 2, 3 );
      
/*      Path << 1.0,  0.0,  1.2,
              1.0,  1.2,  1.8,
              2.5,  3.0,  4.5,
              4.0,  3.7,  5.2,
              5.0,  4.0,  5.0;*/

      Path << 0,       0, 1.98827,
              0,       0, 1.98827,
              0.437104, 1.05237, 2.04592,
              1.44795, 2.12849, 1.9552,
              3.42961,  3.0739,  2.5091,
              5.62685, 4.83824, 2.59805,
              7.88225, 6.97777, 2.43901,
              9.29645, 8.92037, 2.29213,
              10.4435, 9.85909,  1.9188,
              11.7939, 9.89528,   1.851,
              13.5321, 10.3058, 1.88601,
              16.3017,  10.843, 2.15074,
              21.2849, 16.9461, 2.06462,
              26.2746, 21.7823, 2.55642,
              33.1271, 30.3311, 2.46626,
              38.4077, 35.7515,    2.36;
      Path /= coeff;
      
      /*Pos << 1.0, 0.0, 1.2,
             5.0, 4.0, 5.0;*/

      Pos << 0,        0,        1.98827,
             38.4077,  35.7515,  2.36;

      Vel << 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0;
      
      Acc << 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0;
      
      Bernstein _bernstein;
      if(_bernstein.setParam(3, 12, min_order) == -1)
      {
          ROS_ERROR(" The trajectory order is set beyond the library's scope, please re-set ");
      }

      Eigen::MatrixXd MQM = _bernstein.getMQM()[_traj_order];
      Eigen::MatrixXd FM  = _bernstein.getFM()[_traj_order];

      _C   = _bernstein.getC()[_traj_order];
      _Cv  = _bernstein.getC_v()[_traj_order];
      _Ca  = _bernstein.getC_a()[_traj_order];
      
      
      ROS_WARN("traj qcqp solver");
      PolyCoeff = _trajGen.BezierPloyCoeffGenerationQCQP( Path, Radius, _Time, MQM, Pos, Vel, Acc, maxVel, maxAcc, _traj_order, min_order );
      
      /*ROS_WARN("traj socp solver");
      PolyCoeff = _trajGen.BezierPloyCoeffGenerationSOCP( Path, Radius, _Time, FM, Pos, Vel, Acc, maxVel, maxAcc, _traj_order, 3 );*/

      ROS_WARN("traj solver finished");
      ros::spin ();
}

vector<double> getStateFromBezier(const Eigen::MatrixXd & polyCoeff,  double t_now, int seg_now )
{
    vector<double > ret(9, 0);
    Eigen::VectorXd ctrl_now = polyCoeff.row(seg_now);
    int ctrl_num1D = polyCoeff.cols() / 3;

    for(int i = 0; i < 3; i++)
    {   
        for(int j = 0; j < ctrl_num1D; j++){
            ret[i] += _C(j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (_traj_order - j) ); 
          
            if(j < ctrl_num1D - 1 )
                ret[i+3] += _Cv(j) * _traj_order 
                      * ( ctrl_now(i * ctrl_num1D + j + 1) - ctrl_now(i * ctrl_num1D + j))
                      * pow(t_now, j) * pow((1 - t_now), (_traj_order - j - 1) ); 
          
            if(j < ctrl_num1D - 2 )
                ret[i+6] += _Ca(j) * _traj_order * (_traj_order - 1) 
                      * ( ctrl_now(i * ctrl_num1D + j + 2) - 2 * ctrl_now(i * ctrl_num1D + j + 1) + ctrl_now(i * ctrl_num1D + j))
                      * pow(t_now, j) * pow((1 - t_now), (_traj_order - j - 2) );                         
        }
    }

      //ROS_WARN("finish get state from Bezier curve ...");
      return ret;  
}


void visBezierTrajectory(const ros::TimerEvent& evt)
{   
    visualization_msgs::Marker _traj_vis;

    _traj_vis.header.stamp       = ros::Time::now();
    _traj_vis.header.frame_id    = "map";

    _traj_vis.ns = "trajectory/trajectory";
    _traj_vis.id = 0;
    _traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    
    _traj_vis.action = visualization_msgs::Marker::DELETE;

    _traj_vis.action = visualization_msgs::Marker::ADD;
    _traj_vis.scale.x = _vis_traj_width;
    _traj_vis.scale.y = _vis_traj_width;
    _traj_vis.scale.z = _vis_traj_width;
    _traj_vis.pose.orientation.x = 0.0;
    _traj_vis.pose.orientation.y = 0.0;
    _traj_vis.pose.orientation.z = 0.0;
    _traj_vis.pose.orientation.w = 1.0;
    _traj_vis.color.r = 1.0;
    _traj_vis.color.g = 0.0;
    _traj_vis.color.b = 0.0;
    _traj_vis.color.a = 0.6;

    double traj_len = 0.0;
    int count = 0;
    Eigen::Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    _traj_vis.points.clear();

    vector<double> state;
    geometry_msgs::Point pt;

    for(int i = 0; i < _segment_num; i++ ){
        for (double t = 0.0; t < 1.0; t += 0.005 / _Time[i], count += 1){
            state = getStateFromBezier( PolyCoeff, t, i );
            cur(0) = pt.x = _Time[i] * state[0];
            cur(1) = pt.y = _Time[i] * state[1];
            cur(2) = pt.z = _Time[i] * state[2];
            _traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);

    _traj_vis_pub.publish(_traj_vis);
}

void visualize_path(const ros::TimerEvent& evt)
{
      path_vis.markers.clear();
      
      visualization_msgs::Marker mk;
      mk.header.frame_id = "/base_link";
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

      int seg_num = Radius.size();

      for(int i = 0; i < seg_num; i++){
          mk.id = i;
          mk.pose.position.x = Path(i, 0); 
          mk.pose.position.y = Path(i, 1); 
          mk.pose.position.z = Path(i, 2); 
          mk.scale.x = 2 * Radius(i);
          mk.scale.y = 2 * Radius(i);
          mk.scale.z = 2 * Radius(i);
          
          path_vis.markers.push_back(mk);
    }

    path_vis_pub.publish(path_vis);
}
