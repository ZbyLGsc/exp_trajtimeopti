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
#include "pcd_trajectory/trajectory_generator.h"

using namespace std;

void visTrajectory( const ros::TimerEvent& evt );
void visualize_path(const ros::TimerEvent& evt);
void visualize_checkPt(const ros::TimerEvent& evt);

vector<double> getDesiredState( const Eigen::MatrixXd & PolyCoeff, const Eigen::VectorXd & Time, double t_now, int seg_now );

visualization_msgs::Marker _traj_vis;
visualization_msgs::MarkerArray path_vis;
visualization_msgs::MarkerArray checkPt_vis;

ros::Publisher _traj_vis_pub;
ros::Publisher path_vis_pub;
ros::Publisher checkPt_vis_pub;

ros::Timer vis_timer;
ros::Timer vis_timer2;
ros::Timer vis_timer3;

double _vis_traj_width = 0.25; 
double _MAX_Acc = 1.0;
double _MAX_Vel = 2.0;

Eigen::MatrixXd PolyCoeff;
Eigen::VectorXd Time;
Eigen::VectorXd Radius;
Eigen::MatrixXd Path;

int main (int argc, char** argv) {
        
      ros::init (argc, argv, "pcd_RRT");
      ros::NodeHandle n( "~" );
      
      _traj_vis_pub =
            n.advertise<visualization_msgs::Marker>("trajectory_vis", 2);

      path_vis_pub =
            n.advertise<visualization_msgs::MarkerArray>("path", 2);

      checkPt_vis_pub =
            n.advertise<visualization_msgs::MarkerArray>("checkPt", 2);

      vis_timer = 
            n.createTimer(ros::Duration(1.0), &visTrajectory);
      
      vis_timer2 = 
            n.createTimer(ros::Duration(1.0), &visualize_path);

      vis_timer3 = 
            n.createTimer(ros::Duration(1.0), &visualize_checkPt);

      TrajectoryGenerator _trajGen;

      Eigen::MatrixXd Vel(2,3);
      Eigen::MatrixXd Acc(2,3);
      
      double maxVel = _MAX_Vel;
      double maxAcc = _MAX_Acc;

      int seg_num = 18;
      double coeff = 1.0;

      Radius.resize(seg_num); 
      Radius << 14.1476, 11.5186, 3.229, 3.78748, 4.224, 3.54606, 4.61302, 4.91227, 4.38179, 4.46045, 8.20993, 8.38996, 4.18431, 5.53876, 3.54133, 2.55436, 2.45662, 3.04125 ;
      //Radius << 14.1476, 11.5186, 13.229, 13.78748, 14.224, 13.54606, 14.61302, 14.91227, 14.38179, 14.46045, 18.20993, 18.38996, 14.18431, 15.53876, 13.54133, 12.55436, 12.45662, 3.04125 ;
      Radius /= coeff;
      
      Time.resize(seg_num);   
      Time << 8.3883, 15.6335, 3.76884, 4.41039, 4.48409, 4.04422, 6.46784, 6.6743, 5.05134, 5.01699, 12.8626, 13.2236, 4.78522, 9.01453, 4.92612, 2.60324, 2.1643, 3.22738;
      //Time << 8.3883, 15.6335, 3.76884, 4.41039, 4.48409, 4.04422, 6.46784, 6.6743, 5.05134, 5.01699, 12.8626, 13.2236, 4.78522, 9.01453, 4.92612, 3.60324, 3.1643, 4.22738;
      //Time << 8.3883, 15.6335, 13.76884, 14.41039, 14.48409, 14.04422, 16.46784, 16.6743, 15.05134, 15.01699, 12.8626, 13.2236, 14.78522, 19.01453, 14.92612, 12.60324, 12.1643, 13.22738;
      Time = Time / coeff;

      Path.resize( seg_num, 3 );
      
      Path << -30.0000,  40.0000 , 15.0000,
              -18.7172,  31.6776 , 13.1051,
              -9.87129,  24.3038 , 12.8694,
              -5.39006,  22.3146 , 12.1797,
              -1.80731,  21.4359 , 11.3214,
               1.83737,  19.3013 , 11.3657,
               7.88089,  21.2065 , 11.3759,
               13.9057,  20.8183 , 11.9702,
               20.2955,  19.6098 , 11.7409,
               24.5734,  20.3927 , 11.205,
               32.2655,  14.253  , 15.0732,
               32.8509, -0.424517, 14.3333,
               28.7394, -7.71215 , 13.719,
               29.7737, -14.1911 , 10.1755,
               30.9423, -21.545  , 9.15894,
               32.3429, -24.4636 , 7.72331,
               33.3533, -26.5688 , 6.6878,
               34.3249, -28.5934 , 5.69191;
      
      Path /= coeff;
      
      Vel << 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0;
      
      Acc << 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0;
      
      
      PolyCoeff = _trajGen.PloyCoeffGeneration( Path, Path, Radius, Radius, Time, Vel, Acc, maxVel, maxAcc );
      
      //visTrajectory(PolyCoeff, Time);

      ros::spin ();
}

bool pub_once = true;
void visTrajectory( const ros::TimerEvent& evt )
{        
        _traj_vis.header.stamp       = ros::Time::now();
        _traj_vis.header.frame_id    = "/base_link";

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
        _traj_vis.color.g = 0.8;
        _traj_vis.color.b = 0.0;
        _traj_vis.color.a = 1.0;

        const static auto sumT = [](Eigen::VectorXd T){
            double T_sum = 0.0;
            for(int i = 0; i < T.size(); i++)
              T_sum += T(i);

            return T_sum; 
        };

        double traj_len = 0.0;
        int count = 0;
        Eigen::Vector3d cur, pre;
        cur.setZero();
        pre.setZero();
        
        double t_begin = 0.0, t_final = sumT(Time);

        _traj_vis.points.clear();
        _traj_vis.points.reserve(static_cast<int>((t_final - t_begin) * 100 + 0.5));
        vector<double> state;
        geometry_msgs::Point pt;

        //ROS_INFO("[GENERATOR] Trajectory visualization prepared.");

        int seg_num = Time.size();
        for(int i = 0; i < seg_num; i++ ){
            for (double t = t_begin; t < Time(i); t += 0.02, count += 1){
                state = getDesiredState(PolyCoeff, Time, t, i);
/*                for(auto ptr:state)
                  cout<<ptr<<endl;*/
                cur(0) = pt.x = state[0];
                cur(1) = pt.y = state[1];
                cur(2) = pt.z = state[2];
                _traj_vis.points.push_back(pt);

                if (count) traj_len += (pre - cur).norm();
                pre = cur;
            }
        }

        if(pub_once){
          ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);
          pub_once = false;
        }

        _traj_vis_pub.publish(_traj_vis);
}


vector<double> getDesiredState( const Eigen::MatrixXd & PolyCoeff, const Eigen::VectorXd & Time, double t_now, int seg_now )
{

        vector<double > ret(3 * 3, 0);

        t_now = min(max(0.0, t_now), Time(seg_now) );
        int poly_num1D = PolyCoeff.cols() / 3;

        for ( int dim = 0; dim < 3; dim++ ){
            Eigen::VectorXd coeff = (PolyCoeff.row(seg_now)).segment( dim * poly_num1D, poly_num1D );

            Eigen::MatrixXd t = Eigen::MatrixXd::Zero( 3, poly_num1D );
            
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

void visualize_checkPt(const ros::TimerEvent& evt)
{
      checkPt_vis.markers.clear();
      
      visualization_msgs::Marker mk;
      mk.header.frame_id = "/base_link";
      mk.header.stamp = ros::Time::now();
      mk.ns = "pcd_RRT/checkPt";
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

      int seg_num = Radius.size();
      vector<double> state;
      
      for(int i = 0; i < seg_num; i++){
          mk.id = i;
          
          state = getDesiredState(PolyCoeff, Time, Time(i), i);
          mk.pose.position.x = state[0];
          mk.pose.position.y = state[1];
          mk.pose.position.z = state[2];
          mk.scale.x = 0.25;
          mk.scale.y = 0.25;
          mk.scale.z = 0.25;
          
          checkPt_vis.markers.push_back(mk);
    }

    checkPt_vis_pub.publish(checkPt_vis);
}
