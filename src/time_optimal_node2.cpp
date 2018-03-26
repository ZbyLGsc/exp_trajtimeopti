//#define PCL_NO_PRECOMPILE
#include <pcl/filters/voxel_grid.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <iostream>

#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/kdtree.h>

#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <ros/console.h>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/PointCloud2.h>
#include <fstream>

#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <pcl/search/impl/kdtree.hpp>
#include <random>

#include <tf/transform_broadcaster.h>
#include "pcd_trajectory/mosek.h"
#include "tf/tf.h"
#include "tf/transform_datatypes.h"
#include "visualization_msgs/Marker.h"
#include "visualization_msgs/MarkerArray.h"

#include "pcd_trajectory/backward.hpp"
#include "pcd_trajectory/bezier_base.h"
#include "pcd_trajectory/rrgPathFinder.h"
#include "pcd_trajectory/time_optimizer.h"
#include "pcd_trajectory/trajectory_generator.h"
#include "pcd_trajectory/trajectory_generator_lite.h"
#include "pcd_trajectory/trajectory_generator_socp_lite.h"
#include "quadrotor_msgs/PolynomialTrajectory.h"
#include "quadrotor_msgs/PositionCommand.h"

#include "grad_traj_optimizer.h"

namespace backward
{
backward::SignalHandling sh;
}

using namespace std;
using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;

visualization_msgs::MarkerArray path_vis;

ros::Publisher _traj_vis_pub;
ros::Publisher _traj_bezier_vis_pub;
ros::Publisher _checkPt_vis_pub;
ros::Publisher _poly_traj_vis_pub;
ros::Publisher _checkTraj_vis_pub;
ros::Publisher _path_vis_pub;
ros::Publisher _ctrlPt_vis_pub;
ros::Publisher _traj_pub;

ros::Publisher _vis_pos_pub;
ros::Publisher _vis_vel_pub;
ros::Publisher _vis_acc_pub;
ros::Publisher _cmd_pub;

ros::Subscriber _dest_pts_sub, _map_sub, _odom_sub;

nav_msgs::Path _waypoints;
nav_msgs::Odometry _odom;
bool _has_odom = false;
vector<double> _state;

float _vis_traj_width = 0.15;

double _MAX_Vel, _MAX_Acc, _Vel, _Acc, _eps;
int _minimize_order, _poly_order_min, _poly_order_max;
int _segment_num;
vector<int> _poly_orderList;
Eigen::MatrixXd _PolyCoeff_p;
int _poly_num1D;

Eigen::VectorXd _Time;
Eigen::VectorXd _Radius;
Eigen::MatrixXd _Path;

Point startPt = { 0.0, 0.0, 1.0 };  // For outdoor CYT model .
Point startVel = { 0.0, 0.0, 0.0 };
Point startAcc = { 0.0, 0.0, 0.0 };
Point endPt;

vector<MatrixXd> _MQMList, _FMList;
vector<VectorXd> _CList, _CvList, _CaList;

vector<double> getStateFromBezier(const Eigen::MatrixXd &polyCoeff, double t_now, int seg_now);
vector<double> getFullStateFromBezier(const Eigen::MatrixXd &polyCoeff, double t_now, int seg_now);
vector<double> getStateFromPolynomial(const Eigen::MatrixXd &PolyCoeff, const Eigen::VectorXd &Time,
                                      double t_now, int seg_now);
vector<double> getFullStateFromPolynomial(const Eigen::MatrixXd &PolyCoeff,
                                          const Eigen::VectorXd &Time, double t_now, int seg_now);
Vector3d getPosPoly(int k, double s);
Vector3d getVelPoly(int k, double s);
Vector3d getAccPoly(int k, double s);

void visBezierTrajectory(Eigen::MatrixXd polyCoeff, int rgb);
void visCtrlPoint(Eigen::MatrixXd polyCoeff, int rgb);
void visPath(Eigen::MatrixXd path, Eigen::VectorXd radius);
int trajGenerator(Eigen::MatrixXd path, Eigen::VectorXd radius, double time_pre);
void trajGenerationTest(Eigen::MatrixXd path, Eigen::VectorXd radius);
void visPolyTraj(vector<Eigen::MatrixXd> polyCoeffList);
void visCheckPt(vector<Eigen::MatrixXd> polyCoeffList);
void plotData(Allocator *time_allocator);

double x_l = -40.0, x_h = 40.0, y_l = -40.0, y_h = 40.0, z_l = 0.0, z_h = 5.0, bias_l = 0.0,
       bias_h = 1.0;
double mag_coeff;
rrgPathFinder _pathPlaner(x_l, x_h, y_l, y_h, z_l, z_h, bias_l, bias_h);
Allocator *time_allocator = new Allocator();

TrajectoryGeneratorPoly _trajectoryGeneratorPoly;
TrajectoryGeneratorBezier _trajectoryGeneratorBezier;
TrajectoryGeneratorSOCP _trajectoryGeneratorSocp;

bool traj_finish = false;
ros::Time traj_time_start;
ros::Time traj_time_final;
ros::Time odom_time;

visualization_msgs::Marker _vis_pos, _vis_vel, _vis_acc;
Vector3d position, velocity, acceleration;

quadrotor_msgs::PositionCommand _cmd;
double pos_gain[3] = { 5.7, 5.7, 6.2 };
double vel_gain[3] = { 3.4, 3.4, 4.0 };

// ------------------------boyu zhou--------------------------------------
GradTrajOptimizer *_gtop;
ros::Publisher _sdf_pub;
sdf_tools::SignedDistanceField _sdf;
sdf_tools ::CollisionMapGrid _collision_map;
// ------------------------boyu zhou--------------------------------------

void pubCmd()
{
    _cmd.kx[0] = pos_gain[0];
    _cmd.kx[1] = pos_gain[1];
    _cmd.kx[2] = pos_gain[2];

    _cmd.kv[0] = vel_gain[0];
    _cmd.kv[1] = vel_gain[1];
    _cmd.kv[2] = vel_gain[2];

    _cmd.trajectory_id = 1;

    _cmd.yaw_dot = 0.0;
    _cmd.yaw = 0.0;

    if(traj_finish == false)
    {
        return;
        //_cmd.position   = _odom.pose.pose.position;
        _cmd.header.stamp = _odom.header.stamp;
        _cmd.header.frame_id = "/map";
        _cmd.trajectory_flag = 0;

        _cmd.position.x = 0.0;
        _cmd.position.y = 0.0;
        _cmd.position.z = 0.0;

        _cmd.velocity.x = 0.0;
        _cmd.velocity.y = 0.0;
        _cmd.velocity.z = 0.0;

        _cmd.acceleration.x = 0.0;
        _cmd.acceleration.y = 0.0;
        _cmd.acceleration.z = 0.0;
        _cmd_pub.publish(_cmd);
    }
    else if(odom_time > traj_time_final)
    {
        _cmd.header.stamp = _odom.header.stamp;
        _cmd.header.frame_id = "/map";
        _cmd.trajectory_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_COMPLETED;
        ;
        _cmd.position = _odom.pose.pose.position;
        _cmd_pub.publish(_cmd);
    }
    else
    {
        // publish position, velocity and acceleration command according to time bias
        double t = (odom_time - traj_time_start).toSec();

        // ROS_WARN("[Time Optimal Trajectory Node] publish command, time is %f", t);
        MatrixXd time = time_allocator->time;
        MatrixXd time_acc = time_allocator->time_acc;
        // cout<<time.row(0)<<endl;
        int grid_num = time.cols();

        int idx;
        for(idx = 0; idx < _segment_num - 1; idx++)
            if(t > time(idx, grid_num - 1))
                t -= time(idx, grid_num - 1);
            else
                break;

        double t_tmp = t;
        // ROS_WARN("[Time Optimal Trajectory Node] publish command, segm index is %d, segm time is
        // %f", idx, t);
        // now we need to find which grid the time instance belongs to
        int grid_idx;
        for(grid_idx = 0; grid_idx < grid_num; grid_idx++)
        {
            if(t > time(idx, grid_idx))
                continue;
            else
            {
                // cout<<"grid_idx: "<<grid_idx<<", time(idx, grid_idx): "<<time(idx,
                // grid_idx)<<endl;

                if(grid_idx > 0)
                    t -= time(idx, grid_idx - 1);
                else
                    t -= 0.0;

                break;
            }
        }

        // ROS_WARN("[Time Optimal Trajectory Node] publish command, grid index is %d, grid time is
        // %f", grid_idx, t);
        double delta_t;
        if(grid_idx > 0)
            delta_t = (time(idx, grid_idx) - time(idx, grid_idx - 1));
        else
            delta_t = time(idx, grid_idx) - 0.0;
        // seemed has BUG ? shoudl be grid_idx - 1 ?

        double delta_s = t * time_allocator->s_step / delta_t;
        double s = time_allocator->s(grid_idx) + delta_s;

        Vector3d position_s = getPosPoly(idx, s);
        Vector3d position = position_s;

        double s_k = time_allocator->s(grid_idx);
        double s_k_1 = time_allocator->s(grid_idx + 1);
        double b_k = time_allocator->b(idx, grid_idx);
        double b_k_1 = time_allocator->b(idx, grid_idx + 1);

        Vector3d velocity_s1 = getVelPoly(idx, s_k);
        Vector3d velocity_s2 = getVelPoly(idx, s_k_1);

        Vector3d velocity1 = velocity_s1 * sqrt(b_k);
        Vector3d velocity2 = velocity_s2 * sqrt(b_k_1);

        Vector3d velocity = velocity1 + (velocity2 - velocity1) * t / delta_t;

        // reset grid_idx and t for acc time
        t = t_tmp;
        for(grid_idx = 0; grid_idx < grid_num; grid_idx++)
        {
            if(t > time_acc(idx, grid_idx))
                continue;
            else
            {
                // cout<<"grid_idx: "<<grid_idx<<", time(idx, grid_idx): "<<time(idx,
                // grid_idx)<<endl;
                if(grid_idx > 0)
                    t -= time_acc(idx, grid_idx - 1);
                else
                    t -= 0.0;

                break;
            }
        }

        if(grid_idx == grid_num) t -= time_acc(idx, grid_num - 1);

        // ROS_WARN("[Time Optimal Trajectory Node] publish command, grid index is %d, grid time is
        // %f", grid_idx, t);
        Vector3d velocity_s, acceleration_s, acceleration1, acceleration2;
        Vector3d acceleration;

        double a_k;
        if(grid_idx == 0 && idx == 0)
        {
            s_k = time_allocator->s(grid_idx);
            s_k_1 = time_allocator->s(grid_idx + 1);

            a_k = time_allocator->a(idx, grid_idx);
            b_k = time_allocator->b(idx, grid_idx);
            b_k_1 = time_allocator->b(idx, grid_idx + 1);

            velocity_s = getVelPoly(idx, (s_k + s_k_1) / 2.0);
            acceleration_s = getAccPoly(idx, (s_k + s_k_1) / 2.0);
            acceleration = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;
            acceleration << 0.0, 0.0, 0.0;
        }
        else if(grid_idx == grid_num && idx == (_segment_num - 1))
        {
            s_k = time_allocator->s(grid_num - 1);
            s_k_1 = time_allocator->s(grid_num);

            a_k = time_allocator->a(idx, grid_num - 1);
            b_k = time_allocator->b(idx, grid_num - 1);
            b_k_1 = time_allocator->b(idx, grid_num);

            velocity_s = getVelPoly(idx, (s_k + s_k_1) / 2.0);
            acceleration_s = getAccPoly(idx, (s_k + s_k_1) / 2.0);
            acceleration = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;
            acceleration << 0.0, 0.0, 0.0;
        }
        else
        {
            if(grid_idx < grid_num && grid_idx > 0)  // take average accleration in a same segment
            {
                delta_t = (time_acc(idx, grid_idx) - time_acc(idx, grid_idx - 1));

                s_k = time_allocator->s(grid_idx - 1);
                s_k_1 = time_allocator->s(grid_idx + 0);

                a_k = time_allocator->a(idx, grid_idx - 1);
                b_k = time_allocator->b(idx, grid_idx - 1);
                b_k_1 = time_allocator->b(idx, grid_idx + 0);

                velocity_s = getVelPoly(idx, (s_k + s_k_1) / 2.0);
                acceleration_s = getAccPoly(idx, (s_k + s_k_1) / 2.0);
                acceleration1 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;

                // acceleration = acceleration1;

                s_k = time_allocator->s(grid_idx + 0);
                s_k_1 = time_allocator->s(grid_idx + 1);

                a_k = time_allocator->a(idx, grid_idx + 0);
                b_k = time_allocator->b(idx, grid_idx + 0);
                b_k_1 = time_allocator->b(idx, grid_idx + 1);
                // cout<<"a_k: "<<a_k<<" , "<<"s_k: "<<s_k<<" , "<<"s_k+1: "<<s_k_1<<endl;

                velocity_s = getVelPoly(idx, (s_k + s_k_1) / 2.0);
                acceleration_s = getAccPoly(idx, (s_k + s_k_1) / 2.0);
                acceleration2 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;
                acceleration = acceleration1 + (acceleration2 - acceleration1) * t / delta_t;
            }
            else if(grid_idx == grid_num)  // take average accleration between two segments
            {
                delta_t =
                    (time(idx, grid_num - 1) - time_acc(idx, grid_num - 1) + time_acc(idx + 1, 0));

                s_k = time_allocator->s(grid_idx - 1);
                s_k_1 = time_allocator->s(grid_idx);

                a_k = time_allocator->a(idx, grid_idx - 1);
                b_k = time_allocator->b(idx, grid_idx - 1);
                b_k_1 = time_allocator->b(idx, grid_idx);

                velocity_s = getVelPoly(idx, (s_k + s_k_1) / 2.0);
                acceleration_s = getAccPoly(idx, (s_k + s_k_1) / 2.0);
                acceleration1 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;
                s_k = time_allocator->s(0);
                s_k_1 = time_allocator->s(1);

                a_k = time_allocator->a(idx + 1, 0);
                b_k = time_allocator->b(idx + 1, 0);
                b_k_1 = time_allocator->b(idx + 1, 1);

                velocity_s = getVelPoly(idx + 1, (s_k + s_k_1) / 2.0);
                acceleration_s = getAccPoly(idx + 1, (s_k + s_k_1) / 2.0);
                acceleration2 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;
                acceleration = acceleration1 + (acceleration2 - acceleration1) * t / delta_t;
            }
            else if(grid_idx == 0)  // take average accleration between two segments
            {
                delta_t = (time(idx - 1, grid_num - 1) - time_acc(idx - 1, grid_num - 1) +
                           time_acc(idx, 0));

                s_k = time_allocator->s(grid_num - 1);
                s_k_1 = time_allocator->s(grid_num);

                a_k = time_allocator->a(idx - 1, grid_num - 1);
                b_k = time_allocator->b(idx - 1, grid_num - 1);
                b_k_1 = time_allocator->b(idx - 1, grid_num);

                velocity_s = getVelPoly(idx - 1, (s_k + s_k_1) / 2.0);
                acceleration_s = getAccPoly(idx - 1, (s_k + s_k_1) / 2.0);
                acceleration1 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;

                s_k = time_allocator->s(grid_idx);
                s_k_1 = time_allocator->s(grid_idx + 1);

                a_k = time_allocator->a(idx, grid_idx);
                b_k = time_allocator->b(idx, grid_idx);
                b_k_1 = time_allocator->b(idx, grid_idx + 1);

                velocity_s = getVelPoly(idx, (s_k + s_k_1) / 2.0);
                acceleration_s = getAccPoly(idx, (s_k + s_k_1) / 2.0);
                acceleration2 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;
                acceleration = acceleration1 + (acceleration2 - acceleration1) *
                                                   (t + time(idx - 1, grid_num - 1) -
                                                    time_acc(idx - 1, grid_num - 1)) /
                                                   delta_t;
            }
            else
                ROS_BREAK();
        }

        _cmd.header.stamp = odom_time;
        _cmd.header.frame_id = "/map";
        _cmd.trajectory_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_READY;
        ;

        _cmd.position.x = position(0);
        _cmd.position.y = position(1);
        _cmd.position.z = position(2);

        _cmd.velocity.x = velocity(0);
        _cmd.velocity.y = velocity(1);
        _cmd.velocity.z = velocity(2);

        _cmd.acceleration.x = acceleration(0);
        _cmd.acceleration.y = acceleration(1);
        _cmd.acceleration.z = acceleration(2);

        _cmd_pub.publish(_cmd);

        _vis_pos.header.stamp = odom_time;
        _vis_vel.header.stamp = odom_time;
        _vis_acc.header.stamp = odom_time;
        _vis_pos.points.clear();
        _vis_vel.points.clear();
        _vis_acc.points.clear();

        geometry_msgs::Point pt;
        pt.x = position(0);
        pt.y = position(1);
        pt.z = position(2);

        _vis_pos.points.push_back(pt);
        _vis_vel.points.push_back(pt);
        _vis_acc.points.push_back(pt);

        double base_l = 3.0;
        double length;

        length = sqrt(pow(velocity(0), 2.0) + pow(velocity(1), 2.0) + pow(velocity(2), 2.0)) /
                 _MAX_Vel * base_l;
        pt.x = position(0) + length * velocity(0);
        pt.y = position(1) + length * velocity(1);
        pt.z = position(2) + length * velocity(2);

        _vis_vel.points.push_back(pt);

        length = sqrt(pow(acceleration(0), 2.0) + pow(acceleration(1), 2.0) +
                      pow(acceleration(2), 2.0)) /
                 _MAX_Acc * base_l;
        pt.x = position(0) + length * acceleration(0);
        pt.y = position(1) + length * acceleration(1);
        pt.z = position(2) + length * acceleration(2);

        _vis_acc.points.push_back(pt);

        _vis_pos_pub.publish(_vis_pos);
        _vis_vel_pub.publish(_vis_vel);
        _vis_acc_pub.publish(_vis_acc);
    }
}

void rcvWaypointsCallback(const nav_msgs::Path &wp)
{
    if(wp.poses[0].pose.position.z < 0.0) return;

    traj_finish = false;
    endPt.x = wp.poses[0].pose.position.x;
    endPt.y = wp.poses[0].pose.position.y;
    endPt.z = wp.poses[0].pose.position.z;

    _pathPlaner.reset();
    _pathPlaner.setPt(startPt, endPt);
    _pathPlaner.RRGpathFind();

    if(!_pathPlaner.path_find_state)
    {
        ROS_WARN("In planning node, can't find a path");
        return;
    }

    _Path = _pathPlaner.getPath();
    _Radius = _pathPlaner.getRadius();

    visPath(_Path, _Radius);
    cout << "Path:" << _Path << endl;
    cout << "radius:" << _Radius << endl;
    ros::Duration(1.0).sleep();
    trajGenerationTest(_Path, _Radius);
}

double _scale;
void trajGenerationTest(Eigen::MatrixXd path, Eigen::VectorXd radius)
{
    Eigen::MatrixXd pos(2, 3);
    Eigen::MatrixXd vel(2, 3);
    Eigen::MatrixXd acc(2, 3);

    pos << startPt.x, startPt.y, startPt.z, endPt.x, endPt.y, endPt.z;

    vel << startVel.x, startVel.y, startVel.z, 0.0, 0.0, 0.0;

    acc << startAcc.x, startAcc.y, startAcc.z, 0.0, 0.0, 0.0;

    Eigen::MatrixXd path_ = path;
    Eigen::VectorXd radius_ = radius;

    _segment_num = radius_.size();
    _scale = 1.0;

    _Time.resize(_segment_num);
    for(int i = 0; i < _Time.size(); i++)
        _Time(i) = 1.0;

    // ----------------generate trajectory in virtual domain:s:[0,1]--------------------
    vector<Eigen::MatrixXd> polyCoeff_pList;
    polyCoeff_pList =  _trajectoryGeneratorPoly.PloyCoeffGeneration(path_, radius_, _Time, vel,
    acc, _MAX_Vel, _MAX_Acc, 8, _scale, 1.0, 10.0); 
    _PolyCoeff_p = polyCoeff_pList.back();
    _poly_num1D = _PolyCoeff_p.cols() / 3;

    // ------------boyu zhou----------------------
    cout << setprecision(2) << "Poly by TrajectoryGenerator:" << _PolyCoeff_p << endl;

    // ------------use gradient based planner---------------------

    // _poly_num1D = 6;
    // convert path to waypoints format
    std::vector<Eigen::Vector3d> waypoints;
    for(int i = 0; i < path_.rows(); ++i)
    {
        Eigen::Vector3d wp = path_.row(i);
        waypoints.push_back(wp);
    }
    _gtop->setWaypoints(waypoints);

    Eigen::MatrixXd coeff;
    // optimization in time domain
    _gtop->setDomainType(1);
    _gtop->initialize();
    _gtop->optimizeTrajectory(1);
    _gtop->optimizeTrajectory(2);
    _gtop->getCoefficient(coeff);
    polyCoeff_pList.push_back(coeff);

    cout << "Poly by gradtrajoptimizer" << coeff << endl;
    // _PolyCoeff_p = coeff;

    // // optimization in virtual domain
    // _gtop->setDomainType(1);
    // _gtop->initialize();
    // _gtop->optimizeTrajectory(1);
    // _gtop->optimizeTrajectory(2);
    // _gtop->getCoefficient(coeff);
    // polyCoeff_pList.push_back(coeff);

    // ------------boyu zhou----------------------

    if(polyCoeff_pList.back().rows() == 3 && polyCoeff_pList.back().cols() == 3)
        ROS_WARN("Monomial solver failed");
    else
    {
        visPolyTraj(polyCoeff_pList);
        visCheckPt(polyCoeff_pList);
    }

    int K = 30;

    MinimumTimeOptimizer time_optimizer;

    // ##### Here only for comparing the difference when using different Q ###
    // Newly add a coefficient for regularizing the control effort when minizing the total time;
    // re-write part of the solver; consistent with the "Optimal Control" concept.
    // How smart am I !

    time_optimizer.MinimumTimeGeneration(polyCoeff_pList.back(), 2.0, 1.0, 0.1, K, 0.0);
    time_optimizer.MinimumTimeGeneration(polyCoeff_pList.back(), 2.0, 1.0, 0.1, K, 10.0);

    time_allocator = time_optimizer.GetTimeAllcoation();

    // traj_time_start = odom_time;
    traj_time_start = ros::Time::now();
    traj_time_final = traj_time_start;

    //    ROS_WARN("[Time Optimal Trajectory Node] check all time allocated in each segmetn");
    for(int i = 0; i < time_allocator->time.rows(); i++)
    {
        // cout<<time_allocator->time(i, K - 1)<<endl;
        traj_time_final += ros::Duration(time_allocator->time(i, K - 1));
    }

    plotData(time_allocator);

    traj_finish = true;

    // for(int k = 0; k < _segment_num; k++)
    // {
    //     cout << "time in velocity discretization \n" << time_allocator->time.row(k) << endl;
    //     cout << "time in acceleration discretization \n" << time_allocator->time_acc.row(k) << endl;
    // }
    // ROS_BREAK();
}

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map)
{
    pcl::PointCloud<pcl::PointXYZ> CloudIn;
    pcl::fromROSMsg(pointcloud_map, CloudIn);
    _pathPlaner.setInput(CloudIn);

    // ----------- boyu zhou------------------------

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
    _collision_map = sdf_tools ::CollisionMapGrid(origin_transform, frame, resolution, x_size,
                                                  y_size, z_size, oob_cell);
    // sdf_tools ::CollisionMapGrid collision_map(origin_transform, frame, resolution, x_size,
    // y_size,
    //                                            z_size, oob_cell);

    // add point cloud to sdf map
    sdf_tools::COLLISION_CELL obstacle_cell(1.0);
    std::vector<pcl::PointXYZ, Eigen::aligned_allocator<pcl::PointXYZ>> pts =
        cloud_filtered->points;
    for(int i = 0; i < int(pts.size()); ++i)
    {
        if((pts[i].x < 14.9 && pts[i].x > -14.9) && (pts[i].y < 14.9 && pts[i].y > -14.9) &&
           (pts[i].z < 4.5 && pts[i].z > 0.2))
        {
            pcl::PointXYZ addpt = pts[i];
            double add_x = addpt.x;
            double add_y = addpt.y;
            double add_z = addpt.z;
            
            _collision_map.Set(add_x, add_y, add_z, obstacle_cell);
            // std::cout << setprecision(2) << i << " th pt x:" << add_x << " ," << add_y << " ,"
            //           << add_z << std::endl;
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
        _collision_map.ExportForDisplay(collision_color, free_color, unknown_color);
    collision_map_marker.ns = "collision_map";
    collision_map_marker.id = 1;

    _sdf_pub.publish(collision_map_marker);

    // Build the signed distance field
    float oob_value = INFINITY;
    std::pair<sdf_tools::SignedDistanceField, std::pair<double, double>> sdf_with_extrema =
        _collision_map.ExtractSignedDistanceField(oob_value);
    _sdf = sdf_with_extrema.first;

    // set sdf to gradient planner
    _gtop->setSignedDistanceField(&_sdf, 0.2);
    // _gtop->setCollisionMap(&_collision_map, 0.2);

    // ----------- boyu zhou------------------------
}

void rcvOdometryCallbck(const nav_msgs::Odometry odom)
{
    if(odom.child_frame_id == "X" || odom.child_frame_id == "O") return;
    _odom = odom;
    _has_odom = true;

    odom_time = _odom.header.stamp;

    _state = { _odom.pose.pose.position.x,
               _odom.pose.pose.position.y,
               _odom.pose.pose.position.z,
               _odom.twist.twist.linear.x,
               _odom.twist.twist.linear.y,
               _odom.twist.twist.linear.z,
               0.0,
               0.0,
               0.0 };

    startPt.x = _odom.pose.pose.position.x;
    startPt.y = _odom.pose.pose.position.y;
    startPt.z = _odom.pose.pose.position.z;

    startVel.x = _odom.twist.twist.linear.x;
    startVel.y = _odom.twist.twist.linear.y;
    startVel.z = _odom.twist.twist.linear.z;
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "timeOptialTraj");
    ros::NodeHandle n("~");

    n.param("dynamic/vec", _Vel, 2.0);
    n.param("dynamic/acc", _Acc, 1.0);
    n.param("dynamic/max_vec", _MAX_Vel, 3.0);
    n.param("dynamic/max_acc", _MAX_Acc, 1.5);
    n.param("optimization/poly_order_min", _poly_order_min, 5);
    n.param("optimization/poly_order_max", _poly_order_max, 10);
    n.param("optimization/minimize_order", _minimize_order, 3);

    Bernstein _bernstein;
    if(_bernstein.setParam(_poly_order_min, _poly_order_max, _minimize_order) == -1)
    {
        ROS_ERROR(" The trajectory order is set beyond the library's scope, please re-set ");
    }

    _MQMList = _bernstein.getMQM();
    _FMList = _bernstein.getFM();
    _CList = _bernstein.getC();
    _CvList = _bernstein.getC_v();
    _CaList = _bernstein.getC_a();

    // subcribed msgs
    _dest_pts_sub = n.subscribe("waypoints", 1, rcvWaypointsCallback);
    _map_sub = n.subscribe("PointCloud", 1, rcvPointCloudCallBack);
    _odom_sub = n.subscribe("odometry", 50, rcvOdometryCallbck);

    // publish visualize related msgs
    _ctrlPt_vis_pub = n.advertise<visualization_msgs::Marker>("ctrlPt", 1);
    _traj_vis_pub = n.advertise<visualization_msgs::Marker>("trajectory_vis", 1);
    _traj_bezier_vis_pub = n.advertise<visualization_msgs::Marker>("bezier_trajectory_vis", 1);
    _checkPt_vis_pub = n.advertise<visualization_msgs::MarkerArray>("checkpoint_vis", 1);

    _poly_traj_vis_pub = n.advertise<visualization_msgs::Marker>("monomial_trajectory_vis", 1);
    _checkTraj_vis_pub = n.advertise<visualization_msgs::Marker>("check_trajectory_vis", 1);
    _path_vis_pub = n.advertise<visualization_msgs::MarkerArray>("path", 1);
    _traj_pub = n.advertise<quadrotor_msgs::PolynomialTrajectory>("trajectory", 10);

    // for visualize the minimal time control command
    _vis_pos_pub = n.advertise<visualization_msgs::Marker>("desired_position", 50);
    _vis_vel_pub = n.advertise<visualization_msgs::Marker>("desired_velocity", 50);
    _vis_acc_pub = n.advertise<visualization_msgs::Marker>("desired_acceleration", 50);

    _cmd_pub = n.advertise<quadrotor_msgs::PositionCommand>("position_command", 50);

    _vis_pos.ns = "pos";
    _vis_pos.id = 0;
    _vis_pos.header.frame_id = "/map";
    _vis_pos.type = visualization_msgs::Marker::SPHERE;
    _vis_pos.action = visualization_msgs::Marker::ADD;
    _vis_pos.color.a = 1.0;
    _vis_pos.color.r = 0.0;
    _vis_pos.color.g = 0.0;
    _vis_pos.color.b = 0.0;
    _vis_pos.scale.x = 0.2;
    _vis_pos.scale.y = 0.2;
    _vis_pos.scale.z = 0.2;

    _vis_vel.ns = "vel";
    _vis_vel.id = 0;
    _vis_vel.header.frame_id = "/map";
    _vis_vel.type = visualization_msgs::Marker::ARROW;
    _vis_vel.action = visualization_msgs::Marker::ADD;
    _vis_vel.color.a = 1.0;
    _vis_vel.color.r = 0.0;
    _vis_vel.color.g = 1.0;
    _vis_vel.color.b = 0.0;
    _vis_vel.scale.x = 0.1;
    _vis_vel.scale.y = 0.2;
    _vis_vel.scale.z = 0.2;

    _vis_acc.ns = "acc";
    _vis_acc.id = 0;
    _vis_acc.header.frame_id = "/map";
    _vis_acc.type = visualization_msgs::Marker::ARROW;
    _vis_acc.action = visualization_msgs::Marker::ADD;
    _vis_acc.color.a = 1.0;
    _vis_acc.color.r = 1.0;
    _vis_acc.color.g = 1.0;
    _vis_acc.color.b = 0.0;
    _vis_acc.scale.x = 0.1;
    _vis_acc.scale.y = 0.2;
    _vis_acc.scale.z = 0.2;

    // ----------- boyu zhou test------------------------
    _gtop = new GradTrajOptimizer(n);
    _sdf_pub = n.advertise<visualization_msgs::Marker>("sdf_visualization", 1, true);
    // ----------- boyu zhou test------------------------

    ros::Rate rate(100);
    bool status = ros::ok();
    while(status)
    {
        pubCmd();
        ros::spinOnce();
        status = ros::ok();
        rate.sleep();
    }

    return 0;
}

vector<double> getFullStateFromBezier(const Eigen::MatrixXd &polyCoeff, double t_now, int seg_now)
{
    vector<double> ret(9, 0);
    Eigen::VectorXd ctrl_now = polyCoeff.row(seg_now);

    int order = _poly_orderList[seg_now];
    int ctrl_num1D = order + 1;

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < ctrl_num1D; j++)
        {
            ret[i] += _CList[order](j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) *
                      pow((1 - t_now), (order - j));

            if(j < ctrl_num1D - 1)
                ret[i + 3] += _CvList[order](j) * order *
                              (ctrl_now(i * ctrl_num1D + j + 1) - ctrl_now(i * ctrl_num1D + j)) *
                              pow(t_now, j) * pow((1 - t_now), (order - j - 1));

            if(j < ctrl_num1D - 2)
                ret[i + 6] +=
                    _CaList[order](j) * order * (order - 1) *
                    (ctrl_now(i * ctrl_num1D + j + 2) - 2 * ctrl_now(i * ctrl_num1D + j + 1) +
                     ctrl_now(i * ctrl_num1D + j)) *
                    pow(t_now, j) * pow((1 - t_now), (order - j - 2));
        }
    }
    return ret;
}

vector<double> getStateFromBezier(const Eigen::MatrixXd &polyCoeff, double t_now, int seg_now)
{
    vector<double> ret(3, 0);
    Eigen::VectorXd ctrl_now = polyCoeff.row(seg_now);

    int order = _poly_orderList[seg_now];
    int ctrl_num1D = order + 1;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < ctrl_num1D; j++)
        {
            ret[i] += _CList[order](j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) *
                      pow((1 - t_now), (order - j));
        }
    }

    // ROS_WARN("finish get state from Bezier curve ...");
    return ret;
}

vector<double> getStateFromPolynomial(const Eigen::MatrixXd &PolyCoeff, const Eigen::VectorXd &Time,
                                      double t_now, int seg_now)
{
    vector<double> ret(3 * 3, 0);

    t_now = min(max(0.0, t_now), Time(seg_now));
    int poly_num1D = PolyCoeff.cols() / 3;

    for(int dim = 0; dim < 3; dim++)
    {
        Eigen::VectorXd coeff = (PolyCoeff.row(seg_now)).segment(dim * poly_num1D, poly_num1D);

        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(3, poly_num1D);

        for(int j = 0; j < poly_num1D; j++)
        {
            t(0, j) = pow(t_now, j);
            t(1, j) = j * pow(t_now, j - 1);
            t(2, j) = j * (j - 1) * pow(t_now, j - 2);
        }

        for(int i = 0; i < 3; i++)
            ret[dim + i * 3] = coeff.dot(t.row(i));
    }

    return ret;
}

vector<double> getFullStateFromPolynomial(const Eigen::MatrixXd &PolyCoeff,
                                          const Eigen::VectorXd &Time, double t_now, int seg_now)
{
    vector<double> ret(3 * 3, 0);

    t_now = min(max(0.0, t_now), Time(seg_now));
    int poly_num1D = PolyCoeff.cols() / 3;

    for(int dim = 0; dim < 3; dim++)
    {
        Eigen::VectorXd coeff = (PolyCoeff.row(seg_now)).segment(dim * poly_num1D, poly_num1D);

        Eigen::MatrixXd t = Eigen::MatrixXd::Zero(3, poly_num1D);

        for(int j = 0; j < poly_num1D; j++)
        {
            t(0, j) = pow(t_now, j);
            t(1, j) = j * pow(t_now, j - 1);
            t(2, j) = j * (j - 1) * pow(t_now, j - 2);
        }

        for(int i = 0; i < 3; i++)
            ret[dim + i * 3] = coeff.dot(t.row(i));
    }

    return ret;
}

void visPolyTraj(vector<Eigen::MatrixXd> polyCoeffList)
{
    for(int k = 0; k < (int)polyCoeffList.size(); k++)
    {
        visualization_msgs::Marker _traj_vis;

        _traj_vis.header.stamp = ros::Time::now();
        _traj_vis.header.frame_id = "map";

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
        Eigen::Vector3d cur, pre;
        cur.setZero();
        pre.setZero();

        _traj_vis.points.clear();
        vector<double> state;
        geometry_msgs::Point pt;
        Vector3d max_v(0.0, 0.0, 0.0);
        Vector3d max_a(0.0, 0.0, 0.0);

        for(int i = 0; i < _Time.size(); i++)
        {
            for(double t = 0.0; t < _Time(i); t += 0.02 / _scale, count += 1)
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
                        if(v_p > max_v(p)) max_v(p) = v_p;

                        double a_p = fabs(state[6 + p]) / _scale;
                        if(a_p > max_a(p)) max_a(p) = a_p;
                    }

                if(count) traj_len += (pre - cur).norm();
                pre = cur;
            }
        }
        /*ROS_WARN(" max velocity of monomial curve in x, y, z axis are %f, %f, %f ",     max_v(0),
        max_v(1), max_v(2)); ROS_WARN(" max acceleration of monomial curve in x, y, z axis are %f,
        %f, %f ", max_a(0), max_a(1), max_a(2));*/
        ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);
        _poly_traj_vis_pub.publish(_traj_vis);
    }
}

void visBezierTrajectory(Eigen::MatrixXd polyCoeff, int rgb)
{
    visualization_msgs::Marker _traj_vis;

    _traj_vis.header.stamp = ros::Time::now();
    _traj_vis.header.frame_id = "map";

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
    Eigen::Vector3d cur, pre;
    cur.setZero();
    pre.setZero();

    _traj_vis.points.clear();

    vector<double> state;
    geometry_msgs::Point pt;
    Vector3d max_v(0.0, 0.0, 0.0);
    Vector3d max_a(0.0, 0.0, 0.0);

    for(int i = 0; i < _segment_num; i++)
    {
        for(double t = 0.0; t < 1.0; t += 0.01 / _Time(i), count += 1)
        {
            state = getFullStateFromBezier(polyCoeff, t, i);
            cur(0) = pt.x = _Time(i) * state[0];
            cur(1) = pt.y = _Time(i) * state[1];
            cur(2) = pt.z = _Time(i) * state[2];
            _traj_vis.points.push_back(pt);

            if(count) traj_len += (pre - cur).norm();
            pre = cur;

            for(int p = 0; p < 3; p++)
            {
                double v_ = fabs(state[3 + p]);
                double a_ = fabs(state[6 + p]) / _Time(i);

                if(v_ > max_v(p)) max_v(p) = v_;

                if(a_ > max_a(p)) max_a(p) = a_;
            }
        }
    }

    ROS_WARN(" max velocity of the Bezier curve in x, y, z axis are %f, %f, %f ", max_v(0),
             max_v(1), max_v(2));
    ROS_WARN(" max acceleration of the Bezier curve in x, y, z axis are %f, %f, %f ", max_a(0),
             max_a(1), max_a(2));
    ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);

    _traj_bezier_vis_pub.publish(_traj_vis);
}

void visPath(Eigen::MatrixXd path, Eigen::VectorXd radius)
{
    Eigen::MatrixXd path_new = path.block(1, 0, path.rows() - 2, path.cols());

    for(auto &mk : path_vis.markers)
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

    for(int i = 0; i < int(path_new.rows()); i++)
    {
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

void visCtrlPoint(Eigen::MatrixXd polyCoeff, int rgb)
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
}

visualization_msgs::MarkerArray checkPt_vis;
void visCheckPt(vector<Eigen::MatrixXd> polyCoeffList)
{
    MatrixXd poly = polyCoeffList.back();

    for(auto &mk : checkPt_vis.markers)
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

    for(int i = 0; i < _Time.size(); i++)
    {
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
}

Vector3d getPosPoly(int k, double s)
{
    Vector3d ret;

    for(int dim = 0; dim < 3; dim++)
    {
        VectorXd coeff = (_PolyCoeff_p.row(k)).segment(dim * _poly_num1D, _poly_num1D);
        VectorXd t = VectorXd::Zero(_poly_num1D);

        for(int j = 0; j < _poly_num1D; j++)
            if(j == 0)
                t(j) = 1.0;
            else
                t(j) = pow(s, j);

        ret(dim) = coeff.dot(t);
    }

    return ret;
}

Vector3d getVelPoly(int k, double s)
{
    Vector3d ret;

    for(int dim = 0; dim < 3; dim++)
    {
        VectorXd coeff = (_PolyCoeff_p.row(k)).segment(dim * _poly_num1D, _poly_num1D);
        VectorXd t = VectorXd::Zero(_poly_num1D);

        for(int j = 0; j < _poly_num1D; j++)
            if(j == 0)
                t(j) = 0.0;
            else
                t(j) = j * pow(s, j - 1);

        ret(dim) = coeff.dot(t);
    }

    return ret;
}

Vector3d getAccPoly(int k, double s)
{
    Vector3d ret;

    for(int dim = 0; dim < 3; dim++)
    {
        VectorXd coeff = (_PolyCoeff_p.row(k)).segment(dim * _poly_num1D, _poly_num1D);
        VectorXd t = VectorXd::Zero(_poly_num1D);

        for(int j = 0; j < _poly_num1D; j++)
            if(j == 0 || j == 1)
                t(j) = 0.0;
            else
                t(j) = j * (j - 1) * pow(s, j - 2);

        ret(dim) = coeff.dot(t);
    }

    return ret;
}

void plotData(Allocator *time_allocator)
{
}