#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/transforms.h>
#include <pcl/registration/default_convergence_criteria.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/point_cloud.h>
#include <pcl/pcl_conversions.h>
#include <pcl/voxel_grid.h>
#include <pcl/statistical_outlier_removal.h>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/Range.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <fstream>

#include <quadrotor_msgs/Odometry.h>

#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <math.h>

using namespace Eigen;
using namespace std;

int _itr_cycle;

ros::Time last_stamp;
ros::Publisher _pub;
ros::Publisher _pub_ref;
ros::Publisher _pub_target;
//ros::Publisher _pub_match;
ros::Publisher _pub_odom;
ros::Publisher _pub_rawdata;
ros::Publisher _pub_laser_path;
ros::Publisher _pub_input_aftrot;
ros::Publisher _pub_enviro;

ros::Publisher _pub_laser_quad_odom;
ros::Publisher _pub_map;

pcl::PointCloud<pcl::PointXYZ> cloud_target;
pcl::PointCloud<pcl::PointXYZ> cloud_input;
//pcl::PointCloud<pcl::PointXYZ> cloud_input_last;


sensor_msgs::PointCloud2 pointcloud_ref;
sensor_msgs::PointCloud2 pointcloud_env;
bool is_laser_init = true;
bool is_imu_init = true;
bool is_change_k = false;
Matrix3d R_imu_init_inv   = MatrixXd::Identity(3, 3);

double _key_thre;
//double roll, pitch;
double fitness_thre;
double _corres_dis;
//Matrix4d RT  = MatrixXd::Identity(4, 4);
Matrix4d RT_init = MatrixXd::Identity(4, 4);
Matrix4d RT_icp  = MatrixXd::Identity(4, 4);
Matrix4d RT_kframe  = MatrixXd::Identity(4, 4);
//Matrix4d RT_odom_last = MatrixXd::Identity(4, 4);
Matrix3d _R_imu;
    Matrix3d R_imu;
double height = 0.0;

nav_msgs::Odometry laser_odom;
nav_msgs::Path laser_path;
//vector<sensor_msgs::Imu> Imu_list;

Eigen::Vector3d R_to_rpy(const Eigen::Matrix3d& R)
{
  Eigen::Vector3d rpy;
  double r = asin(R(2,1));
  double p = atan2(-R(2,0)/cos(r), R(2,2)/cos(r));
  double y = atan2(-R(0,1)/cos(r), R(1,1)/cos(r));
  rpy(0) = r;
  rpy(1) = p;
  rpy(2) = y;

  return rpy;
}

Vector3d R_to_ypr(const Matrix3d& R)
{
    Vector3d n = R.col(0);
    Vector3d o = R.col(1);
    Vector3d a = R.col(2);

    Vector3d ypr(3);
    double y = atan2(n(1), n(0));
    double p = atan2(-n(2), n(0)*cos(y)+n(1)*sin(y));
    double r = atan2(a(0)*sin(y)-a(1)*cos(y), -o(0)*sin(y)+o(1)*cos(y));
    ypr(0) = y;
    ypr(1) = p;
    ypr(2) = r;

    return ypr;
}

Matrix3d ypr_to_R(const Vector3d& ypr)
{
    double c, s;
    Matrix3d Rz = Matrix3d::Zero();
    double y = ypr(0);
    c = cos(y);
    s = sin(y);
    Rz(0,0) =  c;
    Rz(1,0) =  s;
    Rz(0,1) = -s;
    Rz(1,1) =  c;
    Rz(2,2) =  1;

    Matrix3d Ry = Matrix3d::Zero();
    double p = ypr(1);
    c = cos(p);
    s = sin(p);
    Ry(0,0) =  c;
    Ry(2,0) = -s;
    Ry(0,2) =  s;
    Ry(2,2) =  c;
    Ry(1,1) =  1;

    Matrix3d Rx = Matrix3d::Zero();
    double r = ypr(2);
    c = cos(r);
    s = sin(r);
    Rx(1,1) =  c;
    Rx(2,1) =  s;
    Rx(1,2) = -s;
    Rx(2,2) =  c;
    Rx(0,0) =  1;

    Matrix3d R = Rz*Ry*Rx;
    return R;
}


void range_callback(const sensor_msgs::Range &msg)
{   
    if(msg.range > 10)
        return;
    //if(abs(msg.range - height) > 0.25)
    //    return;
    height = msg.range + 0.28;
    
    //if( abs(height - RT_kframe(2,3)) > 0.2 )
     //   RT_kframe(2,3) = height;
}

/*void odom_callback(const nav_msgs::Odometry &msg)
{   
    if(msg.pose.pose.position.z > 10.0)
        return;

    height = msg.pose.pose.position.z + 0.28;
}*/

void imu_callback(const sensor_msgs::Imu::ConstPtr &msg)
{
    Quaterniond Q_imu;
    Q_imu.w() = msg->orientation.w;
    Q_imu.x() = msg->orientation.x;
    Q_imu.y() = msg->orientation.y;
    Q_imu.z() = msg->orientation.z;

    // = Q_imu.toRotationMatrix();
   // Matrix3d R_imu;
    R_imu = Q_imu;

    if(is_imu_init)
    {
        R_imu_init_inv = R_imu.inverse();
        is_imu_init = false;
    }

    _R_imu =  R_imu_init_inv * R_imu;

    //cout<<"imu attitude : \n"<<_R_imu<<endl;
    
    //Matrix3d Riv  = Quaterniond(0, -1, 0, 0).toRotationMatrix();
    Vector3d ypr = R_to_ypr(_R_imu);
    ROS_WARN("Imu Yaw is : ");
    cout<<ypr(0)<<endl;

    //sensor_msgs::Imu imu_msg = *msg;
   // Imu_list.push_back(imu_msg);
}


int counter = 0;
double Z_THRESH = 2.0;
double Y_THRESH = 20.0;
double X_THRESH = 20.0;
double Z_LOW_THRESH = 0.5;
double sum = 0;
int kfid = 0;

Matrix3d R_imu_kf;
//roslaunch laser_ekf.launch itr:=10 key_thre:=2.0 corres_dis:=5.0

//pcl::PointCloud<pcl::PointXYZ> cloud_map;

pair<double, Matrix4d> IcpFunction(pcl::PointCloud<pcl::PointXYZ> cloud_input, pcl::PointCloud<pcl::PointXYZ> cloud_target)
{

        // Perform GICP
        //pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
        pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ,double> icp;

        icp.setInputTarget(cloud_target.makeShared());
        icp.setInputSource(cloud_input.makeShared());

        // Set the max correspondence distance to 5cm (e.g., correspondences with higher distances will be ignored)
        icp.setMaxCorrespondenceDistance (_corres_dis);// value smaller: slower & result worse
        // Set the maximum number of iterations (criterion 1)
        icp.setMaximumIterations(_itr_cycle);
        // Set the transformation epsilon (criterion 2)
        icp.setTransformationEpsilon (1e-8);
        //icp.setRotationEpsilon(1e-6);
        // Set the euclidean distance difference epsilon (criterion 3)
        icp.setEuclideanFitnessEpsilon (0.001);

        pcl::PointCloud<pcl::PointXYZ> match;

        //Vector3d ypr_rp      = R_to_ypr(_R_imu);
        
        //Vector3d ypr_k_frame = R_to_ypr(RT_kframe.block(0,0,3,3));//RT_odom_last.block(0,0,3,3)
        
        //Vector3d ypr_k_frame = R_to_ypr(RT_icp.block(0,0,3,3));//RT_odom_last.block(0,0,3,3)

        //Vector3d ypr;
        //ypr << ypr_k_frame(0),
         //      ypr_rp(1),
         //      ypr_rp(2);

        //Matrix3d R_k_init = ypr_to_R(ypr);        

        RT_init.block(0,0,3,3) = R_imu_kf.inverse() * R_imu;// _R_imu;
        //RT_kframe.block(0,0,3,3).inverse() * R_k_init;//_R_imu;
        RT_init.block(0,3,3,1) = RT_icp.block(0,3,3,1);

        //if(!is_change_k) 
        //{   
        
        //if(abs(RT_kframe(2,3) - height) > 0.2 )
        //    RT_kframe(2,3) = height;
       //if( abs(height - RT_kframe(2,3)) < 0.2 )
            RT_init(2,3) = height - RT_kframe(2,3);
            //RT_init.block(0,0,3,3) = RT_icp.block(0,0,3,3);
        //}
       // cout<<"RT_init is :\n"<<RT_init<<endl;
        //Vector3d rpy_init = R_to_rpy(RT_init.block(0,0,3,3));
        //cout<<"RT init roll pitch yaw is :\n"<<rpy_init<<endl;
        icp.align(match, RT_init);
        ros::Time time1 = ros::Time::now();
        
       //cout<<"has converged: "<<icp.hasConverged()<<" score: "<<icp.getFitnessScore()<<endl;
        // Convert the pcl/PointCloud to sensor_msgs/PointCloud2
        //double fitness = icp.getFitnessScore();
        //icp.computeTransformation(match, RT_init);
        Matrix4d RT_ = icp.getFinalTransformation();
        //cout<<"RT_icp is :\n"<<RT_icp<<endl;
        ros::Time time2 = ros::Time::now();
        return make_pair(icp.getFitnessScore(), RT_);
}

double pc_msg_time;
ros::Time tkf;
Vector3d T_kf;
Quaterniond Q_kf;
bool is_abort = false;
int counter_good_laser = 0;
int counter_map = 0;
Matrix4d RT_odom;
pcl::PointCloud<pcl::PointXYZ> cloud_map;

void pointcloud_callback(const sensor_msgs::PointCloud2& pc_msg) 
{
        //ROS_WARN("[test_ICP] input_pointcloud Received .");
        ros::Time time0 =  ros::Time::now();
        pcl::PointCloud<pcl::PointXYZ> cloud_msg;
        pcl::PointCloud<pcl::PointXYZ> cloud_msg_near;

        pcl::PointCloud<pcl::PointXYZ> cloud_msg_filter;
        pcl::fromROSMsg( pc_msg, cloud_msg );
        cloud_msg.header.frame_id="/map";

        _pub_rawdata.publish(cloud_msg);

        cloud_msg_near.header = cloud_msg.header;
        
        pcl::PointXYZ newpoint;
        pcl::PointCloud<pcl::PointXYZ>::iterator itr;
        for(itr = cloud_msg.begin();itr < cloud_msg.end(); itr++)
        {
            newpoint.x = itr->x;
            newpoint.y = itr->y;
            newpoint.z = itr->z;
            //newpoint.z < Z_THRESH && newpoint.z >  Z_LOW_THRESH  && 
           if ( newpoint.x > -0.1 )//&& abs(newpoint.y) < Y_THRESH)        
                cloud_msg_near.push_back(newpoint);
        }
        
        ros::Time time1=  ros::Time::now();
        //cloud_msg_near = cloud_msg; 
        //&& abs(newpoint.x) < X_THRESH && newpoint.y < Y_THRESH - 6 && newpoint.y > -Y_THRESH

       // ROS_WARN("[test_ICP] Origin pointcloud size is : \n");
       // cout<< cloud_msg.size()<<endl;
        //cout<< cloud_msg_near.size()<<endl;
        
        cloud_msg_filter = cloud_msg_near;        
        pcl::PointCloud<pcl::PointXYZ> cloud_msg_downsampled;
        pcl::VoxelGrid<pcl::PointXYZ>  voxelSampler;

        //cloud_msg_downsampled = cloud_msg_near;
        
        voxelSampler.setInputCloud(cloud_msg_filter.makeShared());
        voxelSampler.setLeafSize(0.02f, 0.02f, 0.02f);
        voxelSampler.filter(cloud_msg_downsampled);
                ros::Time time2=  ros::Time::now();
        
      //  ROS_WARN("ICP cloud size is : ");
     //   cout<<cloud_msg_downsampled.size()<<endl;
        /*
        pcl::PointCloud<pcl::PointXYZ> cloud_msg_downsampled_remove;
        pcl::StatisticalOutlierRemoval<pcl::PointXYZ> statFilter;
        statFilter.setInputCloud(cloud_msg_downsampled.makeShared());
        statFilter.setMeanK(5);
        statFilter.setStddevMulThresh(5.0);
        statFilter.filter(cloud_msg_downsampled_remove);
        */
  /*      
        Matrix4d RT_rot = Matrix4d::Identity(4,4);
        RT_rot.block(0,0,3,3) = _R_imu.inverse();//RT_kframe.block(0,0,3,3).inverse()*
        //RT_rot = RT_rot;
        cout<<"RT_rot : \n"<<RT_rot<<endl;   
        Vector3d rpy_rot = R_to_rpy(RT_rot.block(0,0,3,3));
        cout<<"rpy_rot is :\n"<<rpy_rot<<endl;
        
        pcl::transformPointCloud (cloud_msg_downsampled, cloud_input_aftrot, RT_rot);
*/
        
        pcl::PointCloud<pcl::PointXYZ> cloud_input_aftrot;
        //cloud_input_aftrot = cloud_msg_downsampled_remove;
        cloud_input_aftrot = cloud_msg_downsampled;
        sensor_msgs::PointCloud2 input_aftrot;
        pcl::toROSMsg( cloud_input_aftrot, input_aftrot );
        _pub_input_aftrot.publish(input_aftrot);

        //pcl::PointCloud<pcl::PointXYZ> cloud_msg_debug;
        //cloud_msg_debug.header = cloud_msg.header;;

       /*
        pcl::PointXYZ tmp_point;
        pcl::PointCloud<pcl::PointXYZ>::iterator tmp_itr;
        Vector3d xyz;

        for(tmp_itr = cloud_msg_downsampled.begin(); tmp_itr < cloud_msg_downsampled.end(); tmp_itr++)
        {
            xyz << tmp_itr->x,
                   tmp_itr->y,
                   tmp_itr->z;

            xyz = _R_imu.inverse() * xyz;
            tmp_point.x = xyz(0);
            tmp_point.y = xyz(1);
            tmp_point.z = xyz(2);
            cloud_msg_debug.push_back(tmp_point);
        } 
        //cout<<"size is :"<<cloud_msg_debug.size()<<endl;

        sensor_msgs::PointCloud2 input_debug;
        pcl::toROSMsg( cloud_msg_debug, input_debug );
        _pub_debug.publish(input_debug);
*/
        if(is_laser_init)
        {   
            pc_msg_time = pc_msg.header.stamp.toSec();
            tkf = pc_msg.header.stamp;
            is_laser_init = false;
            cloud_target = cloud_input_aftrot;   
            pointcloud_env = pc_msg;
            pointcloud_env.header.frame_id = "/map";
            RT_kframe = MatrixXd::Identity(4, 4);// RT_init; //RT_icp * RT_kframe;
            R_imu_kf  = R_imu;//MatrixXd::Identity(3, 3);
            RT_kframe(2,3) = height;
            T_kf = RT_kframe.block(0,3,3,1);
            
            Matrix3d R_ = RT_kframe.block(0,0,3,3);
            Q_kf = R_;
            pointcloud_ref = pc_msg;
            kfid ++;
            tkf = pc_msg.header.stamp;
            return;
        }

        _pub_enviro.publish(pointcloud_env);
        cloud_input = cloud_input_aftrot;

        pair<double, Matrix4d> icp_result = IcpFunction(cloud_input, cloud_target);
        pcl::PointCloud<pcl::PointXYZ> cloud_output;
        RT_icp = icp_result.second;
        ros::Time time3 = ros::Time::now();

        double fitness  = icp_result.first;

        cout<<"Now fitness is : "<<fitness<<endl;
        
        pcl::transformPointCloud (cloud_input, cloud_output, RT_icp);

        sensor_msgs::PointCloud2 output_;
        sensor_msgs::PointCloud2 target_;

        pcl::toROSMsg(cloud_output, output_ );
        pcl::toROSMsg(cloud_target, target_ );
        
        _pub_ref.publish( pointcloud_ref);        
        _pub_target.publish(target_);
        _pub.publish( output_ );

        if(fitness > 2.0)
        {   
             is_abort = true;
             ROS_WARN("ICP Abort ###############################");
        }

        Vector3d T   = RT_odom.block( 0, 3, 3, 1); 
        Matrix3d R   = RT_odom.block( 0, 0, 3, 3);

        Vector3d ypr = R_to_ypr(RT_icp.block(0,0,3,3));
        bool rotation_change = false;

        if(abs(ypr(0)) > M_PI/6.0 )// || ypr(1) > 0.1 || ypr(2) > 0.1)
        {
            ROS_WARN("################################################# FLAG ...");
            rotation_change = true;
            //is_abort = false;
        }
        
        if(!is_abort)
          RT_odom = RT_kframe * RT_icp ;
        
        /*
        {
            ROS_WARN("Laser Odom: Stop ...");
            is_abort = true;
            RT_kframe = Matrix4d::Identity(4,4);
            //RT_kframe(2,3) = height;
            counter_good_laser = 0;
            cloud_target = cloud_input;
            kfid ++ ;
        }
        else if(is_abort && (fitness < 3.0) )
        {
            counter_good_laser ++ ;
            if(counter_good_laser > 5)
            {   
                ROS_WARN("Laser Odom: Start ... ");
                is_abort = false;
            }
        }*/

        double x_key_thre = _key_thre;
        double y_key_thre = _key_thre;
        double z_key_thre = _key_thre/8.0;

        /*
        if(cloud_msg_downsampled.size() < 1000 )
           {
                z_key_thre = z_key_thre/2;
                x_key_thre = x_key_thre/2;
                y_key_thre = y_key_thre/2;
           }
        */
        quadrotor_msgs::Odometry laser_quad_odom;
        laser_quad_odom.curodom.header.frame_id         =  string("/map");
        
        if(is_abort)
            laser_quad_odom.curodom.header.frame_id     =  string("null");    

        if(is_abort || abs(T(0) - RT_kframe(0,3)) > x_key_thre || abs(T(1) - RT_kframe(1,3)) > y_key_thre ||abs(T(2) - RT_kframe(2,3)) > z_key_thre || rotation_change)// ||fitness > fitness_thre)//counter % 5 ==0)
        //if( T_dif.norm() > _key_thre || rotation_change || abs(T_dif(2)) > 0.3 || (fitness > fitness_thre))
        {   
            pc_msg_time = pc_msg.header.stamp.toSec();
            cloud_target = cloud_input;

            RT_kframe = RT_odom;
            //if(abs(RT_kframe(2,3) - height) < 0.2)
                RT_kframe(2,3) = height;
            //RT_kframe = Matrix4d::Identity(4,4);
            RT_icp = Matrix4d::Identity(4,4);
            RT_init = Matrix4d::Identity(4,4);
            
            ROS_WARN("###################### CHANGE K FRAME #####################");
            is_change_k = true;
            kfid ++;
            tkf = pc_msg.header.stamp;
            T_kf = RT_kframe.block(0,3,3,1);
            Matrix3d R_ = RT_kframe.block(0,0,3,3);
            Q_kf = R_;
            R_imu_kf = R_imu;//_R_imu;
        }

        Quaterniond Q;
        Q = R;
        Vector3d ypr_vis = R_to_ypr(R);
        ROS_WARN("Icp Yaw is : ");
         cout<<ypr_vis(0)<<endl;

        laser_odom.header.frame_id = "/map";
        laser_odom.pose.pose.orientation.x =Q.x();
        laser_odom.pose.pose.orientation.y =Q.y();
        laser_odom.pose.pose.orientation.z =Q.z();
        laser_odom.pose.pose.orientation.w =Q.w();
        
        laser_odom.pose.pose.position.x = T(0);
        laser_odom.pose.pose.position.y = T(1);
        laser_odom.pose.pose.position.z = T(2);

        laser_odom.header.stamp = pc_msg.header.stamp;

        laser_path.header.frame_id = string("/map");
        laser_path.header.stamp = pc_msg.header.stamp;
        
        cloud_map.header.frame_id = "/map";
/*    
    counter_map++;
    if(counter_map > 20)
    {
    pcl::PointCloud<pcl::PointXYZ> cloud_msg;
    
    pcl::fromROSMsg(pc_msg, cloud_msg);
    cloud_msg.header.frame_id = "/map";
    pcl::PointXYZ tmp_point;
    pcl::PointCloud<pcl::PointXYZ>::iterator tmp_itr;

        pcl::PointCloud<pcl::PointXYZ> cloud_build_map;
        cloud_build_map.header = cloud_msg.header;
        
       // Eigen::Vector3d ypr_pre;
       // double yaw = -0.0/180.0 * M_PI;
        //ypr_pre << yaw,
        //           0,
       //            0;

       // Eigen::Matrix3d R_pre = ypr_to_R(ypr_pre);

        //Matrix4d RT_odom_ = RT_odom;
        //RT_odom_.block(0,0,3,3) = RT_odom.block(0,0,3,3) * R_pre;
        
        pcl::transformPointCloud (cloud_msg, cloud_build_map, RT_odom);
        
        for(tmp_itr = cloud_build_map.begin(); tmp_itr < cloud_build_map.end(); tmp_itr++)
        {
            tmp_point.x = tmp_itr->x;
            tmp_point.y = tmp_itr->y;
            tmp_point.z = tmp_itr->z;
            if(tmp_point.z == 0 || tmp_point.z > 5.0)
              continue;
 
            cloud_map.push_back(tmp_point);
        } 

        pcl::PointCloud<pcl::PointXYZ> cloud_msg_filter;
        cloud_msg_filter = cloud_map;        
        
        pcl::PointCloud<pcl::PointXYZ> cloud_msg_downsampled;
        pcl::VoxelGrid<pcl::PointXYZ>  voxelSampler;

        voxelSampler.setInputCloud(cloud_msg_filter.makeShared());
        voxelSampler.setLeafSize(0.5f, 0.5f, 0.5f);
        voxelSampler.filter(cloud_msg_downsampled);

        pcl::PointCloud<pcl::PointXYZ> cloud_msg_downsampled_remove;
        pcl::StatisticalOutlierRemoval<pcl::PointXYZ> statFilter;
        statFilter.setInputCloud(cloud_msg_downsampled.makeShared());
        statFilter.setMeanK(10);
        statFilter.setStddevMulThresh(1.0);
        statFilter.filter(cloud_msg_downsampled_remove);
        
        sensor_msgs::PointCloud2 cloud_map_;
        pcl::toROSMsg( cloud_msg_downsampled_remove, cloud_map_ );
        _pub_map.publish(cloud_map_);
        counter_map = 0;
    }
*/
        geometry_msgs::PoseStamped pose;
        pose.pose.position.x = T(0);
        pose.pose.position.y = T(1);
        pose.pose.position.z = T(2);
        pose.pose.orientation.w = Q.w();
        pose.pose.orientation.x = Q.x();
        pose.pose.orientation.y = Q.y();
        pose.pose.orientation.z = Q.z();
        pose.header = laser_path.header;
        laser_path.poses.push_back(pose);
       
        laser_quad_odom.kfid = kfid;
        laser_quad_odom.status = quadrotor_msgs::Odometry::STATUS_ODOM_VALID;
        // Current odom 
        laser_quad_odom.curodom.header.stamp            =  pc_msg.header.stamp;
        //laser_quad_odom.curodom.header.frame_id         =  string("/map");
        laser_quad_odom.curodom.pose.pose.position      =  laser_odom.pose.pose.position;
        laser_quad_odom.curodom.pose.pose.orientation   =  laser_odom.pose.pose.orientation;

    Matrix3d cov_base;
    cov_base<<    0.0001,     0.00,       0.00,
                  0.00,       0.0001,     0.00,
                  0.00,       0.00,       0.00001;
    
    //fitness = min(fitness, 2.5);
    Matrix3d cov = sqrt(fitness/10.0) * cov_base;
    cov(2,2) = cov_base(2,2);

    //Matrix3d cov = cov_base;
    for (int j = 0; j < 2; j++)
      for (int i = 0; i < 2; i++)
        laser_quad_odom.curodom.pose.covariance[i+j*6] = cov(i,j);
    
    laser_quad_odom.curodom.pose.covariance[0+3*6] = cov(0,2);
    laser_quad_odom.curodom.pose.covariance[1+3*6] = cov(1,2);
    laser_quad_odom.curodom.pose.covariance[3+0*6] = cov(2,0);
    laser_quad_odom.curodom.pose.covariance[3+1*6] = cov(2,1);
    laser_quad_odom.curodom.pose.covariance[3+3*6] = cov(2,2);

    // Reference odom
    laser_quad_odom.kfodom.header.stamp            = tkf;
    laser_quad_odom.kfodom.header.frame_id         = string("/map");
    laser_quad_odom.kfodom.pose.pose.position.x    = T_kf[0];
    laser_quad_odom.kfodom.pose.pose.position.y    = T_kf[1];
    laser_quad_odom.kfodom.pose.pose.position.z    = T_kf[2];
    laser_quad_odom.kfodom.pose.pose.orientation.x   = Q_kf.x();
    laser_quad_odom.kfodom.pose.pose.orientation.y   = Q_kf.y();
    laser_quad_odom.kfodom.pose.pose.orientation.z   = Q_kf.z();
    laser_quad_odom.kfodom.pose.pose.orientation.w   = Q_kf.w();
    // Publish
    
    _pub_odom.publish(laser_odom);
    _pub_laser_path.publish(laser_path);
    _pub_laser_quad_odom.publish(laser_quad_odom);

     
        ros::Time time4 = ros::Time::now();

        
        
       // cout<<"All Time is:"<<(time4 - time0).toSec()<<endl;
        /*
        cout<<"Filter1 Time is:"<<(time1 - time0).toSec()<<endl;
        cout<<"Filter2 Time is:"<<(time2 - time1).toSec()<<endl;
        cout<<"Align Time is:"<<(time3 - time2).toSec()<<endl;
        cout<<"Other Time is:"<<(time4 - time3).toSec()<<endl;
        */
        is_abort = false;
}

int main(int argc, char** argv) 
{

        // Initialize ROS
        ros::init (argc, argv, "icp_node");
        ros::NodeHandle n("~");

        // Create ROS subscriber for target_PointCloud
        ros::Subscriber sub_lidar = n.subscribe("/velodyne_points", 1, pointcloud_callback );
        ros::Subscriber sub_imu   = n.subscribe("/djiros/imu", 100, imu_callback );
        ros::Subscriber sub_range   = n.subscribe("/teraranger/timeofflight", 10, range_callback );
        //ros::Subscriber odom_range   = n.subscribe("/quadrotor_ukf/odom", 10, odom_callback );

        // Create ROS publisher for transformed pointcloud
        _pub = n.advertise<sensor_msgs::PointCloud2>("transformed_point_cloud", 1 );
        _pub_ref = n.advertise<sensor_msgs::PointCloud2>("reference_point_cloud", 1 );
        _pub_target = n.advertise<sensor_msgs::PointCloud2>("target_point_cloud", 1 );
        _pub_input_aftrot = n.advertise<sensor_msgs::PointCloud2>("input_aftrot", 1 );
        _pub_rawdata = n.advertise<sensor_msgs::PointCloud2>("cloud_rawdata", 1 );
        _pub_enviro = n.advertise<sensor_msgs::PointCloud2>("cloud_env", 1 );
        _pub_map = n.advertise<sensor_msgs::PointCloud2>("cloud_map", 1 );
        //_pub_match = n.advertise<sensor_msgs::PointCloud2>("~matched_point_cloud", 1 );
        _pub_odom = n.advertise<nav_msgs::Odometry>("laser_odom", 10 );
        _pub_laser_path = n.advertise<nav_msgs::Path>("laser_path", 10 );
        _pub_laser_quad_odom = n.advertise<quadrotor_msgs::Odometry>("laser_quad_odom",10);
        //nh.setParam("/global_param", 5);
        n.param("itr_cycle", _itr_cycle, 5);
        n.param("fitness_thre", fitness_thre, 0.25);
        n.param("corres_dis", _corres_dis,8.0);
        n.param("key_thre", _key_thre, 3.0);
        // Spin
        ros::spin ();
}
