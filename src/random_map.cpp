#include <iostream>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/kdtree.h>
#include <pcl/search/impl/kdtree.hpp>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>

#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <math.h>

#include <time.h>
#include <sys/time.h>
#include <random>

#include "visualization_msgs/MarkerArray.h"
#include "visualization_msgs/Marker.h"
#include "quadrotor_msgs/PositionCommand.h"

#define inf 99999999.0

using namespace std;

pcl::search::KdTree<pcl::PointXYZ> kdtreeLocalMap;
vector<int>     pointIdxRadiusSearch;
vector<float>   pointRadiusSquaredDistance;        

random_device rd;
default_random_engine eng(rd());
uniform_real_distribution<double>  rand_x;
uniform_real_distribution<double>  rand_y;
uniform_real_distribution<double>  rand_w;
uniform_real_distribution<double>  rand_h;

ros::Publisher _all_map_pub, _local_map_pub, _local_map_nofloor_pub;
ros::Publisher _raw_map_pub;
ros::Publisher _delta_map_pub;

ros::Subscriber _map_sub;
ros::Subscriber _odom_sub;
ros::Subscriber _cmd_sub;

deque<nav_msgs::Odometry> _odom_queue;
vector<double> _state;
const size_t _odom_queue_size = 200;
nav_msgs::Odometry _odom;

pcl::PointXYZ LstMapCenter;

// Map model : 
int _obsNum;
double _x_l, _x_h, _y_l, _y_h, _w_l, _w_h, _h_l, _h_h;
double _resolution;
double _sense_rate;
double _z_limit;
double _sensing_range;
double _safe_margin = 0.1;

// Sensor model : 
double _field_of_view_horizontal, _field_of_view_vertical, _angular_resolution;

//ros::Timer vis_map;
bool map_ok = false;
bool _has_odom = false;
bool _has_ground;

sensor_msgs::PointCloud2 globalMap_pcd, localMap_pcd, localMap_nofloor_pcd, deltaMap_pcd;
sensor_msgs::PointCloud2 rawMap_pcd;
sensor_msgs::PointCloud2 localMapInflate_pcd;
pcl::PointCloud<pcl::PointXYZ> cloudMap, cloudMap_nofloor;

void RandomMapGenerate()
{    
#define ADD_GROUND_PLANE _has_ground

      pcl::PointXYZ pt_random;
      cloudMap.clear();

      // #####################################################
      // Generate a random map, with vertical obstacles
      rand_x = uniform_real_distribution<double>(_x_l, _x_h );
      rand_y = uniform_real_distribution<double>(_y_l, _y_h );
      rand_w = uniform_real_distribution<double>(_w_l, _w_h);
      rand_h = uniform_real_distribution<double>(_h_l, _h_h);

      for(int i = 0; i < _obsNum; i ++){
         double x, y; 
         x    = rand_x(eng);
         y    = rand_y(eng);

         double w, h;
         w    = rand_w(eng);

         int widNum = ceil(w/_resolution);

         for(int r = -widNum / 2; r < widNum / 2; r ++ )
            for(int s = -widNum / 2; s < widNum / 2; s ++ ){
               h    = rand_h(eng);  
               //if(h < 1.0) continue;
               int heiNum = ceil(h/_resolution);
               for(int t = 0; t < heiNum; t ++ ){
                  pt_random.x = x + r * _resolution;
                  pt_random.y = y + s * _resolution;
                  pt_random.z = t * _resolution;
                  cloudMap.points.push_back( pt_random );
               }
            }
      }
      cloudMap_nofloor = cloudMap;
      cloudMap_nofloor.width = cloudMap_nofloor.points.size();
      cloudMap_nofloor.height = 1;
      cloudMap_nofloor.is_dense = true;

if(ADD_GROUND_PLANE)
{     
      ROS_WARN("[Map Generator] Has ground plane");
      uniform_real_distribution<double>  rand_g(0.0, 1.0);
      int x_all = (_x_h - _x_l) / _resolution;
      int y_all = (_y_h - _y_l) / _resolution;

      for(int i = -x_all; i < x_all; i ++)
      {
         for(int j = -y_all; j < y_all; j ++)
         {
            double p_g = rand_g(eng);
            
            if(p_g < 0.45) continue;
            pt_random.x = i * _resolution;
            pt_random.y = j * _resolution;
            pt_random.z = 0.0;
            cloudMap.points.push_back( pt_random );
         }
      }
}

      cloudMap.width = cloudMap.points.size();
      cloudMap.height = 1;
      cloudMap.is_dense = true;

      ROS_WARN("[Map Generator] Finished generate random map ");
      
      cout<<cloudMap.size()<<endl;
      
      kdtreeLocalMap.setInputCloud( cloudMap.makeShared() ); 

      map_ok = true;
}

void rcvOdometryCallbck(const nav_msgs::Odometry odom)
{
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

        _odom_queue.push_back(odom);
        while (_odom_queue.size() > _odom_queue_size) _odom_queue.pop_front();       
}

void rcvPosCmdCallBack(const quadrotor_msgs::PositionCommand cmd)
{
      _has_odom = true;
      _odom.pose.pose.position.x = cmd.position.x;
      _odom.pose.pose.position.y = cmd.position.y;
      _odom.pose.pose.position.z = cmd.position.z;
      
      _odom.twist.twist.linear.x = cmd.velocity.x;
      _odom.twist.twist.linear.y = cmd.velocity.y;
      _odom.twist.twist.linear.z = cmd.velocity.z;
      _odom.header = cmd.header;

      _state = {
          _odom.pose.pose.position.x, 
          _odom.pose.pose.position.y, 
          _odom.pose.pose.position.z, 
          _odom.twist.twist.linear.x,
          _odom.twist.twist.linear.y,
          _odom.twist.twist.linear.z,
          0.0, 0.0, 0.0
      };

      _odom_queue.push_back(_odom);
      while (_odom_queue.size() > _odom_queue_size) _odom_queue.pop_front();
}

#if 0
void SensorModelFilter(pcl::PointCloud<pcl::PointXYZ>::Ptr ScanPtr)
{     
      //return;
      ROS_INFO("The points number before the filter in the scan is %d", int(ScanPtr->points.size()));
      int ray_num = int(_field_of_view_horizontal / _angular_resolution);
      pair<int, pair<int, double> > angleState[ray_num]; // The first element is the state whether this angle has a point, the second element is the idx of the nearest point in pointcloud at this angle, and the nearest poiint's distance

      for(int i = 0; i < ray_num; i ++ )
         angleState[i].first = 0; // 0 indicates no point at this direction

      pcl::PointXYZ _pt;
      double angle = 0.0, dis = 0.0;
      double dx, dy;

      for(int i = 0; i < int(ScanPtr->points.size()); i ++ ){
         _pt = ScanPtr->points[i];
         dx = _pt.x - _state[0];
         dy = _pt.y - _state[1];
         //dz = _pt.z - _state[2];
         
         if(dx > 0.0 && dy >= 0.0)       // Case 1: the point lays in the 1st orphant
            angle = atan(dy / dx)   * 180.0 / M_PI;
         else if(dx <= 0.0 && dy > 0.0)  // Case 2: the point lays in the 2nd orphant
            angle = atan(dy / -dx)  * 180.0 / M_PI + 90.0;
         else if(dx < 0.0  && dy <= 0.0) // Case 3: the point lays in the 3rd orphant
            angle = atan(-dy / -dx) * 180.0 / M_PI + 180.0;
         else if(dx >= 0.0 && dy < 0.0)  // Case 4: the point lays in the 4th orphant
            angle = atan(-dy / dx)  * 180.0 / M_PI + 270.0;

         dis = sqrt(pow(dx, 2) + pow(dy, 2));
         double ag;       
/*         if(dis < 2.5) ag = _angular_resolution * 4;
         else if(dis < 5.0) ag = _angular_resolution * 2;*/
         //else 
            ag = _angular_resolution;
         
         int ray_idx = int(angle / ag);

         //dis = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );

         //ROS_INFO("The angle of the ray is %f, the index of the ray is : %d, the distance is: %f", angle, ray_idx, dis);

         if( angleState[ray_idx].first == 0 ){
            angleState[ray_idx].first  = 1;
            angleState[ray_idx].second = make_pair(i, dis);
         }else if(dis < angleState[ray_idx].second.second)
            angleState[ray_idx].second = make_pair(i, dis);     
      }

      pcl::PointCloud<pcl::PointXYZ>::Ptr ScanAftFilter(new pcl::PointCloud<pcl::PointXYZ>());

      for(int i = 0; i < ray_num; i ++){
         if(angleState[i].first == 1)
            ScanAftFilter->points.push_back(ScanPtr->points[angleState[i].second.first]);
      }

      ScanPtr->points.clear();
      ScanPtr->points = ScanAftFilter->points;

      ROS_INFO("The points number now in the scan is %d", int(ScanPtr->points.size()) );

}
#endif

void PostFilter(pcl::PointCloud<pcl::PointXYZ>::Ptr ScanPtr)
{
   pcl::PointCloud<pcl::PointXYZ>::Ptr ScanAftFilter(new pcl::PointCloud<pcl::PointXYZ>()); 
   ScanAftFilter->points.clear();
   pcl::PointXYZ robot(_state[0], _state[1], _state[2]);
   //cout<<_state[0]<<","<<_state[1]<<","<<_state[2]<<endl;

   // Basic idea : maintain a occupancy map locally, do ray tracing on the map in a stupid way .
   double length  = _sensing_range;
   double reso    = _resolution;

   int num_grid   = round(2 * length / reso);
   double grid_h[num_grid][num_grid]; // records the higher bound of the height at each grid .

   //ROS_WARN("Grid number is: %d", num_grid);
   for(int i = 0; i < num_grid; i ++ )
      for(int j = 0; j < num_grid; j ++){
         grid_h[i][j] = -1.0;
      }

/*   for(int i = 0; i < num_grid; i ++ ){
      for(int j = 0; j < num_grid; j ++){
         cout<<grid_h[i][j];
      }
      cout<<"\n"<<endl;
   }*/

   /*ROS_ERROR("Debug");
   cout<<int(-19.2)<<endl;*/
   pcl::PointXYZ pt;
   for(int i = 0; i < int(ScanPtr->points.size()); i ++ ){
      pt = ScanPtr->points[i];
      int coord_x = round( (pt.x - round(robot.x)) / reso) + (num_grid/2 - 1), coord_y = round( (pt.y - round(robot.y) )/ reso) + (num_grid/2 - 1), coord_z = round(pt.z/ reso);
      //ROS_INFO("coordinate of point is %d, %d", coord_x, coord_y);
      double height = coord_z;
      if(height > grid_h[coord_x][coord_y]) // formulate the entire grid map
         grid_h[coord_x][coord_y] = height;
   } 
   
/*   ROS_ERROR("height in the cornor");
   cout<<grid_h[99][99]<<endl;
   cout<<grid_h[98][99]<<endl;
   cout<<grid_h[99][98]<<endl;*/

   int ray_num = _field_of_view_horizontal / _angular_resolution;
   int ray_num_orphant = ray_num / 4;
   for(int ray_idx = 0; ray_idx < ray_num_orphant; ray_idx ++ ){ // just in the 1st orphant 
      double pt_x, pt_y, pt_z;
      double k = tan( double(ray_idx) * _angular_resolution / 180.0 * M_PI);

      for(int x_id = num_grid/2 - 1; x_id < num_grid; x_id ++ ){
         int y_id = round((x_id - (num_grid/2 - 1) ) * k + (num_grid/2 - 1));
         
         y_id = min(num_grid, y_id);
         /*cout<<"slope is: "<<k<<endl;
         cout<<"x_id: "<<x_id<<" , "<<"y_id"<<y_id<<endl;*/
         
         if(grid_h[x_id][y_id] >= 1){
            for(int z_id = 0; z_id < grid_h[x_id][y_id]; z_id ++){
                  pt_x = double(x_id - (num_grid/2 - 1) ) * reso + robot.x;
                  pt_y = double(y_id - (num_grid/2 - 1) ) * reso + robot.y;
                  pt_z = double(z_id) * reso;
                  ScanAftFilter->points.push_back(pcl::PointXYZ(pt_x, pt_y, pt_z));
               }
            break;
         }
      }
   }

   ScanPtr->points.clear();
   ScanPtr->points = ScanAftFilter->points;
  /* // Now prepare for doing ray casting
   for(int x_id = 0; x_id < num_grid; x_id ++ )
      for(int y_id = 0; y_id < num_grid; y_id ++ ){
         bool discard = false;
         double pt_x, pt_y, pt_z;        
         if( (x_id > (num_grid/2 - 1) ) && (y_id > (num_grid/2 - 1) ) ){
            if(grid_h[x_id][y_id] <= 0.0) continue;
            pt_x = (x_id - (num_grid/2 - 1) ) * reso;
            pt_y = (y_id - (num_grid/2 - 1) ) * reso;
            double k; // k is the slope rate of the ray casting from the point to the robot .
            k = double(pt_y - robot.y) / double(pt_x - robot.x);
            //k = double((y_id - (num_grid/2 - 1) ) ) / double((x_id - (num_grid/2 - 1) ));
            cout<<"slope is: "<<k<<endl;
            cout<<"x_id: "<<x_id<<" , "<<"y_id"<<y_id<<endl;
            ROS_ERROR("Begin checking");
            for(int x_r = (x_id - 1); x_r >= (num_grid/2 - 1); x_r--){
               cout<<x_r<<endl;
               int y_r =  int( double(x_r - (num_grid/2 - 1) ) * k + (num_grid/2 - 1));
               if( grid_h[x_r][y_r] >= 0.0 ){ // ... something indicate that some point block the ray
                  discard = true; // discard this point
                  break;
               }
            }
            ROS_ERROR("End checking, discard status is %d", discard);

            //ROS_ERROR("discard state: %d", discard );
            if(!discard){
               for(int z_id = 0; z_id < grid_h[x_id][y_id]; z_id ++){
                  pt_x = (x_id - (num_grid/2 - 1) ) * reso + robot.x;
                  pt_y = (y_id - (num_grid/2 - 1) ) * reso + robot.y;
                  pt_z = z_id * reso;
                  ScanAftFilter->points.push_back(pcl::PointXYZ(pt_x, pt_y, pt_z));
               }
            }
         }
      }*/

    // Now prepare for doing ray casting


  // ROS_BREAK();

}

int i = 0;
void pubSensedPoints()
{     
      if(i < 10 ){
         pcl::toROSMsg(cloudMap_nofloor, globalMap_pcd);
         globalMap_pcd.header.frame_id = "map";
         _all_map_pub.publish(globalMap_pcd);
      }

      i ++;
      {
            if(!map_ok || !_has_odom)
               return;

            ros::Time time_bef_sensing = ros::Time::now();
            
            pcl::PointCloud<pcl::PointXYZ>::Ptr deltaMap(new pcl::PointCloud<pcl::PointXYZ>());           
            pcl::PointCloud<pcl::PointXYZ>::Ptr localMap(new pcl::PointCloud<pcl::PointXYZ>());
            pcl::PointCloud<pcl::PointXYZ>::Ptr localMap_nofloor(new pcl::PointCloud<pcl::PointXYZ>());
            pcl::PointXYZ searchPoint(_state[0], _state[1], _state[2]);

            pointIdxRadiusSearch.clear();
            pointRadiusSquaredDistance.clear();
            
            pcl::PointXYZ ptInNoflation;

            if ( kdtreeLocalMap.radiusSearch (searchPoint, _sensing_range, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0 ){
               for (size_t i = 0; i < pointIdxRadiusSearch.size (); ++i){
                  ptInNoflation = cloudMap.points[pointIdxRadiusSearch[i]];      
                  if( abs(ptInNoflation.z - searchPoint.z ) >  _sensing_range * sin(_field_of_view_vertical/ 2.0 / 180.0 * M_PI) )
                     continue;
                  /*if(sqrt(pow(ptInNoflation.x - searchPoint.x, 2) + pow(ptInNoflation.y - searchPoint.y, 2) + pow(ptInNoflation.z - searchPoint.z, 2) ) < 0.5 )
                     continue;*/

                  localMap->points.push_back(ptInNoflation);
                  
                  if(ptInNoflation.z > 0.05) // filter all points on the ground, just for visualization
                     localMap_nofloor->points.push_back(ptInNoflation);

                  if( sqrt(pow(ptInNoflation.x - LstMapCenter.x, 2) + pow(ptInNoflation.y - LstMapCenter.y, 2) + pow(ptInNoflation.z - LstMapCenter.z, 2) ) < _sensing_range )
                     continue;
                  deltaMap->points.push_back(ptInNoflation);
               }

            }
            else{
               ROS_ERROR("[Map Generator] No obstacles .");
               cout<<searchPoint.x<<" , "<<searchPoint.y<<" , "<<searchPoint.z<<endl;
               return;
            }

            localMap->width = localMap->points.size();
            localMap->height = 1;
            localMap->is_dense = true;
               
            pcl::toROSMsg(*localMap, localMap_pcd);
            localMap_pcd.header.frame_id = "map";
            _local_map_pub.publish(localMap_pcd);

            localMap_nofloor->width = localMap_nofloor->points.size();
            localMap_nofloor->height = 1;
            localMap_nofloor->is_dense = true;
               
            pcl::toROSMsg(*localMap_nofloor, localMap_nofloor_pcd);
            localMap_nofloor_pcd.header.frame_id = "map";
            _local_map_nofloor_pub.publish(localMap_nofloor_pcd);

            deltaMap->width = deltaMap->points.size();
            deltaMap->height = 1;
            deltaMap->is_dense = true;
               
            pcl::toROSMsg(*deltaMap, deltaMap_pcd);
            deltaMap_pcd.header.frame_id = "map";
            _delta_map_pub.publish(deltaMap_pcd);

            ros::Time time_aft_sensing = ros::Time::now();
            LstMapCenter = pcl::PointXYZ(_state[0], _state[1], _state[2]);

            /*ROS_ERROR("Time consume in pubsensing is : ");
            cout<<time_aft_sensing - time_bef_sensing<<endl;*/
         }
}

// ---------------boyu----------------------------------
void nmapCallback(const geometry_msgs::PoseStamped::Ptr ps)
{
    if(fabs(ps->pose.position.x - 1) < 1e-3)
        _obsNum = 30;
    else if(fabs(ps->pose.position.x - 2) < 1e-3)
        _obsNum = 50;
    else if(fabs(ps->pose.position.x - 3) < 1e-3)
        _obsNum = 60;
    else if(fabs(ps->pose.position.x - 4) < 1e-3)
        _obsNum = 20;
    else if(fabs(ps->pose.position.x - 5) < 1e-3)
        _obsNum = 40;
    else if(fabs(ps->pose.position.x - 6) < 1e-3)
        _obsNum = 70;

    _obsNum *= 4; 
    RandomMapGenerate();
    i = 0;
}
// ---------------boyu----------------------------------


int main (int argc, char** argv) {
        
      ros::init (argc, argv, "random map generator");
      ros::NodeHandle n( "~" );

      _local_map_pub =
            n.advertise<sensor_msgs::PointCloud2>("RandomMap", 1);                            
      
      _local_map_nofloor_pub =
            n.advertise<sensor_msgs::PointCloud2>("RandomMap_noFloor", 1);                            
      
      _all_map_pub =
            n.advertise<sensor_msgs::PointCloud2>("AllMap", 1);  

      _raw_map_pub   =
            n.advertise<sensor_msgs::PointCloud2>("RawMap", 1);                            
      
      _delta_map_pub =
            n.advertise<sensor_msgs::PointCloud2>("DeltaMap", 1);                            

      _odom_sub      = 
            n.subscribe( "odometry", 50, rcvOdometryCallbck );

      _cmd_sub      = 
            n.subscribe( "position_cmd", 50, rcvPosCmdCallBack );

      n.param("mapBoundary/lower_x", _x_l,       0.0);
      n.param("mapBoundary/upper_x", _x_h,     100.0);
      n.param("mapBoundary/lower_y", _y_l,       0.0);
      n.param("mapBoundary/upper_y", _y_h,     100.0);
      
      n.param("ObstacleShape/lower_rad", _w_l,   0.3);
      n.param("ObstacleShape/upper_rad", _w_h,   0.8);
      n.param("ObstacleShape/lower_hei", _h_l,   3.0);
      n.param("ObstacleShape/upper_hei", _h_h,   7.0);
      n.param("ObstacleShape/z_limit", _z_limit, 5.0);
      n.param("ObstacleShape/has_ground", _has_ground, true);
      
      n.param("LocalBoundary/radius",  _sensing_range, 10.0);
      
      n.param("ObstacleNum", _obsNum,  30);
      n.param("Resolution",  _resolution, 0.2);
      n.param("SensingRate", _sense_rate, 10.0);
      
      n.param("SensorModel/fov_horizon",   _field_of_view_horizontal, 360.0);
      n.param("SensorModel/fov_vertical",  _field_of_view_vertical,   30.0 ); // 15 threads vertically
      n.param("SensorModel/angular_res",   _angular_resolution, 0.25);        // 0.25 degree / ray

      ros::Time time_beg = ros::Time::now();
      RandomMapGenerate();
      ros::Time time_aft = ros::Time::now();
      ROS_WARN("[Map Generator] Time consume is");
      cout<<time_aft - time_beg<<endl;

      // -----------------boyu--------------------
      ros::Subscriber nmap_sub = n.subscribe("trajtest/nmap", 1, nmapCallback);
      // -----------------boyu--------------------

      ros::Rate loop_rate(_sense_rate);
      
      while (ros::ok())
      {
        pubSensedPoints();
        ros::spinOnce();
        loop_rate.sleep();
      }
}