#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <fstream>
#include <geometry_msgs/PoseStamped.h>
#include <pcl_conversions/pcl_conversions.h>

#include <Eigen/Eigen>

#include <iostream>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <random>

#include "quadrotor_msgs/PositionCommand.h"
#include "quadrotor_msgs/PolynomialTrajectory.h"
#include "visualization_msgs/MarkerArray.h"
#include "visualization_msgs/Marker.h"
#include "pcd_trajectory/gridNode.h"
#include "pcd_trajectory/qp_generator.h"
#include "pcd_trajectory/backward.hpp"

#include <arc_utilities/voxel_grid.hpp>
#include "sdf_tools/SDF.h"
#include "sdf_tools/sdf.hpp"
#include "arc_utilities/voxel_grid.hpp"
#include "arc_utilities/pretty_print.hpp"
#include "sdf_tools/collision_map.hpp"
#include "arc_utilities/dynamic_spatial_hashed_voxel_grid.hpp"
#include "sdf_tools/dynamic_spatial_hashed_collision_map.hpp"

namespace backward {
backward::SignalHandling sh;
}

using namespace Eigen;

double _vis_traj_width = 0.1;
ros::Publisher  _traj_vis_pub;
ros::Publisher  _path_vis_pub;
ros::Publisher  _map_vis_pub;

ros::Subscriber _odom_sub;
ros::Subscriber _dest_pts_sub;
ros::Subscriber _map_sub;

double map_x, map_y, map_z;          // length(size) of the map in x, y, z directions, unit is m
int    num_x, num_y, num_z, num_obs; // voxel num of the the map in x, y, z directions
double resolution = 0.2;
double cell_time;

Vector3d startPt, endPt;
State startState;
double _start_x, _start_y, _start_z; 
double _end_x  , _end_y  , _end_z; 

double _vel, _acc, _max_vel, _max_acc;
double _init_vel_x, _init_vel_y, _init_vel_z;
double _init_acc_x, _init_acc_y, _init_acc_z;
bool targetReceive = false;
bool has_map  = false;
bool has_path = false;

int  _traj_order = 5;
int  _his_length = 3; // how many histotic nodes a node would remember
int  _segment_num = _his_length - 1;

void visPath( vector<GridNodePtr> grid_path );
void visTraj( vector<CellTraj> trajectory);

Vector3d mapOrigin(-10.0, -10.0, 0.0);
double _x_size = 20.0, _y_size = 20.0, _z_size = 2.0;
    
Translation3d origin_translation( 0.0 + mapOrigin(0), 0.0 + mapOrigin(1), 0.0);
Quaterniond origin_rotation(1.0, 0.0, 0.0, 0.0);

Affine3d origin_transform = origin_translation * origin_rotation;
sdf_tools::COLLISION_CELL oob_cell(0.0);
sdf_tools::CollisionMapGrid collision_map(origin_transform, "map", resolution, _x_size, _y_size, _z_size, oob_cell);

int64_t max_x = collision_map.GetNumXCells();
int64_t max_y = collision_map.GetNumYCells();
int64_t max_z = collision_map.GetNumZCells();

TrajectoryGenerator _trajGenerator;

pair< vector<GridNodePtr>, vector<CellTraj> > KinoPathSearch(Vector3d start_pt, Vector3d end_pt);
pair< vector<GridNodePtr>, vector<CellTraj> > testQpSolver();

vector<double> getDesiredState( const CellTraj & traj, double t_now );
vector<double> getDesiredFullState( const CellTraj & traj, double t_now );

GridNodePtr GridNodeMap[200][200][40]; 
void initGridMap();

double heuristicCostEstimate(GridNodePtr newNodePtr, GridNodePtr endNodePtr)
{   
    // estimate the planning cost by using one node and the target node
    double h = 0.0;

    return h;
}

MatrixXd PathStitch(GridNodePtr curNodePtr, vector<GridNodePtr> hisNodeList)
{
    // reformulate the historic nodes into a matrix format
    // curNode  is the target of the tiny trajectory
    // hisNodes are previous waypoints for the tiny trajectory

    MatrixXd path(_his_length, 3);
    reverse(hisNodeList.begin(), hisNodeList.end());

    for(int i = 0; i < (int)hisNodeList.size(); i ++ )
    {
        auto ptr = hisNodeList[i];
        //cout<<"ptr->state.pos: \n"<<ptr->state.pos<<endl;
        path.row(i) = ptr->state.pos;
    }

    path.row(_his_length - 1) = curNodePtr->state.pos;

    return path;
}

inline Idx pos2idx(Vector3d location)
{
    vector<int64_t> pt_idx = collision_map.LocationToGridIndex(location(0), location(1), location(2));
    Idx index(pt_idx[0], pt_idx[1], pt_idx[2]);

    return index;
}

inline GridNodePtr pos2grid(Vector3d pos)
{
    vector<int64_t> idx = collision_map.LocationToGridIndex( min(max(-_x_size/2.0, pos(0)), _x_size/2.0), min(max(-_y_size/2.0, pos(1)), _y_size/2.0), min(max(0.0, pos(2)), _z_size/2.0) );
    vector<double> round_pos = collision_map.GridIndexToLocation(idx[0], idx[1], idx[2]);

    Idx index(idx[0], idx[1], idx[2]);
    pos << round_pos[0], round_pos[1], round_pos[2];

    GridNodePtr grid_ptr = new GridNode(index, pos);

    return grid_ptr;
}

Vector3d idx2pos(Idx index)
{   
    vector<int64_t> pt_idx = {(int64_t)index.idx, (int64_t)index.idy, (int64_t)index.idz};
    vector<double> pos = collision_map.GridIndexToLocation(pt_idx[0], pt_idx[1], pt_idx[2]);

    Vector3d location;
    location << pos[0], pos[1], pos[2];

    return location;
}

inline Vector3d grid2pos(GridNodePtr gridPtr)
{
    vector<double> round_pos = collision_map.GridIndexToLocation(gridPtr->index.idx, gridPtr->index.idy, gridPtr->index.idz);

    Vector3d pos;
    pos << round_pos[0], round_pos[1], round_pos[2];
    return pos;
}

double dynamicCostEstimate(GridNodePtr curNodePtr, vector<GridNodePtr> hisNodeList)
{
    // get the planning cost by using the hostoric nodes list
    double qp_cost;
    State init_state = hisNodeList.back()->state;

    MatrixXd path = PathStitch(curNodePtr, hisNodeList);
    MatrixXd poly = _trajGenerator.PolyQPGeneration( path, init_state, cell_time, qp_cost);

    return qp_cost;
}

void Task()
{       
    // auto result = testQpSolver();
    if( has_map == false || has_path == true )
        return;

    auto result = KinoPathSearch(startPt, endPt);
    vector<GridNodePtr> path = result.first;
    vector<CellTraj>    traj = result.second;

    ROS_WARN("Path size is %d", (int)path.size());
    visPath( path );
    visTraj( traj );

    has_path = true;
}

visualization_msgs::MarkerArray cube_vis; 
void visPath( vector<GridNodePtr> grid_path )
{
    for(auto & mk: cube_vis.markers) 
        mk.action = visualization_msgs::Marker::DELETE;

    _path_vis_pub.publish(cube_vis);
    cube_vis.markers.clear();

    visualization_msgs::Marker mk;
    mk.header.frame_id = "map";
    mk.header.stamp = ros::Time::now();
    mk.ns = "traj/path";
    mk.type = visualization_msgs::Marker::CUBE;
    mk.action = visualization_msgs::Marker::ADD;

    mk.pose.orientation.x = 0.0;
    mk.pose.orientation.y = 0.0;
    mk.pose.orientation.z = 0.0;
    mk.pose.orientation.w = 1.0;
    mk.color.a = 0.2;
    mk.color.r = 1.0;
    mk.color.g = 0.0;
    mk.color.b = 0.0;

    int idx = 0;
    for(int i = 0; i < int(grid_path.size()); i++)
    {
        mk.id = idx;

        Vector3d pos = idx2pos(grid_path[i]->index);
        mk.pose.position.x = pos(0); 
        mk.pose.position.y = pos(1); 
        mk.pose.position.z = pos(2);  

        mk.scale.x = resolution;
        mk.scale.y = resolution;
        mk.scale.z = resolution;

        idx ++;
        cube_vis.markers.push_back(mk);
    }

    _path_vis_pub.publish(cube_vis);
}

void visTraj(vector<CellTraj> trajectory)
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
    _traj_vis.points.clear();

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    vector<double> state;
    geometry_msgs::Point pt;

    for(int i = 0; i < (int)trajectory.size(); i++ ){
        for (double t = 0.0; t < cell_time; t += 0.01, count += 1){
            state = getDesiredState(trajectory[i], t);
            cur(0) = pt.x = state[0];
            cur(1) = pt.y = state[1];
            cur(2) = pt.z = state[2];
            _traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);
    _traj_vis_pub.publish(_traj_vis);
}

pcl::PointCloud<pcl::PointXYZ> cloud;
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 & pointcloud_map)
{    
    pcl::fromROSMsg(pointcloud_map, cloud);
    
    if((int)cloud.points.size() == 0)
        return;

    vector<pcl::PointXYZ> inflatePts;
    for (int idx = 0; idx < (int)cloud.points.size(); idx++)
    {   
        auto mk = cloud.points[idx];
        pcl::PointXYZ pt(mk.x, mk.y, mk.z);

        Vector3d addPt(pt.x, pt.y, pt.z);
        sdf_tools::COLLISION_CELL obstacle_cell(1.0); // Occupancy values > 0.5 are obstacles
        collision_map.Set3d(addPt, obstacle_cell);            
    }

    std_msgs::ColorRGBA filled_color, free_color, unknown_color; 
    filled_color.a  = 1.0;   filled_color.b  = 0.0; filled_color.g  = 1.0; filled_color.r  = 1.0;    
    free_color.a    = 0.001; free_color.b    = 1.0; free_color.g    = 1.0; free_color.r    = 1.0;
    unknown_color.a = 0.5;   unknown_color.b = 1.0; unknown_color.g = 0.0; unknown_color.r = 0.0;

    visualization_msgs::Marker collision_map_marker = collision_map.ExportForDisplay(filled_color, free_color, unknown_color);
    collision_map_marker.ns = "collision_map";
    collision_map_marker.id = 1;
    _map_vis_pub.publish(collision_map_marker);
    
    initGridMap();
    has_map = true;
}

void rcvWaypointsCallback(const nav_msgs::Path & wp)
{     
      if(wp.poses[0].pose.position.z < 0.0)
        return;

      endPt <<  wp.poses[0].pose.position.x;
                wp.poses[0].pose.position.y;
                wp.poses[0].pose.position.z;

      targetReceive = true;
}

void rcvOdometryCallbck(const nav_msgs::Odometry odom)
{

}

int main (int argc, char** argv) 
{        
    ros::init (argc, argv, "Kinodynamic_Astar");
    ros::NodeHandle n( "~" );
    
    n.param("KinoPlanning/start_x", _start_x,  0.0);
    n.param("KinoPlanning/start_y", _start_y,  0.0);
    n.param("KinoPlanning/start_z", _start_z,  0.0);

    n.param("KinoPlanning/end_x", _end_x,  5.0);
    n.param("KinoPlanning/end_y", _end_y,  5.0);
    n.param("KinoPlanning/end_z", _end_z,  0.0);

    n.param("KinoPlanning/init_vel_x", _init_vel_x,  0.0);
    n.param("KinoPlanning/init_vel_y", _init_vel_y,  0.0);
    n.param("KinoPlanning/init_vel_z", _init_vel_z,  0.0);
    n.param("KinoPlanning/init_acc_x", _init_acc_x,  0.0);
    n.param("KinoPlanning/init_acc_y", _init_acc_y,  0.0);
    n.param("KinoPlanning/init_acc_z", _init_acc_z,  0.0);

    n.param("KinoPlanning/vel", _vel,    1.0);
    n.param("KinoPlanning/acc", _acc,    0.5);
    n.param("KinoPlanning/max_vel", _max_vel,  2.0);
    n.param("KinoPlanning/max_acc", _max_acc,  1.0);

    startPt << _start_x, _start_y, _start_z;
    endPt   << _end_x,   _end_y,   _end_z;

    startState.pos = Vector3d(_start_x,    _start_y,    _start_z   );
    startState.vel = Vector3d(_init_vel_x, _init_vel_y, _init_vel_z);
    startState.acc = Vector3d(_init_acc_x, _init_acc_y, _init_acc_z);

    cell_time = resolution / _vel;

    // subcribed msgs
    _dest_pts_sub = n.subscribe( "waypoints",  1, rcvWaypointsCallback  );
    _map_sub      = n.subscribe( "map",        1, rcvPointCloudCallBack );
    _odom_sub     = n.subscribe( "odometry",  50, rcvOdometryCallbck    );

    // publish visualize related msgs
    _traj_vis_pub  = n.advertise<visualization_msgs::Marker>("trajectory_vis", 1);
    _path_vis_pub  = n.advertise<visualization_msgs::MarkerArray>("grid_path_vis", 1);
    _map_vis_pub   = n.advertise<visualization_msgs::Marker>("display_distance_field", 1);

    ros::Rate rate(10);
    bool status = ros::ok();
    while (status) {
      ros::spinOnce();
      Task();
      status = ros::ok();
      rate.sleep();
    }
}

vector<double> getDesiredState( const CellTraj & traj, double t_now )
{
    vector<double > ret(3, 0);
    int poly_num1D = _traj_order + 1;

    for ( int dim = 0; dim < 3; dim++ ){
        VectorXd coeff = traj.coeff.row(dim);
        VectorXd t = VectorXd::Zero( poly_num1D );
        
        for(int j = 0; j < poly_num1D; j ++)
            t(j) = pow(t_now, j);
            
            ret[dim] = coeff.dot(t);
    }
  return ret;
}

vector<double> getDesiredFullState( const CellTraj & traj, double t_now )
{
    vector<double > ret(3 * 3, 0);
    int poly_num1D = _traj_order + 1;  

    for ( int dim = 0; dim < 3; dim++ ){
        VectorXd coeff = traj.coeff.row(dim);
        MatrixXd t     = MatrixXd::Zero( 3, poly_num1D );
        
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

pair< vector<GridNodePtr>, vector<CellTraj> > testQpSolver()
{
    // This function is only used for testing the QP solution
    // We construct a small program by hand
    // 
    GridNodePtr cur_node_ptr = pos2grid(startPt);

    Vector3d his_node_pos1, his_node_pos2;
    
    his_node_pos1 = cur_node_ptr->coord;
    his_node_pos1(0) -= resolution;

    his_node_pos2 = his_node_pos1;
    his_node_pos2(1) -= resolution;

    vector<GridNodePtr> his_node_ptrs;

    GridNodePtr his_node_ptr1 = pos2grid(his_node_pos1);
    his_node_ptrs.push_back(his_node_ptr1);

    GridNodePtr his_node_ptr2 = pos2grid(his_node_pos2);
    his_node_ptrs.push_back(his_node_ptr2);
    
    // MatrixXd PathStitch(GridNodePtr curNodePtr, vector<GridNodePtr> hisNodeList)
    MatrixXd test_path = PathStitch(cur_node_ptr, his_node_ptrs);

    //ROS_WARN("Enter in QP solver");
    double qp_cost;
    //cout<<"test_path \n"<<test_path<<endl;
    //cout<<"cell_time: "<<cell_time<<endl;
    MatrixXd polyCoeff = _trajGenerator.PolyQPGeneration(test_path, startState, cell_time, qp_cost);
    //ROS_WARN("Finish QP solver");
    //cout<<"polyCoeff: \n"<<polyCoeff<<endl;

    vector<CellTraj> cell_traj_list;
    int poly_num = _traj_order + 1;
    for(int i = 0; i < _segment_num; i ++)
    {
        CellTraj cell_traj(_traj_order);

        for(int m = 0; m < 3; m ++)
            cell_traj.coeff.row(m) = polyCoeff.block(i, m * poly_num, 1, poly_num);

        cell_traj_list.push_back(cell_traj);
    }

    ROS_WARN("Test Finish");
    ROS_WARN("The QP cost is :%f", qp_cost);
    his_node_ptrs.push_back(cur_node_ptr);
    return make_pair(his_node_ptrs, cell_traj_list);
}

vector<GridNodePtr> retrievePath(GridNodePtr current)
{   
    vector<GridNodePtr> path;
    path.push_back(current);

    while(current->cameFrom != NULL)
    {
        current = current -> cameFrom;
        path.push_back(current);
    }

    return path;
}

vector<CellTraj> retrieveTraj(GridNodePtr current)
{   
    vector<CellTraj> traj;
    //

    return traj;
}

void initGridMap()
{   
    if(has_map == true)
        return;

    GridNodePtr tmpNodePtr = NULL;
    for(int64_t i = 0; i < max_x; i++){
        for(int64_t j = 0; j < max_y; j++){
            for(int64_t k = 0; k < max_z; k++){
                Idx tmpIdx(i,j,k);
                Vector3d pos = idx2pos(tmpIdx);
                tmpNodePtr = new GridNode(tmpIdx, pos);
                tmpNodePtr->id = 0; //initilize id with 0 (neither openset nor closedset)
                tmpNodePtr->cameFrom = NULL;
                tmpNodePtr->gScore = inf;
                tmpNodePtr->fScore = inf;
                tmpNodePtr->occupancy = collision_map.Get(i, j, k).first.occupancy;
                GridNodeMap[i][j][k] = tmpNodePtr;
            }
        }
    }
}

#if 1
pair< vector<GridNodePtr>, vector<CellTraj> > KinoPathSearch(Vector3d start_pt, Vector3d end_pt)
{   
    ros::Time time_1 = ros::Time::now();    
    vector<GridNodePtr> path;
    vector<CellTraj>    traj;
    GridNodePtr startPtr = pos2grid(start_pt);
    GridNodePtr endPtr   = pos2grid(end_pt);
    multimap<double, GridNodePtr> openSet;

    GridNodePtr neighborPtr = NULL;
    GridNodePtr current = NULL;

    startPtr -> gScore = 0;
    startPtr -> fScore = heuristicCostEstimate(startPtr, endPtr);
    startPtr -> id = 1; //put start node in open set
    startPtr -> coord = start_pt;
    openSet.insert( make_pair(startPtr -> fScore, startPtr) ); //put start in open set

    double tentative_gScore;
    /*cout<<"coord: \n"<<startPtr->coord<<endl;
    cout<<"Index: \n"<<startPtr->index.idx<<" , "<<startPtr->index.idy<<" , "<<startPtr->index.idz<<endl;
    cout<<"coord: \n"<<endPtr->coord<<endl;
    cout<<"Index: \n"<<endPtr->index.idx<<" , "<<endPtr->index.idy<<" , "<<endPtr->index.idz<<endl;*/

    int num_iter = 0;
    while ( !openSet.empty() )
    {   
        num_iter ++;
        current = openSet.begin() -> second;

        cout<<"iter: "<<num_iter<<endl;
        if(current->index.idx == endPtr->index.idx
        && current->index.idy == endPtr->index.idy
        && current->index.idz == endPtr->index.idz )
        {
            //cout << "goal coord: " << endl << current->real_coord << endl; 
            ROS_WARN("[Astar]Reach goal..");
            cout << "total number of  iteration used in Astar: " << num_iter  << endl;
            ros::Time time_2 = ros::Time::now();
            ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
            path = retrievePath(current);
            traj = retrieveTraj(current);
            return make_pair(path, traj);
        }         

        openSet.erase(openSet.begin());
        current -> id = -1; //move current node from open set to closed set.

        for(int dx = -1; dx < 2; dx++)
            for(int dy = -1; dy < 2; dy++)
                for(int dz = -1; dz < 2; dz++){
                    if(dx == 0 && dy == 0 && dz ==0){
                        continue; //skip the current point itself
                    }

                    Idx neighborIdx;
                    neighborIdx.idx = (current -> index).idx + dx;
                    neighborIdx.idy = (current -> index).idy + dy;
                    neighborIdx.idz = (current -> index).idz + dz;

                    if(    neighborIdx.idx < 0 || neighborIdx.idx >= max_x 
                        || neighborIdx.idy < 0 || neighborIdx.idy >= max_y
                        || neighborIdx.idz < 0 || neighborIdx.idz >= max_z){
                        continue;
                    }

                    neighborPtr = GridNodeMap[neighborIdx.idx][neighborIdx.idy][neighborIdx.idz];
                    /*cout<<"Index: \n"<<neighborIdx.idx<<" , "<<neighborIdx.idy<<" , "<<neighborIdx.idz<<endl;
                    cout<<"Index: \n"<<neighborPtr->index.idx<<" , "<<neighborPtr->index.idy<<" , "<<neighborPtr->index.idz<<endl;
                    cout<<"coord: \n"<<neighborPtr->coord<<endl;*/

                    //if(collision_map.Get(neighborIdx.idx, neighborIdx.idy, neighborIdx.idz).first.occupancy > 0.5){
                    if(neighborPtr->occupancy > 0.5){
                        continue;
                    }

                    if(neighborPtr -> id == -1){
                        continue; //ignore neighbor which is already evaluated (in closed set).
                    }

                    /*double static_cost = sqrt(dx * dx + dy * dy + dz * dz);                    
                    tentative_gScore = current -> gScore + static_cost; */

                    GridNodePtr ptr = current->cameFrom;
                    vector<GridNodePtr> his_nodes;

                    int i = 0;
                    while(ptr != NULL && i < _his_length - 1)
                    {   
                        his_nodes.push_back(ptr);
                        ptr = ptr->cameFrom;
                        i ++;
                    }
                    //cout<<"his_nodes size: "<<his_nodes.size()<<endl;
                    MatrixXd test_path = PathStitch(current, his_nodes);
                    
                    State init_state = his_nodes[0]->state;
                    double qp_cost;
                    MatrixXd polyCoeff = _trajGenerator.PolyQPGeneration(test_path, init_state, cell_time, qp_cost);

                    double dynamic_cost = qp_cost;                    
                    double static_cost = sqrt(dx * dx + dy * dy + dz * dz);

                    tentative_gScore = current -> gScore + dynamic_cost + static_cost; 

                    if(neighborPtr -> id != 1){
                        //discover a new node
                        neighborPtr -> id        = 1;
                        neighborPtr -> cameFrom  = current;
                        neighborPtr -> gScore    = tentative_gScore;
                        neighborPtr -> fScore    = neighborPtr -> gScore + heuristicCostEstimate(neighborPtr, endPtr); 
                        neighborPtr -> state     =  

                        neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                        continue;
                    }
                    else if(tentative_gScore <= neighborPtr-> gScore){ //in open set and need update
                        neighborPtr -> cameFrom = current;
                        neighborPtr -> gScore   = tentative_gScore;
                        neighborPtr -> fScore   = tentative_gScore + heuristicCostEstimate(neighborPtr, endPtr); 
                        neighborPtr -> state    = 
                        openSet.erase(neighborPtr -> nodeMapIt);
                        neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                    }
                        
                }
    }

    ros::Time time_2 = ros::Time::now();
    ROS_WARN("Time consume in A* path finding is %f", (time_2 - time_1).toSec() );
    return make_pair(path, traj);
}
#endif