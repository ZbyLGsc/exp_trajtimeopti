#ifndef _DATA_TYPE_
#define _DATA_TYPE_

#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <vector>

#define inf 1>>30
#define _PI M_PI

using namespace std;

struct GridNode;

typedef GridNode* GridNodePtr;

struct Idx
{
      int idx, idy, idz;

      Idx(int id_x, int id_y, int id_z)
      {
         idx = id_x;
         idy = id_y;
         idz = id_z;
      }

      Idx(){};
      ~Idx(){};
};

struct State
{
   Eigen::Vector3d  pos;
   Eigen::Vector3d  vel;
   Eigen::Vector3d  acc;

   State( Eigen::Vector3d _pos,  Eigen::Vector3d _vel,  Eigen::Vector3d _acc)
   {
      pos = _pos;
      vel = _vel;
      acc = _acc;
   }
   
   State(){};
   ~State(){};
};

struct CellTraj
{
   int order;        // the order of the polynomial cell trajectory asigned in this cell
   double cell_time; // the time allocated in each cell
   Eigen::MatrixXd coeff;

   CellTraj( int _order )
   {
      order = _order;
      coeff.resize(3, order + 1);
   }

   ~CellTraj(){}; 
};

struct GridNode
{     
   int id;        // 1--> open set, -1 --> closed set
   Idx index;
   Eigen::Vector3d coord;

   double gScore, fScore;
   GridNodePtr cameFrom;
   std::multimap<double, GridNodePtr>::iterator nodeMapIt;
   double occupancy; 

   vector<GridNodePtr> hisNodeList; // use a list to record nodes in its history
   State  state;
   
   /*CellTraj tinyTraj;
   GridNode(Idx _index, int _order)
   {
      id = 0;
      grid_index = _index;
      gScore = inf;
      fScore = inf;
      cameFrom = NULL;

      tinyTraj = CellTraj(_order);
   }*/

   GridNode(Idx _index, Eigen::Vector3d _coord)
   {  
      id = 0;
      index = _index;
      coord = _coord;
      
      state.pos = _coord;
      state.vel = Eigen::VectorXd::Zero(3);
      state.acc = Eigen::VectorXd::Zero(3);;

      gScore = inf;
      fScore = inf;
      cameFrom = NULL;
   }

   GridNode(){};
   
   ~GridNode(){};
};

#endif