#ifndef _DATA_TYPE_
#define _DATA_TYPE_

#define inf 9999999.0
#define _PI M_PI
#include <unordered_map>  

using namespace std;

struct Node;
typedef Node * NodePtr;

struct Point
{
      float x, y, z;

      Point(float x_, float y_, float z_)
      {
	      	x = x_;
	      	y = y_;
	      	z = z_;
      }
      
      Point(){}

      ~Point(){}
};
  
struct Node
{     
      float radius; // radius of this node
      bool valid;
      bool best;
      bool change;

      int id;

      // temporary variables, only for speedup the tree rewire procedure
      int rel_id;
      float rel_dis;

      Point  coord;
      vector<NodePtr> nxtNode_ptr;
	  
      //unordered_map<int, NodePtr> connect_ptr;
      NodePtr preNode_ptr;
      float g; // total cost of the shortest path from this node to the root
      float f; // heuristic value of the node to the target point
      
      Node( Point xyz_, float radius_, int id_, float g_, float f_)
      {		
      		coord  = xyz_;
      		radius = radius_;
      		id     = id_;
      		g      = g_;  
      		f      = f_;
			
			rel_id  = - 2;   // status undefined
			rel_dis = - 1.0; // distance undifined

      		valid  = true;
      		best   = false; 
	      	change = false;

      		preNode_ptr = NULL;
	      	nxtNode_ptr.clear();
      }

      Node(){}

      ~Node(){}
};

typedef vector<NodePtr> NodeVec;

#endif