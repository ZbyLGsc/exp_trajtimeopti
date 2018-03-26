#ifndef _DATA_TYPE_
#define _DATA_TYPE_

#define inf 9999999.0
#define _PI M_PI

using namespace std;

struct Node;
typedef Node * NodePtr;

struct Point
{
      double x, y, z;
};
  
struct Node
{     
      double radius; // radius of this node
      
      int id;
      Point  coord;

      vector<NodePtr> father_ptr;
      vector<double> dis; // distance to father nodes

      NodePtr preNode_ptr;

      double g, h, f;
};

typedef vector<NodePtr> NodeVec;

#endif