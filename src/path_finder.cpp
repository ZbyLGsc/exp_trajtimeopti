#include "pcd_trajectory/path_finder.h"
#include <unordered_set>
#include <queue>
#include <algorithm>

#define fasle false

rrtPathFinder::rrtPathFinder( ){ 
      cach_size  = 100;
}

rrtPathFinder::~rrtPathFinder(){ }

void rrtPathFinder::setParam( float safety_margin_, float search_margin_, float max_radius_, float sample_range_ )
{     
      safety_margin = safety_margin_;
      search_margin = search_margin_;
      max_radius    = max_radius_;
      sample_range  = sample_range_;
}

void rrtPathFinder::reset()
{     
      TreeDestruct();
      
      NodeList.clear();
      EndList.clear();
      sampleSet.clear();
      invalidSet.clear();
      reCheckList.clear();
      PathList.clear();

      best_end_ptr    = new Node();
      rootNode        = new Node();

      path_find_state = true;
      refine_status   = false;
      inform_sample   = false;
      commitEnd       = false;
      best_distance   = inf;  

      node_id    = 0;
      nr_nodes   = 0;
}

void rrtPathFinder::setStartPt( Point startPt, Point endPt)
{
      start_pt = startPt;
      end_pt   = endPt; 

/*      inform_centroid.x = (start_pt.x + end_pt.x) / 2.0;
      inform_centroid.y = (start_pt.y + end_pt.y) / 2.0;

      ctheta = (end_pt.x - start_pt.x) / min_distance; 
      stheta = (end_pt.y - start_pt.y) / min_distance; 
      min_distance = sqrt( pow(start_pt.x - end_pt.x, 2) + pow(start_pt.y - end_pt.y, 2) + pow(start_pt.z - end_pt.z, 2) );*/

      rand_x_in = uniform_real_distribution<float>(start_pt.x - sample_range, start_pt.x + sample_range);
      rand_y_in = uniform_real_distribution<float>(start_pt.y - sample_range, start_pt.y + sample_range);
}

void rrtPathFinder::setPt( Point startPt, Point endPt, float xl, float xh, float yl, float yh, float zl, float zh, float biasl, float biash, 
                           float localRange, int max_iter, float sample_portion, float goal_portion )
{
      start_pt = startPt;
      end_pt   = endPt; 

      eng = default_random_engine(rd()) ;
      
      x_l = xl; x_h = xh;
      y_l = yl; y_h = yh;
      z_l = zl; z_h = zh;

      bias_l = biasl; bias_h = biash;
      
      sample_range = localRange;

      rand_x    = uniform_real_distribution<float>(x_l, x_h);
      rand_y    = uniform_real_distribution<float>(y_l, y_h);
      rand_z    = rand_z_in = uniform_real_distribution<float>(z_l + safety_margin, z_h); //lower bound: z_l --> z_l + safety_margin
      rand_bias = uniform_real_distribution<float>(bias_l, bias_h);

      rand_x_in = uniform_real_distribution<float>(start_pt.x - sample_range, start_pt.x + sample_range);
      rand_y_in = uniform_real_distribution<float>(start_pt.y - sample_range, start_pt.y + sample_range);

      min_distance = sqrt( pow(start_pt.x - end_pt.x, 2) + pow(start_pt.y - end_pt.y, 2) + pow(start_pt.z - end_pt.z, 2) );
      //ROS_ERROR("[Path Finder] min_distance is %f meters", min_distance);
      inform_centroid.x = (start_pt.x + end_pt.x) / 2.0;
      inform_centroid.y = (start_pt.y + end_pt.y) / 2.0;

      ctheta = (end_pt.x - start_pt.x) / min_distance; 
      stheta = (end_pt.y - start_pt.y) / min_distance; 

      max_samples  = max_iter;
      inlier_ratio = sample_portion;
      goal_ratio   = goal_portion; 
}

void rrtPathFinder::setInput(pcl::PointCloud<pcl::PointXYZ> CloudIn)
{     
      pcl::VoxelGrid<pcl::PointXYZ>  VoxelSampler;
      pcl::PointCloud<pcl::PointXYZ> Cloud_DS = CloudIn;
      pcl::PointCloud<pcl::PointXYZ> Cloud_DS_NaN;
      
      ros::Time time1 = ros::Time::now();
      /*VoxelSampler.setLeafSize(0.1f, 0.1f, 0.1f);
      VoxelSampler.setInputCloud( CloudIn.makeShared() );      
      VoxelSampler.filter( Cloud_DS );    
*/
      std::vector<int> indices;
      pcl::removeNaNFromPointCloud(Cloud_DS,Cloud_DS_NaN, indices);

      //ROS_WARN("Input Set OK .. ");
      kdtreeForMap.setInputCloud( Cloud_DS_NaN.makeShared() );    
      /*ROS_WARN("Check points number of pointCloud : %lu", Cloud_DS_NaN.points.size());
      ros::Time time2 = ros::Time::now();
      
      cout<<"time input is: "<<time2-time1<<endl;*/
}

void rrtPathFinder::rcvAddMap(pcl::PointCloud<pcl::PointXYZ> CloudAdd){
      kdtreeAddMap.setInputCloud( CloudAdd.makeShared() );
}

void rrtPathFinder::rcvDelMap(pcl::PointCloud<pcl::PointXYZ> CloudDel){
      kdtreeDelMap.setInputCloud( CloudDel.makeShared() );
}

inline float rrtPathFinder::getDis(const NodePtr node1, const NodePtr node2){
      return sqrt(pow(node1->coord.x - node2->coord.x, 2) + pow(node1->coord.y - node2->coord.y, 2) + pow(node1->coord.z - node2->coord.z, 2) );
}

inline float rrtPathFinder::getDis(const NodePtr node1, const Point pt){
      return sqrt(pow(node1->coord.x - pt.x, 2) + pow(node1->coord.y - pt.y, 2) + pow(node1->coord.z - pt.z, 2) );
}

inline float rrtPathFinder::getDis(const Point p1, const Point p2){
      return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2) );
}

inline float rrtPathFinder::getDisL1(const NodePtr node1, const Point pt){
      return fabs( node1->coord.x - pt.x) + fabs(node1->coord.y - pt.y) + fabs(node1->coord.z - pt.z) ;
}

float rrtPathFinder::radius_search( Point search_Pt)//, pcl::search::KdTree<pcl::PointXYZ>::Ptr Tree)
{     
      if(getDis(search_Pt, start_pt) > sample_range + max_radius )
         return max_radius - search_margin;

      pcl::PointXYZ searchPoint;
      searchPoint.x = search_Pt.x;
      searchPoint.y = search_Pt.y;
      searchPoint.z = search_Pt.z;

      pointIdxRadiusSearch.clear();
      pointRadiusSquaredDistance.clear();

      kdtreeForMap.nearestKSearch(searchPoint, 1, pointIdxRadiusSearch, pointRadiusSquaredDistance);
      float radius = sqrt(pointRadiusSquaredDistance[0]) - search_margin; // The height of flight yard is about 2.6 m
      return min(radius, float(max_radius));
}

float rrtPathFinder::radius_search_add( Point search_Pt)
{     
      pcl::PointXYZ searchPoint;
      searchPoint.x = search_Pt.x;
      searchPoint.y = search_Pt.y;
      searchPoint.z = search_Pt.z;

      pointIdxRadiusSearch.clear();
      pointRadiusSquaredDistance.clear();

      kdtreeAddMap.nearestKSearch(searchPoint, 1, pointIdxRadiusSearch, pointRadiusSquaredDistance);
      float radius = sqrt(pointRadiusSquaredDistance[0]) - search_margin; 
      return min(radius, float(max_radius));
}


void rrtPathFinder::clearBranchS(NodePtr node_delete) // Strong branch cut: no matter how, cut all nodes in this branch
{     
      for( auto nodeptr: node_delete->nxtNode_ptr ){
          if( nodeptr->valid)
              invalidSet.push_back(nodeptr);
              
              nodeptr->valid = false;
              clearBranchS(nodeptr);
      }
}

void rrtPathFinder::TreeSparsify(NodePtr newPtr)
{     
      NodePtr ptr = newPtr;
      if( ptr->g + ptr->f > best_distance ){ // delete it and all its branches
          // cout<<"estimate cost: "<<ptr->g + ptr->f<<endl;
          // cout<<"lower bound: "<<best_distance<<endl;
          ptr->valid = false;
          invalidSet.push_back(ptr);
          clearBranchS(ptr);
      }    
}

void rrtPathFinder::removeInvalid()
{     
      vector<NodePtr> UpdateNodeList;
      vector<NodePtr> UpdateEndList;
      kd_clear(kdTree_);
     // ROS_WARN("[remove invalid] Check nodes num : %d", int( NodeList.size() ) );
/*      int count_end_false = 0;
      for(auto ptr:EndList){
        if(ptr->valid == false){
          //cout<<"cost : "<<ptr->g<<endl;
          count_end_false ++;
        }
      }

      ROS_WARN("[remove invalid] Best distance is %f", best_distance);
      ROS_WARN("[remove invalid] Before remove, EndList size is %d", int(EndList.size() ));
      ROS_WARN("[remove invalid] invalid end pointer num is %d", count_end_false);
*/
      for(auto nodeptr:NodeList){ // Go through all nodes, record all valid ones
          if(nodeptr->valid){

            float pos[3] = {nodeptr->coord.x, nodeptr->coord.y, nodeptr->coord.z};
            kd_insertf(kdTree_, pos, nodeptr);

            UpdateNodeList.push_back(nodeptr);

            if( check_end(nodeptr))
                UpdateEndList.push_back(nodeptr);

          }
      }

      NodeList.clear();
      EndList.clear();

      NodeList = UpdateNodeList;
      EndList  = UpdateEndList;

      // Now we should deal with all invalid nodes, broken all its related connections

      for(int i = 0; i < int(invalidSet.size()); i++ ){
          NodePtr nodeptr = invalidSet[i];

          if(nodeptr->preNode_ptr != NULL){ // let this node's father forget it
              nodeptr->change = true;
              
              vector<NodePtr> child = nodeptr->preNode_ptr->nxtNode_ptr;        
              nodeptr->preNode_ptr->nxtNode_ptr.clear();
              for(auto ptr: child){
                  if( ptr->change == true ) 
                      continue;
                  else
                      nodeptr->preNode_ptr->nxtNode_ptr.push_back(ptr);
              }
          }
      }

      vector<NodePtr> deleteList;
      for(int i = 0; i < int(invalidSet.size()); i++ ){
          NodePtr nodeptr = invalidSet[i];

          for(auto childptr: nodeptr->nxtNode_ptr) // let this node's children forget it
          {   
              if(childptr -> valid) //ROS_BREAK();
                  childptr->preNode_ptr = NULL;
          }

          deleteList.push_back(nodeptr);
      }

      invalidSet.clear();
      for(int i = 0; i < int(deleteList.size()); i++)
      {
          NodePtr ptr = deleteList[i];
          delete ptr;
      }

      //ROS_WARN("[remove invalid] Check nodes num : %d", int( NodeList.size() ) );
      //ROS_WARN("[remove invalid] EndList size is %d", int(EndList.size() ));
}

void rrtPathFinder::clearBranchW(NodePtr node_delete) // Weak branch cut: if a child of a node is on the current best path, keep the child
{     
      for( auto nodeptr: node_delete->nxtNode_ptr ){
          if( nodeptr->best ) continue;
          else{

              if(nodeptr->valid)
                invalidSet.push_back(nodeptr);
              
              nodeptr->valid = false;
              clearBranchW(nodeptr);
          }

      }
}

void rrtPathFinder::resetRoot(Point commitTarget)
{     
      NodePtr lstNode = PathList.front();
      
      if(getDis(lstNode, commitTarget) < lstNode->radius){
          ROS_ERROR("[Path Finder] almost reach the final target, return ");
          commitEnd = true;
          return;
      }
      ROS_WARN("[Path Finder] start reset the ROOT node");
      bool deleteRoot = false;
      double costReduction = 0;

      commit_root = commitTarget;
      vector<NodePtr> cutList;
      
      for(auto nodeptr:NodeList)
          nodeptr->best = false;

      for(auto nodeptr: PathList){
          if( (!deleteRoot) && (getDis(nodeptr, commitTarget) < (nodeptr->radius - 0.1) ) ){ // now the drone is contained in this node's sphere
              ROS_WARN("[Path Finder] This node is the last one we should keep");
              cout<<"x: "<<nodeptr->coord.x<<" , y: "<<nodeptr->coord.y<<" , z: "<<nodeptr->coord.z<<endl;
              deleteRoot = true;
              costReduction = nodeptr->g; 
              nodeptr->best   = true;
              //nodeptr->valid  = true;
              nodeptr->preNode_ptr = NULL;
              rootNode = nodeptr;
              continue;
          }
          if( deleteRoot ){
              ROS_WARN("[Path Finder] This node should be deleted");
              cout<<"x: "<<nodeptr->coord.x<<" , y: "<<nodeptr->coord.y<<" , z: "<<nodeptr->coord.z<<endl;
              cout<<"R: "<<nodeptr->radius<<endl;
              nodeptr->best  = false;
              nodeptr->valid = false;
              //nodeptr->preNode_ptr = NULL;
              cutList.push_back(nodeptr);
          }
      }
      
      costRecast(costReduction, commitTarget);
      
      for(auto nodeptr:cutList){
          invalidSet.push_back(nodeptr);
          clearBranchW(nodeptr);
      }

      ROS_WARN("[path finder] finish label the deleted nodes due to ROOT reset");

      removeInvalid();

      ROS_WARN("[path finder] finish reset the ROOT node");
      ROS_WARN("[path finder] now nodes num is %d", int(NodeList.size()) );

      /*for( int i = 0; i < int(NodeList.size()); i++)
      {   
        NodePtr ptr = NodeList[i];
        cout<<"coordinate: "<<ptr->coord.x<<" , "<<ptr->coord.y<<" , "<<ptr->coord.z<<"R: "<<ptr->radius<<endl;
      }*/
}

void rrtPathFinder::costRecast(double costReduction, Point target)
{
      for(auto nodeptr: NodeList){
          //cout<<"coordinate: "<<nodeptr->coord.x<<" , "<<nodeptr->coord.y<<" , "<<nodeptr->coord.z<<endl;
          nodeptr->g -= costReduction;
      }
      min_distance   = getDis(target, end_pt);
      
      inform_centroid.x = (target.x + end_pt.x) / 2.0;
      inform_centroid.y = (target.y + end_pt.y) / 2.0;

      ctheta = (end_pt.x - target.x) / min_distance; 
      stheta = (end_pt.y - target.y) / min_distance; 

      best_distance -= costReduction;
}

void rrtPathFinder::UpdateHeuristicRegion(NodePtr update_end_node)
{   
    // Now ONLY consider to be used in a static environment 
    // This function update the heuristic hype-ellipsoid sampling region once and better path has been found.
    // Update the up-to-date traversal and conjugate diameter of the ellipsoid.
    // If there is no improvement in the path, maintain the heuristic unchange.

    float update_cost = update_end_node->g + getDis(update_end_node, end_pt) + getDis(rootNode, commit_root);
    //ROS_WARN("Update heuristic, update_end_node->g is : %f", update_end_node->g);

    if(update_cost < best_distance){
  /*      ROS_WARN("best_distance is : %f", update_cost);
        ROS_WARN("minimum distance is: %f", min_distance);
  */      
        best_distance = update_cost;

        elli_l = best_distance;
        elli_s = sqrt(best_distance * best_distance - min_distance * min_distance);

        // release the old best path, free the best status
        if(inform_sample){
            for(auto ptr:NodeList)
                ptr->best = false;
        }

        //ROS_WARN("half updated");
        // update the nodes in the new best path to be marked as best
        NodePtr ptr = update_end_node;
        while( ptr != NULL  ) {  
            ptr->best = true;
            ptr = ptr->preNode_ptr;
        }

        best_end_ptr = update_end_node;
  //      ROS_WARN("Update Once");
    }    
}

Point rrtPathFinder::genSample()
{     
      float x, y, z, bias; 
      Point pt;
      
      bias = rand_bias(eng);
      if( bias <= goal_ratio ){
          pt = end_pt;
          //cout<<"end_pt bias: "<<pt.x<<", "<<pt.y<<" , "<<pt.z<<endl;
          return pt;
      }

      // To generate samples in a heuristic hype-ellipsoid region.
      // Basic idea : 1. generate samples according to (rho, phi), which is a unit circle
      //              2. scale the unit circle to a ellipsoid
      //              3. rotate the ellipsoid

      if(!inform_sample){ 
          if( bias > goal_ratio && bias <= (goal_ratio + inlier_ratio) ){ 
          // sample inside the local map's boundary
              pt.x    = rand_x_in(eng);
              pt.y    = rand_y_in(eng);
              pt.z    = rand_z_in(eng);  
          }
          else{           
          // uniformly sample in all region
              pt.x    = rand_x(eng);
              pt.y    = rand_y(eng);
              pt.z    = rand_z(eng);  
          }
      }
      else{ 
          // Sample in the informed region
          float angle  = rand_phi(eng);
          float radius = rand_rho(eng);

          x = sqrt(radius) * cos(angle) * elli_l / 2.0;
          y = sqrt(radius) * sin(angle) * elli_s / 2.0;
          z = rand_z(eng);  
          
          pt.x = x * ctheta - y * stheta + inform_centroid.x;
          pt.y = x * stheta + y * ctheta + inform_centroid.y;
          pt.z = z;
      }

      /*Generate a random point inside a circle of radius 1. This can be done by taking a random angle phi in the interval [0, 2*pi) 
      and a random value rho in the interval [0, 1) and compute

      x = sqrt(rho) * cos(phi)
      y = sqrt(rho) * sin(phi)

      The square ROOT in the formula ensures a uniform distribution inside the circle. Scale x and y to the dimensions of the ellipse
      
      x = x * width/2.0
      y = y * height/2.0 */

      sampleSet.push_back(pt);
      return pt;
}

NodePtr rrtPathFinder::genNewNode( Point pt_sample, NodePtr node_nearst_ptr )
{
      float dis       = getDis(node_nearst_ptr, pt_sample);

      Point center;
      if(dis > node_nearst_ptr->radius)
      {
          float steer_dis = node_nearst_ptr->radius / dis;
          center.x = node_nearst_ptr->coord.x + (pt_sample.x - node_nearst_ptr->coord.x) * steer_dis;
          center.y = node_nearst_ptr->coord.y + (pt_sample.y - node_nearst_ptr->coord.y) * steer_dis;
          center.z = node_nearst_ptr->coord.z + (pt_sample.z - node_nearst_ptr->coord.z) * steer_dis;
      }
      else
      {
          center.x = pt_sample.x;
          center.y = pt_sample.y;
          center.z = pt_sample.z;
      }
      float radius_  = radius_search( center );
      float h_dis_   = getDis(center, end_pt);

      node_id ++;
      NodePtr node_new_ptr = new Node( center, radius_, node_id, inf, h_dis_ ); //Node( Point xyz_, float radius_, int id_, float g_, gloat f_)
      
      //cout<<"check the id"<<node_new_ptr->id<<endl;
      if(node_new_ptr->id == 0)
        ROS_BREAK();

      return node_new_ptr;
}

bool rrtPathFinder::check_boundary_reach( NodePtr ptr )
{
      float distance = getDis(ptr, start_pt);

      if( (distance + ptr->radius > sample_range) ) // && (distance - ptr->radius < sample_range)
        return true;

      return false;      
}

bool rrtPathFinder::check_end( NodePtr ptr )
{    
      float distance = getDis(ptr, end_pt);
      
      if(distance + 0.1 < ptr->radius)
          return true;      
     
      return false;
}

bool rrtPathFinder::checkTrajInvalid(Point traj_pt)
{     
      //if(radius_search(traj_pt) + search_margin < safety_margin )
      
      if(radius_search(traj_pt) < 0.0 )
        return true;

      return false;
}

inline bool rrtPathFinder::CheckConnect( float dis, NodePtr node_1, NodePtr node_2 )
{     
      if( ((dis + 0.15) < 0.90 * (node_1->radius + node_2->radius)) && (dis > node_1->radius) && (dis > node_2->radius) ) 
          return true;
      else 
          return false;
}

#if 1
NodePtr rrtPathFinder::findNearstVertex( Point pt_sample )
{
    //kdres * nearest = kd_nearest3( kdTree_, pt_sample.x, pt_sample.y, pt_sample.z );
    float pos[3] = {pt_sample.x, pt_sample.y, pt_sample.z};
    kdres * nearest = kd_nearestf( kdTree_, pos );
    
    /*if (kd_res_size(nearest) <= 0) {
      kd_res_free(nearest);
      ROS_ERROR("[RRT path finder] Error in nearest vertex searching"); 
      return NULL;
    }*/

    NodePtr node_nearst_ptr = (NodePtr) kd_res_item_data( nearest );
    kd_res_free(nearest);

    return node_nearst_ptr;
}
#endif

inline int  rrtPathFinder::CheckRelation( float dis, NodePtr node_1, NodePtr node_2 )
{
      // -1 indicate two nodes are connected good
      //  0 indicate two nodes are not connected
      //  1 indicate node_2 contains node_1, node_1 should be deleted
      //  2 indicate node_1 contains node_2, node_2 should be deleted. Or should be checked and rejected in samling / node_generation steps ?

      int status;
      
      /*if( (dis + node_1->radius) <= node_2->radius)
          status = 1;*/
      if( (dis + node_2->radius) == node_1->radius)
          status = 2;
      else if( (dis + 0.1) < 0.95 * (node_1->radius + node_2->radius) ) // && (dis > node_1->radius) && (dis > node_2->radius)
          status = -1;
      else
          status = 0;
      
      return status;
}

#if 1
void rrtPathFinder::TreeRewire( NodePtr node_new_ptr, NodePtr node_nearst_ptr)
{     
      NodePtr newPtr     = node_new_ptr;      
      NodePtr nearestPtr = node_nearst_ptr;

      float range = newPtr->radius * 2.0f;
      float pos[3] = {newPtr->coord.x, newPtr->coord.y, newPtr->coord.z};
      struct kdres *presults  = kd_nearest_rangef( kdTree_, pos, range);

      vector<NodePtr> nearPtrList;
      bool isInvalid = false;
      while( !kd_res_end( presults ) ) { // Pull all the nodes outside the result data structure from the kd-tree
                                         // And go through all the nearby vertex, check the relations, to see whether the newPtr should be discarded
            NodePtr nearPtr = (NodePtr)kd_res_item_data( presults );

            float dis = getDis( nearPtr, newPtr );
            int   res = CheckRelation( dis, nearPtr, newPtr );
            nearPtr->rel_id  = res;    // temporary variables to record the local realtions with nearby nodes, to avoid repeated calculation in this two terms
            nearPtr->rel_dis = dis;

            if( res == 2 ){
                //cout<<"res is 2 : "<<newPtr->coord.x<<" , "<<newPtr->coord.y<<" , "<<newPtr->coord.z<<endl;
                newPtr->valid = false;
                isInvalid = true;
                nearPtrList.push_back(nearPtr);
                break;
            }
            else{
                nearPtrList.push_back(nearPtr);
                kd_res_next( presults );
            }
      }
      kd_res_free( presults );

      if(isInvalid){
          for(auto nodeptr: nearPtrList){ // release all the temporary variables
              nodeptr->rel_id  = -2;
              nodeptr->rel_dis = -1.0;
          }
          return;
      }
      float min_cost = nearestPtr->g +  getDis( nearestPtr, newPtr );
      newPtr->preNode_ptr = nearestPtr;
      newPtr->g = min_cost;
      nearestPtr->nxtNode_ptr.push_back(newPtr);
      NodePtr lstParentPtr = nearestPtr; // a pointer records the last best father of the new node

      vector<NodePtr> nearVertex;

      // Go through all the nearby vertex again, to deal with another two case. Note, better not mix this go through with the previous one
      for(auto nodeptr:nearPtrList){
        NodePtr nearPtr = nodeptr;
        
        // Choose the parent
        int res   = nearPtr->rel_id;
        float dis = nearPtr->rel_dis;
        float cost = nearPtr->g + dis;
        
        if( res == -1 ){
            if( cost < min_cost ){ // has a shorter path if new_ptr's father is the new node
                min_cost = cost;
                newPtr->preNode_ptr = nearPtr;
                newPtr->g = min_cost;
                lstParentPtr->nxtNode_ptr.pop_back();
                lstParentPtr = nearPtr; // change a father, record the last father
                lstParentPtr->nxtNode_ptr.push_back(newPtr);
            }
            nearVertex.push_back(nearPtr);
        }      
        nearPtr->rel_id  = -2;
        nearPtr->rel_dis = -1.0;
      }

      //ROS_WARN("Finish choose the parent node");
      // *** Rewire the neighbor *** //
      for(int i = 0; i < int(nearVertex.size()); i++ ){
          NodePtr nodeptr = nearVertex[i];
          
          NodePtr nearPtr = nodeptr;
          if(nearPtr->valid == false) continue;

          float dis = getDis( nearPtr, newPtr );
          float cost = dis + newPtr->g;
          
          if( cost < nearPtr->g){ // a near Node changed a father, delete this node from its parent's childer list
              
              if(isSuccessor(nearPtr, newPtr->preNode_ptr)) 
                  continue; // By experiments, this function does work here
              
              if(nearPtr->preNode_ptr == NULL ){
                  nearPtr->preNode_ptr = newPtr;
                  nearPtr->g = cost;
              }
              else{
                  NodePtr lstNearParent = nearPtr->preNode_ptr;
                  nearPtr->preNode_ptr = newPtr;
                  nearPtr->g = cost;
                  
                  nearPtr->change = true; // use a temporary flag to indicate this pointer should be deleted
                  vector<NodePtr> child = lstNearParent->nxtNode_ptr;
                      
                  lstNearParent->nxtNode_ptr.clear();
                  for(auto ptr: child){
                      if( ptr->change == true ) continue;
                      else lstNearParent->nxtNode_ptr.push_back(ptr);
                  }
                  nearPtr->change = false; // release the flag after the delete  
              }
            
              newPtr->nxtNode_ptr.push_back(nearPtr);
          }
      }
}
#endif

void rrtPathFinder::AddtoGraph(NodePtr new_node){
      //if(new_node->preNode_ptr == NULL) {ROS_WARN("FUck"); ROS_BREAK();}
      NodeList.push_back(new_node);
      nr_nodes = int(NodeList.size());
}

void rrtPathFinder::RRTpathFind( float path_time_limit )
{     
      double duration = 0.0;
      ros::Time time_bef_pathFind = ros::Time::now();
      
      kdTree_ = kd_create(3); //2
      commit_root  = start_pt;
      NodePtr root_ = new Node(start_pt, radius_search( start_pt ), node_id, 0.0, min_distance);  // Node( Point xyz_, float radius_, int id_, float g_, gloat f_)

      AddtoGraph(root_);

      float pos[3] = {root_->coord.x, root_->coord.y, root_->coord.z};
      kd_insertf(kdTree_, pos, root_);
      rootNode = root_;
      int iter_count;

      for( iter_count = 0; iter_count < max_samples; iter_count ++)
      {     
          ros::Time time_in_loop = ros::Time::now();
          duration = (time_in_loop - time_bef_pathFind).toSec();

          if( duration > (double)path_time_limit ) 
              break;
          
          Point   pt_sample       =  genSample();
              
          NodePtr node_nearst_ptr = findNearstVertex(pt_sample);
          
          if(!node_nearst_ptr->valid || node_nearst_ptr == NULL ) 
            continue;
                    
          NodePtr node_new_ptr  =  genNewNode(pt_sample, node_nearst_ptr); 
          
          if( node_new_ptr->coord.z < z_l  || node_new_ptr->radius < safety_margin ) 
            continue;
             
          TreeRewire(node_new_ptr, node_nearst_ptr);  
          
          if( ! node_new_ptr->valid) continue;

          if(check_end( node_new_ptr ))
          {

              if( !inform_sample ) 
                  best_end_ptr = node_new_ptr;
              
              EndList.push_back(node_new_ptr);
              UpdateHeuristicRegion(node_new_ptr);
              inform_sample = true;  
          }
          
          pos[0] = node_new_ptr->coord.x;
          pos[1] = node_new_ptr->coord.y;
          pos[2] = node_new_ptr->coord.z;
          kd_insertf(kdTree_, pos, node_new_ptr);

          AddtoGraph(node_new_ptr);
          TreeSparsify(node_new_ptr);      

          if( int(invalidSet.size()) >= cach_size) removeInvalid();
      }

      removeInvalid();

      ros::Time time_aft_pathFind = ros::Time::now();
      ROS_WARN("[path finder] Time in path finding function, check time consumption");
      cout<<time_aft_pathFind - time_bef_pathFind<<endl;

      ROS_WARN("[path finder] Check nodes num : %d", int( NodeList.size() ) );
      ROS_WARN("[path finder] Check feasible path num : %d", int( EndList.size() ) );

/*      int count_null = 0;
      for(auto ptr: NodeList)
        if(ptr->preNode_ptr == NULL)
          count_null ++;

      ROS_WARN("[path refine] pre null nodes num : %d", count_null);
      if(count_null > 1) ROS_BREAK();*/
      /*ROS_WARN("[path find] check the best_distance");
      cout<<best_distance<<endl;*/

      tracePath();
      ROS_WARN("[path refine] Check feasible path num : %d", int( EndList.size() ) );
      /*for(auto ptr:EndList)
      {
        cout<<"coordinate: "<<ptr->coord.x<<" , "<<ptr->coord.y<<" , "<<ptr->coord.z<<"R: "<<ptr->radius<<endl;
        cout<<"id: "<<ptr->id<<endl;
      }*/
}

void rrtPathFinder::tracePath()
{   
      vector<NodePtr> feasibleEndList;      
          
      for(auto endPtr: EndList)
      {
          if( (CheckValidEnd(endPtr) == false) || (check_end(endPtr) == false) || (endPtr->valid == false) )
              continue;
          else
              feasibleEndList.push_back(endPtr);
      }

      if( int(feasibleEndList.size()) == 0 )
      {
          ROS_WARN("[trace path] can't find a feasible path, ### ");
          path_find_state = false;
          best_distance = inf;
          inform_sample = false;
          EndList.clear();
          Path   = Eigen::MatrixXd::Identity(3,3);
          Radius = Eigen::VectorXd::Zero(3);
          return;
      }
  
      EndList = feasibleEndList;

      best_end_ptr = feasibleEndList[0];
      float best_cost = inf;
      for(auto nodeptr: feasibleEndList)
      {   
          float cost = (nodeptr->g + getDis(nodeptr, end_pt) + getDis(rootNode, commit_root) );
          if( cost < best_cost )
          {
              best_end_ptr  = nodeptr;
              best_cost = cost;
              best_distance = best_cost;
          }
      }

      //best_distance = inf;
      NodePtr ptr = best_end_ptr;

/*      cout<<"coordinate: "<<ptr->coord.x<<" , "<<ptr->coord.y<<" , "<<ptr->coord.z<<"R: "<<ptr->radius<<endl;
      cout<<"check end of this node: "<<check_end(ptr)<<endl;*/

      int idx = 0;
      PathList.clear();

      while( ptr != NULL ) { 
          PathList.push_back( ptr );
          ptr = ptr->preNode_ptr;
          idx ++;
      }

      Path.resize(idx + 2, 3) ;
      Radius.resize( idx );
      idx = 0;
      
      ptr = best_end_ptr;
      Path.block< 1, 3 >( idx, 0 ) << start_pt.x, start_pt.y, start_pt.z;      
      idx ++ ;

      while( ptr != NULL ){      
            Path.block< 1, 3 >( Path.rows() - idx - 1 , 0 ) << ptr->coord.x, ptr->coord.y, ptr->coord.z;
            Radius[ Radius.size() - idx] = ptr->radius;
            
            ptr = ptr->preNode_ptr;
            idx ++;
      }

      Path.block< 1, 3 >( idx, 0 ) << end_pt.x, end_pt.y, end_pt.z;
      
      path_find_state = true;
      ROS_WARN("[trace path] now have %d feasible path", int(feasibleEndList.size()) );
 /*     cout<<"Path:\n"<<Path<<endl;
      cout<<"cost: "<<best_end_ptr->g<<endl;*/
}

pair<Eigen::MatrixXd, Eigen::VectorXd> rrtPathFinder::getPath()
{
    return make_pair(Path, Radius);
}

vector<NodePtr> rrtPathFinder::getTree()
{
    /*int conut_child = 0;
    for(auto nodeptr:NodeList){
        //nodeptr->id = 0;
        conut_child += nodeptr->nxtNode_ptr.size();
    }*/    
    
    /*ROS_ERROR("NodeList Size is : %d", int(NodeList.size()));
    ROS_ERROR("Children Number is : %d", conut_child);*/
    return NodeList;
}

bool rrtPathFinder::RefineStatus()
{
    return refine_status;
}

void rrtPathFinder::TreeDestruct()
{
    //ROS_WARN("prepare to destruct the tree");
    kd_free(kdTree_);
    
    for( int i = 0; i < int(NodeList.size()); i++)
    {   
        NodePtr ptr = NodeList[i];
        //Nodefree(ptr);
        delete ptr;
    }

    //ROS_WARN("finish to destruct the tree");
}

float rrtPathFinder::reCheckRadius(Point pt)
{
      return radius_search( pt );
}

float rrtPathFinder::reCheckAddRadius(Point pt)
{
      return radius_search_add( pt );
}

int rrtPathFinder::ReCheckConnection(float new_radius, float old_radius)
{
      if( new_radius < safety_margin)  
          return -1; // -1 for not qualified, this node should be deleted
      else if( new_radius < old_radius )
          return 0;  // 0 for node radius shrink, this node should be inspected for connection with others 
      else 
          return 1;  // 1 for no difference, nothing should be changed
}


bool rrtPathFinder::isSuccessor(NodePtr curPtr, NodePtr nearPtr) // check if curPtr is father of nearPtr
{     
      NodePtr prePtr = nearPtr->preNode_ptr;
      
      while(prePtr != NULL)
      {
          if(prePtr == curPtr){
            //ROS_BREAK(); // seemed this function works
            return true;
          }
          else
            prePtr = prePtr->preNode_ptr;
      }

      return false;
}

void rrtPathFinder::RRTpathReEvaluate( float revaluate_limit, bool isPrint)
{     
      //if(isPrint) ROS_WARN("[path re-evaluate] Start re-evaluate some nodes");
      
      if(path_find_state == false)
      {
          ROS_WARN("[path re-evaluate] no path exists. ");
          return;
      }
      
      ros::Time time_bef_evaluate = ros::Time::now();
      vector< pair<Point, double> > evaFailList;
      while(true)
      {   
          for(int i = 0; i < int(PathList.size()); i++ ) 
          {   
              NodePtr ptr = PathList[i];
/*              if( ptr->valid == false ) //continue;
                ROS_BREAK();  // here may has BUGS, don't know why there will be invalid node in the path*/
              NodePtr pre_ptr =  ptr->preNode_ptr;
              
              if( pre_ptr != NULL )
              {
                  float update_radius = reCheckRadius(ptr->coord);
                  int ret = ReCheckConnection(update_radius, ptr->radius); // -1: not qualified, 0: shrink, 1: no difference, continue
                  float old_radius = ptr->radius;
                  ptr->radius = update_radius; // the radius of a node may shrink or remain no change, but can not enlarge

                  if( ret == -1)
                  {   // The node ptr now does't have enough volume fot the robot to pass, the node should be deleted from the tree; 
                      // all its child nodes should be put in the rewire list 
                      //ROS_WARN("Node invalid now, delete");
                      ptr->valid = false;          // Delete this node
                      //ptr->best  = false;
                      invalidSet.push_back(ptr);
                      clearBranchS(ptr);
                      evaFailList.push_back(make_pair(ptr->coord, old_radius));
                  }
                  else 
                  {   // Test whether the node should be disconnected with its parent and all children,
                      // If disconnect with parent, delete it, if discoonect with a child, delete the child
                      //ROS_WARN("node shrinked, check it and its children");
                      vector<NodePtr> childList = ptr->nxtNode_ptr;
                      for(auto nodeptr: childList)
                      {  // inspect each child of ptr, to see whether they are still connected 
                          float dis = getDis( ptr, nodeptr );
                          int res = CheckRelation( dis, ptr, nodeptr );
                          if( res != -1 ) // the child is disconnected with its parent
                          {   
                              if(nodeptr->valid == true)
                              {
                                  nodeptr->valid = false;      
                                  invalidSet.push_back(nodeptr);
                                  clearBranchS(nodeptr);
                                  evaFailList.push_back(make_pair(nodeptr->coord, nodeptr->radius));
                              }
                          } 
                      }
                  }
              }
          }

          bool isBreak = true;

          for(auto ptr:PathList)
          { 
            //cout<<ptr->valid<<endl;
            isBreak = isBreak && ptr->valid;
          }

          //removeInvalid();

          if(isBreak)
            break;

          //if(isPrint) ROS_WARN("[path re-evaluate] finish inspect all nodes along the previous best path");
          ROS_WARN("[path re-evaluate] now %d nodes reach target", int(EndList.size()) );

          vector<NodePtr> feasibleEndList;                
          for(auto endPtr: EndList)
          {
              if( (endPtr->valid == false)  || (check_end(endPtr) == false) )
                  continue;
              else
                  feasibleEndList.push_back(endPtr);
          }
          EndList = feasibleEndList;

          ROS_WARN("[path re-evaluate] now %d path feasible", int(feasibleEndList.size()) );
          //double repair_limit = revaluate_limit - (time_aft_evaluate - time_bef_evaluate).toSec();
          ros::Time time_in_evaluate = ros::Time::now();
          if(feasibleEndList.size() == 0 || (time_in_evaluate - time_bef_evaluate).toSec() > revaluate_limit )
          {
              path_find_state = false;
              inform_sample = false;
              best_distance = inf;
              break;
          }
          else
          {
              best_end_ptr = feasibleEndList[0];
              float best_cost = inf;
              for(auto nodeptr: feasibleEndList)
              {
                  float cost = (nodeptr->g + getDis(nodeptr, end_pt) + getDis(rootNode, commit_root) );
                  if( cost < best_cost )
                  {
                      best_end_ptr  = nodeptr;
                      best_cost = cost;
                      best_distance = best_cost;
                  }
              }
          
              PathList.clear();
              NodePtr ptrr = best_end_ptr;
              while( ptrr != NULL ){
                PathList.push_back( ptrr );
                ptrr = ptrr->preNode_ptr;
              }
          }
          //ROS_WARN("[path re-evaluate] finish one iteration ");
        } 

        removeInvalid();
        if(path_find_state == true)
        {
            NodePtr ptr = PathList[0];

            int idx = 0;
            PathList.clear();

            while( ptr != NULL ) { 
                PathList.push_back( ptr );
                ptr = ptr->preNode_ptr;
                idx ++;
            }

            Path.resize(idx + 2, 3) ;
            Radius.resize( idx );
            idx = 0;
            
            ptr = best_end_ptr;
            Path.block< 1, 3 >( idx, 0 ) << start_pt.x, start_pt.y, start_pt.z;      
            idx ++ ;

            while( ptr != NULL ){      
                  Path.block< 1, 3 >( Path.rows() - idx - 1 , 0 ) << ptr->coord.x, ptr->coord.y, ptr->coord.z;
                  Radius[ Radius.size() - idx] = ptr->radius;
                  
                  ptr = ptr->preNode_ptr;
                  idx ++;
            }
            Path.block< 1, 3 >( idx, 0 ) << end_pt.x, end_pt.y, end_pt.z;
        }            
        //if(isPrint) ROS_WARN("[path re-evaluate] finish re-evaluate the path");
        
        treeRepair(0.0, evaFailList);

        tracePath();
}


void rrtPathFinder::treeRepair(double repair_limit, vector< pair<Point, double> > repairList)
{
        ros::Time time_bef_treeRepair = ros::Time::now();
        // check nearby nodes in all the failure nodes' neighbour, 
        // try to repair the tree, since where the nodes failed previously should be most possible fail again

        int repairSize = int(repairList.size());

        //ROS_WARN("[tree repair] the size of the repair list is %d", repairSize );
        for(int i = 0; i < repairSize; i++)
        {   
            Point center = repairList[i].first;
            //cout<<"coordinate: "<<center.x<<" , "<<center.y<<" , "<<center.z<<" , R: "<<repairList[i].second<<endl;

            float range  = repairList[i].second * 2.0f;
            float pos[3] = {center.x, center.y, center.z};
            struct kdres *presults  = kd_nearest_rangef( kdTree_, pos, range);

            while( !kd_res_end( presults ) ) 
            { 
                NodePtr ptr = (NodePtr)kd_res_item_data( presults );
                
                if(ptr->valid == false)
                {   
                    kd_res_next( presults );
                    continue;
                }
                
                NodePtr pre_ptr = ptr->preNode_ptr;
                
                if(pre_ptr == rootNode || ptr == rootNode)
                {
                    kd_res_next( presults );
                    continue;
                }
                
                float update_radius = reCheckRadius(ptr->coord);
                int ret = ReCheckConnection(update_radius, ptr->radius); // -1: not qualified, 0: shrink, 1: no difference, continue
                ptr->radius = update_radius; // the radius of a node may shrink or remain no change, but can not enlarge

                if( ret == -1)
                {   
                    if(ptr->valid == true)
                    {
                        ptr->valid = false;          // Delete this node
                        invalidSet.push_back(ptr);
                        clearBranchS(ptr);
                    }
                }
                else 
                {   
                    float dis = getDis( pre_ptr, ptr );
                    int   res = CheckRelation( dis, pre_ptr, ptr );
                    if( res != -1 ) // ptr and pre_ptr are not connected anymore
                    {   
                        if(pre_ptr->valid == true)
                        {
                              pre_ptr->valid = false;      
                              invalidSet.push_back(pre_ptr);
                              clearBranchS(pre_ptr);
                              kd_res_next( presults );
                              continue;
                        }
                    }

                    vector<NodePtr> childList = ptr->nxtNode_ptr;
                    for(auto childptr: childList)
                    {  // inspect each child of ptr, to see whether they are still connected 
                        float dis = getDis( ptr, childptr );
                        int res = CheckRelation( dis, ptr, childptr );
                        if( res != -1 ) // the child is disconnected with its parent
                        {   
                            if(childptr->valid == true)
                            {
                                childptr->valid = false;      
                                invalidSet.push_back(childptr);
                                clearBranchS(childptr);
                            }
                        } 
                    }
                }

                kd_res_next( presults );
            }

            kd_res_free( presults );
        }
        
        removeInvalid();

        ROS_ERROR("[path repair] Finish repair the path");
        ROS_WARN("[path repair] now %d nodes reach target", int(EndList.size()) );
        ROS_WARN("[path repair] now %d nodes num", int(NodeList.size()) );
}

/* code for repair the tree from the very beginning node in the NodeList
int TreeSize = int(NodeList.size());
int EndListSize = int(EndList.size());
for(int i = 0; i < TreeSize; i ++ )
{     
    ros::Time time_in_loop = ros::Time::now();
    duration = (time_in_loop - time_bef_treeRepair).toSec();

    if( duration > (double)repair_limit ) 
        break;

    NodePtr nodeptr = NodeList[i];
    
        if(nodeptr->valid == false)      
            continue;

        float update_radius = reCheckRadius(nodeptr->coord);
        //float update_radius = reCheckAddRadius(nodeptr->coord);
        if(update_radius< nodeptr->radius)
        {   
            int ret = ReCheckConnection(update_radius, nodeptr->radius); // -1: not qualified, 0: shrink, 1: no difference, continue
            nodeptr->radius = update_radius; 
            if( ret == -1)
            {  
                nodeptr->valid = false;         
                invalidSet.push_back(nodeptr);
                clearBranchS(nodeptr);
            }
            else 
            {   
                vector<NodePtr> childList = nodeptr->nxtNode_ptr;
                for(auto childptr: childList)
                {   
                    if(childptr->valid == true)
                    {
                        float dis = getDis( nodeptr, childptr );
                        int res = CheckRelation( dis, nodeptr, childptr );
                        if( res != -1 ) 
                        {   
                            childptr->valid = false;      
                            invalidSet.push_back(childptr);
                            clearBranchS(childptr);
                        } 
                    }
                }
            }
        }
}
*/

void rrtPathFinder::parentForgetChild(NodePtr ptrParent, NodePtr nodeptr)
{     
      /*ROS_WARN("[parent forget child] start");
      ROS_WARN("[parent forget child] children list size is %d", int(ptrParent->nxtNode_ptr.size()));*/
      vector<NodePtr> child = ptrParent->nxtNode_ptr;
      nodeptr->change = true;  // Let its parent forget it
      ptrParent->nxtNode_ptr.clear();

      for(int i = 0; i < int(child.size()); i++ ){
          NodePtr ptr  = child[i];
          //cout<<ptr->change<<endl;

          if( ptr->change == true ) 
              continue;
          else
              ptrParent->nxtNode_ptr.push_back(ptr);
      }

      nodeptr->change = false; // release the flag after the delete          
}

bool rrtPathFinder::CheckValidEnd(NodePtr endPtr)
{
      NodePtr ptr = endPtr;

      while( ptr != NULL )
      {   
          if(!ptr->valid) 
            return false;
    
          if( getDis( ptr, rootNode->coord ) < ptr->radius )
            return true;

          ptr = ptr->preNode_ptr;
      }

      return false;  

}
// ***************************************************************************************************
// ** Every time the refine function is called, new samples are continuously generated and the tree is continuously rewired, hope to get a better solution.
#if 1
void rrtPathFinder::RRTpathRefine(float refine_limit, bool isPrint)
{     
      //ROS_WARN("[Path Finder] prepare to refine the path");
      double duration = 0.0;
      // /ROS_WARN("[path refine] check the best_distance, %f", best_distance);

      ros::Time time_bef_pathFind = ros::Time::now();

      float pos[3];
      while( true )
      {     
          ros::Time time_in_loop = ros::Time::now();
          duration = (time_in_loop - time_bef_pathFind).toSec();

          if( duration > (double)refine_limit ) 
              break;
          
          Point   pt_sample        =  genSample();
              
          NodePtr node_nearst_ptr = findNearstVertex(pt_sample);
          
          if(!node_nearst_ptr->valid || node_nearst_ptr == NULL ) 
            continue;
                    
          NodePtr node_new_ptr  =  genNewNode(pt_sample, node_nearst_ptr); 
          
          if( node_new_ptr->coord.z < z_l  || node_new_ptr->radius < safety_margin ) 
            continue;
             
          TreeRewire(node_new_ptr, node_nearst_ptr);  
          
          if( node_new_ptr->valid == false ) 
          {
              //ROS_BREAK();
              continue;
          }

          if(check_end( node_new_ptr )){
              
              if( !inform_sample ) 
                  best_end_ptr = node_new_ptr;
              
              EndList.push_back(node_new_ptr);
              UpdateHeuristicRegion(node_new_ptr);
              inform_sample = true;  
          }
          
          pos[0] = node_new_ptr->coord.x;
          pos[1] = node_new_ptr->coord.y;
          pos[2] = node_new_ptr->coord.z;
          kd_insertf(kdTree_, pos, node_new_ptr);

          AddtoGraph(node_new_ptr);
          TreeSparsify(node_new_ptr);      

          if( int(invalidSet.size()) >= cach_size) removeInvalid();
    }

      //ROS_ERROR("[path Finder] check duration %f",duration);
      removeInvalid();
      
      tracePath();

      /*int childNum = 0;
      for(auto nodeptr: NodeList)
        childNum += int(nodeptr->nxtNode_ptr.size());
      ROS_WARN("[path refine] Check children num : %d", childNum );*/

      ROS_WARN("[path refine] Check nodes num : %d", int( NodeList.size() ) );
      //ROS_WARN("[path refine] Check truly feasible path num : %d", int( EndList.size() ) );

/*      int count_null = 0;
      for(auto ptr: NodeList)
        if(ptr->preNode_ptr == NULL)
          count_null ++;
*/
      //ROS_WARN("[path refine] pre null nodes num : %d", count_null);
      //ROS_WARN("[path refine] Finish refine the path");

      ROS_WARN("[path refine] Check feasible path num : %d", int( EndList.size() ) );
 /*     for(auto ptr:EndList)
      {
        cout<<"coordinate: "<<ptr->coord.x<<" , "<<ptr->coord.y<<" , "<<ptr->coord.z<<"R: "<<ptr->radius<<endl;
        cout<<"id: "<<ptr->id<<endl;
      }*/
}
#endif

int rrtPathFinder::Repair()
{
    return -1;
}