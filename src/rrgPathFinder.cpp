#include "pcd_trajectory/rrgPathFinder.h"

rrgPathFinder::rrgPathFinder(double xl, double xh, double yl, double yh, double zl, double zh, double biasl, double biash)
{     
      eng = default_random_engine(rd()) ;
      x_l = xl;
      x_h = xh;
      y_l = yl;
      y_h = yh;
      z_l = zl;
      z_h = zh;
      bias_l = biasl;
      bias_h = biash;

      path_find_state  = true;
      
      rand_x = uniform_real_distribution<double>(x_l, x_h);
      rand_y = uniform_real_distribution<double>(y_l, y_h);
      rand_z = uniform_real_distribution<double>(z_l, z_h);
      rand_bias = uniform_real_distribution<double>(bias_l, bias_h);

      /*kdtreeForMap = new pcl::search::KdTree<pcl::PointXYZ> ();
      CloudIn = new pcl::PointCloud<pcl::PointXYZ>();*/
}

rrgPathFinder::~rrgPathFinder(){ }

void rrgPathFinder::reset()
{
      path_find_state = true;
      NodeList.clear();
      Clear_NodeMap();
      EndList.clear();
      pointIdxRadiusSearch.clear();
      pointRadiusSquaredDistance.clear();        

      rootPtr = NULL;
      p_Astar_target = NULL;
      p_Astar_start = NULL;
}

void rrgPathFinder::setPt( Point startPt, Point endPt)
{
      start_pt = startPt;
      end_pt   = endPt; 
      
      //ROS_WARN("check start pt");
      
      cout<<start_pt.x<<endl;
      cout<<start_pt.y<<endl;
      cout<<start_pt.z<<endl;

      //ROS_WARN("check end pt");
      
      cout<<end_pt.x<<endl;
      cout<<end_pt.y<<endl;
      cout<<end_pt.z<<endl;
}

void rrgPathFinder::setInput(pcl::PointCloud<pcl::PointXYZ> CloudIn)
{     
      pcl::VoxelGrid<pcl::PointXYZ>  VoxelSampler;
      pcl::PointCloud<pcl::PointXYZ> Cloud_DS;
      pcl::PointCloud<pcl::PointXYZ> Cloud_DS_NaN;
      
      ros::Time time1 = ros::Time::now();
      VoxelSampler.setLeafSize(0.1f, 0.1f, 0.1f);
      VoxelSampler.setInputCloud( CloudIn.makeShared() );      
      VoxelSampler.filter( Cloud_DS );    

      std::vector<int> indices;
      pcl::removeNaNFromPointCloud(Cloud_DS,Cloud_DS_NaN, indices);

      ROS_WARN("Input Set OK .. ");
      kdtreeForMap.setInputCloud( Cloud_DS_NaN.makeShared() );    
      ROS_WARN("Check points number of pointCloud : %lu", Cloud_DS_NaN.points.size());
      ros::Time time2 = ros::Time::now();
      
      cout<<"time input is: "<<time2-time1<<endl;
}

double rrgPathFinder::getDis(const NodePtr node1, const NodePtr node2){
      return sqrt(pow(node1->coord.x - node2->coord.x, 2) + pow(node1->coord.y - node2->coord.y, 2) + pow(node1->coord.z - node2->coord.z, 2) );
}

double rrgPathFinder::getDis(const NodePtr node1, const Point pt){
      return sqrt(pow(node1->coord.x - pt.x, 2) + pow(node1->coord.y - pt.y, 2) + pow(node1->coord.z - pt.z, 2) );
}

double rrgPathFinder::getDisL1(const NodePtr node1, const Point pt){
      return fabs( node1->coord.x - pt.x) + fabs(node1->coord.y - pt.y) + fabs(node1->coord.z - pt.z) ;
}


double rrgPathFinder::radius_search( Point search_Pt)//, pcl::search::KdTree<pcl::PointXYZ>::Ptr Tree)
{    
      pcl::PointXYZ searchPoint;
      searchPoint.x = search_Pt.x;
      searchPoint.y = search_Pt.y;
      searchPoint.z = search_Pt.z;

      kdtreeForMap.nearestKSearch(searchPoint, 1, pointIdxRadiusSearch, pointRadiusSquaredDistance);
      double radius = sqrt(pointRadiusSquaredDistance[0]); // The height of flight yard is about 2.6 m
     // double rad = min(radius, 10.0);
      return radius - 0.15;
}

Point rrgPathFinder::genSample()
{     
      double x, y, z, bias; 
      x    = rand_x(eng);
      y    = rand_y(eng);
      z    = rand_z(eng);  
      bias = rand_bias(eng);

      Point pt;

      pt.x = x;
      pt.y = y;
      pt.z = z;

      if( bias > 0.50 && bias < 0.525 )
          pt = end_pt;

      return pt;

}

NodePtr rrgPathFinder::findNearst( Point pt_sample )
{
      NodePtr node_nearst_ptr = NULL;
      
      double min_dis = inf;
      double distance;
      
      int row = int( pt_sample.x * 0.1 );
      int col = int( pt_sample.y * 0.1 );

      if(row <= 0) row += 6;
      else row += 7;

      if(col <= 0) col += 6;
      else col += 7;

    for(int i = -1; i < 2; i++ )
      for(int j = -1; j < 2; j++ ){
          
          int _row = min( max(0, row + i), 13 );
          int _col = min( max(0, col + j), 13 );
          
          if(NodeMap[_row][_col].empty()) continue;
          
          for(vector<NodePtr>::iterator iter = NodeMap[_row][_col].begin(); iter != NodeMap[_row][_col].end(); ++iter){   
              distance = getDis(*iter, pt_sample);
              if(distance < min_dis){ 
                  min_dis = distance;
                  node_nearst_ptr = *iter;     
              }      
          }
      }
      return node_nearst_ptr;
}

NodePtr rrgPathFinder::findNearst( Point pt_sample, vector<NodePtr> NodeList )
{
      NodePtr node_nearst_ptr = NULL;
      double min_dis = inf;
      double distance;
      
      for(vector<NodePtr>::iterator iter = NodeList.begin(); iter != NodeList.end(); ++iter){   

          distance = getDis(*iter, pt_sample);// -  (*iter) ->radius;
          if(distance < min_dis){ 
              min_dis = distance;
              node_nearst_ptr = *iter;     
          }      
      }

      return node_nearst_ptr;
}

NodePtr rrgPathFinder::genNewNode( Point pt_sample, NodePtr node_nearst_ptr )
{
      NodePtr node_new_ptr = new Node;
      double dis = getDis(node_nearst_ptr, pt_sample);
      Point node_new_center;
      node_new_center.x = node_nearst_ptr->coord.x + (pt_sample.x - node_nearst_ptr->coord.x) * node_nearst_ptr->radius / dis;
      node_new_center.y = node_nearst_ptr->coord.y + (pt_sample.y - node_nearst_ptr->coord.y) * node_nearst_ptr->radius / dis;
      node_new_center.z = node_nearst_ptr->coord.z + (pt_sample.z - node_nearst_ptr->coord.z) * node_nearst_ptr->radius / dis;

      node_new_ptr->coord = node_new_center;
      node_new_ptr->id    = -1;

      node_new_ptr->f     = inf; 
      node_new_ptr->g     = inf; 
      node_new_ptr->h     = inf; 
      
      node_new_ptr->preNode_ptr = NULL;

      return node_new_ptr;
}

bool rrgPathFinder::check_end( NodePtr ptr )
{    
      double distance = getDis(ptr, end_pt);
      
      if(distance + 0.1 < ptr->radius)
          return true;      
     
     return false;
}

bool rrgPathFinder::check_end_AStar(NodePtr ptr)
{
      double distance = getDis(ptr, start_pt);

      if(distance + 0.1 < ptr->radius)
          return true;

      return false;
}

void rrgPathFinder::Add_new_node(NodePtr ptr_node_new)
{   
    NodePtr ptr = ptr_node_new;
    int row = int(ptr->coord.x * 0.1 );
    int col = int(ptr->coord.y * 0.1 );
    
    if(row <= 0) row += 6;
    else row += 7;

    if(col <= 0) col += 6;
    else col += 7;
    
    NodeMap[row][col].push_back(ptr);

}

void rrgPathFinder::Clear_NodeMap()
{
    for(int i = 0; i < 14; i++ )
      for(int j =0; j < 14; j++ )
        NodeMap[i][j].clear();
}

bool rrgPathFinder::check_no_Empty(Point pt_sample)
{   
      int row = int( pt_sample.x * 0.1 );
      int col = int( pt_sample.y * 0.1 );

      if(row <= 0) row += 6;
      else row += 7;

      if(col <= 0) col += 6;
      else col += 7;

    for(int i = -1; i < 2; i++)
      for(int j = -1; j < 2; j++){   
          int _row = min( max(0, row + i), 13 );
          int _col = min( max(0, col + j), 13 );
          
          if( !NodeMap[_row][_col].empty() )
              return true;
      }

    return false;
    //
}

bool rrgPathFinder::CheckConnect( double dis, NodePtr node_1, NodePtr node_2 )
{     
      for(vector<NodePtr>::iterator iter = node_2->father_ptr.begin(); iter != node_2->father_ptr.end(); ++iter){   
          if( (*iter)->radius == node_1->radius && (*iter)->coord.x == node_1->coord.x && (*iter)->coord.y == node_1->coord.y && (*iter)->coord.z == node_1->coord.z )         
              return false;
      }

/*      if( dis < sqrt(abs(node_1->radius * node_1->radius - node_2->radius * node_2->radius)) + 1.0 )
          return false;
      
      if((dis + 1.0) < (node_1->radius + node_2->radius))
          return true;
*/
/*      if( dis < (1.2 * sqrt(abs(node_1->radius * node_1->radius - node_2->radius * node_2->radius)) + 0.25) )//0.25 ?
          return false;*/
      
      if( (dis + 0.3 < (node_1->radius + node_2->radius)) && (dis > node_1->radius) && (dis > node_2->radius) )
          return true;

      return false;
}

void rrgPathFinder::AddtoGraph( NodePtr node_new_ptr, NodePtr node_nearst_ptr)
{ 
      
      NodePtr ptr     = node_nearst_ptr;
      NodePtr ptr_new = node_new_ptr;      
      
      ptr_new->father_ptr.push_back( ptr );
      ptr->father_ptr.push_back( ptr_new );

      int row = int( ptr_new->coord.x * 0.1 );
      int col = int( ptr_new->coord.y * 0.1 );

      if(row <= 0) row += 6;
      else row += 7;

      if(col <= 0) col += 6;
      else col += 7;

    /*for(int i = -2; i < 3; i++ )
      for(int j = -2; j < 3; j++ ){*/
    int sub = ceil(ptr_new->radius * 0.1) + 1;
    int add = ceil(ptr_new->radius * 0.1) + 2;

/*    cout<<"sub: "<<sub<<endl;
    cout<<"add: "<<add<<endl;*/

    for(int i = -sub; i < add; i++ )
      for(int j = -sub; j < add; j++ ){
          int _row = min( max(0, row + i), 13 );
          int _col = min( max(0, col + j), 13 );

          if(NodeMap[_row][_col].empty())
              continue;
          
          for(vector<NodePtr>::iterator iter = NodeMap[_row][_col].begin(); iter != NodeMap[_row][_col].end(); ++iter){   
              double dis = getDis( *iter, ptr_new );
              if( CheckConnect( dis, *iter, ptr_new ) )
                  ptr_new->father_ptr.push_back( *iter );

              if( CheckConnect( dis, ptr_new, *iter ) )
                  (*iter)->father_ptr.push_back( ptr_new );
          }
      }

}

void rrgPathFinder::RRGpathFind()
{
      Node root;
      root.dis.clear();
      root.coord = start_pt;
      root.preNode_ptr = NULL;
      root.radius = radius_search( root.coord );//, kdtreeForMap );
      root.id = -1;
      root.f = inf;
      root.g = inf;
      root.h = inf;

      rootPtr = &root;

      NodeList.push_back(&root);
      Add_new_node(&root);

      double accu_time_find_nearst = 0.0;
      double accu_time_radius_search = 0.0;
      int iter_count;

      for( iter_count = 0; iter_count < 30000; iter_count ++)
      {     
         
         ros::Time time_bef_findNearest = ros::Time::now();
         
         Point   pt_sample        =  genSample();
         NodePtr node_nearst_ptr;

         if(check_no_Empty(pt_sample))
            node_nearst_ptr   =  findNearst(pt_sample);
         else
            node_nearst_ptr  =  findNearst(pt_sample, NodeList);
         
         ros::Time time_aft_findNearest = ros::Time::now();
         
         accu_time_find_nearst += (time_aft_findNearest - time_bef_findNearest).toSec();

         NodePtr node_new_ptr     =  genNewNode(pt_sample, node_nearst_ptr);

         ros::Time time_bef_radiusSearch = ros::Time::now();

         node_new_ptr->radius     =  radius_search( node_new_ptr->coord );//, kdtreeForMap );
         
         ros::Time time_aft_radiusSearch = ros::Time::now();

         accu_time_radius_search += (time_aft_radiusSearch - time_bef_radiusSearch).toSec();
         

         if(fabs(node_new_ptr->coord.x) + fabs(node_new_ptr->radius) > 40.0 || fabs(node_new_ptr->coord.y) + fabs(node_new_ptr->radius) > 40.0)
            continue;
         
         //if(node_new_ptr->radius < 1.5)
         if( node_new_ptr->radius < 1.0 )
            continue;
         
         AddtoGraph(node_new_ptr, node_nearst_ptr);

         if(check_end( node_new_ptr ))
            {
              EndList.push_back(node_new_ptr);
              Add_new_node( node_new_ptr );
              NodeList.push_back(node_new_ptr);
              
              NodeList.clear();
              NodeList.push_back(&root);
              Clear_NodeMap();
              Add_new_node(&root);

            }
          else{   
              Add_new_node( node_new_ptr );
              NodeList.push_back(node_new_ptr);
          }
      }

      ROS_WARN("Check feasible path num : %d", int( EndList.size() ) );
      
      if( int( EndList.size() ) < 1){
        ROS_WARN("can't find a feasible path");
        path_find_state = false;
        return;
      }

      /*ROS_WARN("check accu nearesrt find time is ");
      cout<<accu_time_find_nearst<<endl;
      ROS_WARN("check accu radius search time is ");
      cout<<accu_time_radius_search<<endl;
      */
      /*cout<<"NodeList size : "<<NodeList.size()<<endl;*/
      
      ROS_WARN("iter_count = %d",iter_count);
      NodePtr ptr_for_Astar = EndList[EndList.size() - 1];

/*      ROS_WARN(" What the fuck ?");
      for(auto ptr: EndList)
        cout<<(ptr->father_ptr).size()<<endl;*/

      double min_dist = inf;
      for(vector<NodePtr>::iterator iter = EndList.begin(); iter != EndList.end() - 1; ++iter){
        
        if( ((*iter)->father_ptr).empty())
          continue;
        
        double dist = getDis( (*iter), start_pt);
        if( dist < min_dist ){
              min_dist = dist;
              ptr_for_Astar = *iter;
          }
      }        

      for(vector<NodePtr>::iterator iter = EndList.begin(); iter != EndList.end(); ++iter)
         if( ptr_for_Astar != (*iter) )
            ptr_for_Astar->father_ptr.push_back(*iter);

      p_Astar_target = ptr_for_Astar;
      
      ROS_WARN(" Time Consume in A Star ");
      ros::Time A_star_beg = ros::Time::now();
      AStarPath( &root, p_Astar_target);
      ros::Time A_star_end = ros::Time::now();
      cout<<A_star_end - A_star_beg<<endl;

      NodePtr ptr = p_Astar_start;
      int idx = 0;
      while( ptr != NULL) 
      {     
            if(getDis( ptr, end_pt) < ptr->radius ){
               ptr->preNode_ptr = NULL;
               idx ++;
               break;
            }

            ptr = ptr->preNode_ptr;
            idx ++;
      }
      
      Path.resize(idx + 2, 3) ;
      Radius.resize( idx );
      idx = 0;
      
      ptr = p_Astar_start;
      Path.block< 1, 3 >( idx, 0 ) << start_pt.x, start_pt.y, start_pt.z;
      
      idx ++ ;
      
      while( ptr != NULL) {      
            Path.block< 1, 3 >( idx, 0 ) << ptr->coord.x, ptr->coord.y, ptr->coord.z;
            Radius[ idx - 1 ] = ptr->radius;
            
            ptr = ptr->preNode_ptr;
            idx ++;
      }

      Path.block< 1, 3 >( idx, 0 ) << end_pt.x, end_pt.y, end_pt.z;
      
      /*cout <<  "Path : "  <<Path.rows()   << "\n"  <<  Path   << endl;
      cout <<  "Radius : "<<Radius.size() << "\n"  <<  Radius << endl;    */

      double last_x   = Path( idx - 1, 0 );
      double last_y   = Path( idx - 1, 1 );
      double last_z   = Path( idx - 1, 2 );

      double last_rad = Radius(Radius.size() - 1 ); 

      /*cout<<"last ball center : "<<last_x<<" , "<<last_y<<" , "<<last_z<<endl;
      cout<<"end point : "<<end_pt.x<<" , "<<end_pt.y<<" , "<<end_pt.z<<endl;
      cout<<"last_rad: "<<last_rad<<endl;*/

      if( (pow( last_x - end_pt.x, 2 ) + pow( last_y - end_pt.y, 2 ) + pow( last_z - end_pt.z, 2 ) ) > last_rad * last_rad ){
        ROS_WARN("last ball path did't hold the end_point");
        path_find_state = false;
      }

      NodeList.clear();
      Clear_NodeMap();
      EndList.clear();
}

void rrgPathFinder::AStarPath(NodePtr p_start, NodePtr p_target)
{       
        // -1 for notKnown, 1 for OpenList, 0 for CloseList
        multimap<double, NodePtr> weight_node;         
        NodePtr u = p_target;
      
        u->id = 1;
        u->h = getDis( p_start, p_target);
        u->g  = 0.0; 
        u->f  = u->g + u->h; 
        weight_node.insert(make_pair( u->f, u));

        //Update ...
        auto update = [&](NodePtr p_now, NodePtr p_pre) -> void
        {
            if (p_now == p_pre) return ;
            double cost = getDis(p_now, p_pre);

            if (p_pre->id == -1){   
                p_pre->h = getDis(p_pre, p_start);
                p_pre->g = p_now->g + cost;
                p_pre->f = p_pre->g + p_pre->h;

                p_pre->preNode_ptr = p_now;

                p_pre->id = 1;
                weight_node.insert(make_pair(p_pre->f, p_pre));
            
            } else if ( cost + p_now->g < p_pre->g ){
                    p_pre->g  = cost + p_now->g;
                    p_pre->id = 1;
                    p_pre->preNode_ptr = p_now;
                
                    weight_node.insert(make_pair(p_pre->f, p_pre));
                }
        };

         // A * begin here
        while ( !weight_node.empty() )
        {   
            u = weight_node.begin() -> second;
            weight_node.erase(weight_node.begin());

            if ( !( u->father_ptr.empty() ) ){
                for (auto & pr: u->father_ptr ){
                      update( u, pr );
                  }
            }
            u->id = 0;

            if( check_end_AStar(u) ){     
                  p_Astar_start = u;
                  break;
            }
            u->id = 0;
        }

}

Eigen::MatrixXd rrgPathFinder::getPath()
{
    return Path;
}

Eigen::VectorXd rrgPathFinder::getRadius()
{
    return Radius;
}