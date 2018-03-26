#ifndef _TIME_DATA_TYPE_
#define _TIME_DATA_TYPE_

#include <eigen3/Eigen/Dense>
#define inf 9999999.0
#define _PI M_PI

using namespace std;
using namespace Eigen;

struct Allocator
{     
      double s_step; // step size of the parameter in the virtual domain
      int K;         // discretion number of the grid in s from 0 ~ 1 
      int seg_num;   // segment number of all pieces of trajectory
      double vel_m, acc_m, jer_m; // maximum limit of velocity, acceleration and jerk, in x,y,z axis 

      MatrixXd a, b, c, d; // result of the minimum time optimization program
      MatrixXd time;
      MatrixXd time_acc; // time divided by a in middle of each pair of b
      VectorXd s;

      Allocator( int _seg_num, double _s_step, int _K, double _vel_m, double _acc_m, double _jer_m)
      {     
            seg_num = _seg_num;
            s_step  = _s_step;
            K = _K;

            vel_m = _vel_m;
            acc_m = _acc_m;
            jer_m = _jer_m;

            a.resize(seg_num, K);
            b.resize(seg_num, K + 1);
            c.resize(seg_num, K + 1);
            d.resize(seg_num, K);
            
            time.resize(seg_num, K);
            time_acc.resize(seg_num, K);
            s.resize(K + 1);
      }
      
      Allocator(){}
      ~Allocator(){}
};

struct PhysicalLimiter
{     
      pair<double, double> vel_limit;
      pair<double, double> acc_limit;
      pair<double, double> jer_limit;

      PhysicalLimiter(pair<double, double> _vel_limit, pair<double, double> _acc_limit, pair<double, double> _jer_limit)
      {
            vel_limit = _vel_limit;
            acc_limit = _acc_limit;
            jer_limit = _jer_limit;
      }

      PhysicalLimiter(double max_vel, double max_acc, double max_jerk)
      {
            vel_limit = make_pair( - max_vel,   max_vel);
            acc_limit = make_pair( - max_acc,   max_acc);
            jer_limit = make_pair( - max_jerk,  max_jerk);
      }

      PhysicalLimiter(){};
      ~PhysicalLimiter(){};
};
#endif