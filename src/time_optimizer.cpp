#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include "pcd_trajectory/timeAllocator.h"
#include "pcd_trajectory/time_optimizer.h"
#include "pcd_trajectory/mosek.h"
#include "pcd_trajectory/bezier_base.h"

using namespace std;    
using namespace Eigen;

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
}

MinimumTimeOptimizer::MinimumTimeOptimizer(){};

MinimumTimeOptimizer::~MinimumTimeOptimizer(){};

void MinimumTimeOptimizer::MinimumTimeGeneration( const MatrixXd & polyCoeff, const double & maxVel, const double & maxAcc, const double & maxJer, const int & K, const double & w_a )
{
/* minimum time physical feasible trajectory time allocator. */
/* objective is to generate motion as fast as possible within the physical limitaion (vel, acc and jerk). */
/* TODO List : try monomial-basis polynomial first, the only different to Bezier-basis polynomial is a evaluating function in evaluating 
               the vel and acc in virtual domain*/

/* add constraints on acceleration seems not correct, can not name it "jerk constraint"
    do not want to acceleration change drastically, may add penalty term about acceleration on the objective,
    like Optimal Control */

    //ROS_ERROR("[Time Optimizer] weight_a is %f", w_a);
    _P          = polyCoeff;
    _seg_num    = polyCoeff.rows();
    _poly_num1D = polyCoeff.cols() / 3;

    double d_s = 1.0 / K;
    VectorXd s(K + 1);

    // construct s vector, containing step-increasing discrete s value
    for(int i = 0; i < K + 1; i++)
        s(i) = i * d_s;

    // pre-define all optimized variables and slack variables
    int num_a, num_b, num_c, num_d;
    int num_t, num_e, num_f, num_g, num_h;

    num_a = num_d = K;
    num_b = num_c = K + 1;

    num_t = num_e = K;
    num_f = num_g = num_h = K + 1;

    int num_x = num_a + num_b + num_c + num_d + num_t + num_e + num_f + num_g + num_h;

/*    cout<<"num_x: "<<num_x<<endl;
    cout<<"num_a: "<<num_a<<endl;
    cout<<"num_b: "<<num_b<<endl;*/

    int num_a_n, num_b_n, num_c_n, num_d_n, num_x_n;
    int num_t_n, num_e_n, num_f_n, num_g_n, num_h_n;

    num_a_n = _seg_num * num_a;
    num_b_n = _seg_num * num_b;
    num_c_n = _seg_num * num_c;
    num_d_n = _seg_num * num_d;
    num_t_n = _seg_num * num_t;
    num_e_n = _seg_num * num_e;
    num_f_n = _seg_num * num_f;
    num_g_n = _seg_num * num_g;
    num_h_n = _seg_num * num_h;

    num_x_n = _seg_num * num_x;
    
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    
    int _equ_equal_num = K * _seg_num;      // equality constraints b(i+1) - b(i) == 2 * a(i) * d_s;
    int _equ_conti_num = _seg_num - 1;      // continuity constraints for continuous velocity, note only in one axis continuous means continuous in x,y,z
    int _equ_slack_num = K * _seg_num + 3 * (K + 1) * _seg_num; // map e,f,g,h to original variables

    int _inequ_acc_conti_num = 3 * (_seg_num - 1); // continuity constraints for acceleration by enforcing jerk_max bound deviation constraints

    int _inequ_vel_num = 3 * (K+1) * _seg_num;     // linear bounding constraints for velocity in x, y, z axis
    int _inequ_acc_num = 3 *  K    * _seg_num;     // linear bounding constraints for acceleration
    int _inequ_jer_num = 3 * (K-1) * _seg_num;     // linear bounding constraints for jerk

    int _equ_con_num   = _equ_equal_num + _equ_conti_num; 
    int _inequ_con_num = _inequ_acc_conti_num + _inequ_vel_num + _inequ_acc_num + _inequ_jer_num;
    int _cone_num      = _seg_num * K  + _seg_num * (K + 1);

    int _var_num = num_x_n + 2;
    int _con_num = _equ_con_num + _inequ_con_num + _equ_slack_num;

    double x_var[_var_num];
    MSKrescodee  r; 
    double primalobj;

    for(int i = 0; i < _equ_con_num; i ++)
    {   
        pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair(  0.0, 0.0 ) ); 
        con_bdk.push_back(cb_eq);
    }

    for(int i = 0; i < _equ_slack_num; i ++)
    {   
        pair<MSKboundkeye, pair<double, double> > cb_eq;

        if( i < K * _seg_num + ( K + 1 ) * _seg_num )
            cb_eq = make_pair( MSK_BK_FX, make_pair( 0.0, 0.0 ) ); 
        
        else if( i < K * _seg_num + 2 * ( K + 1 ) * _seg_num)
            cb_eq = make_pair( MSK_BK_FX, make_pair(-1.0, -1.0 ) ); 
        
        else
            cb_eq = make_pair( MSK_BK_FX, make_pair( 1.0, 1.0 ) ); 
        
        //cout<<"cb_eq: "<<cb_eq.second.first<<" , "<<cb_eq.second.second<<endl;
        con_bdk.push_back(cb_eq);
    }

    for(int i = 0; i < _inequ_acc_conti_num; i++)
    {
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxJer, + maxJer ) ); 
        con_bdk.push_back(cb_ie);   
    }

    for(int i = 0; i < _inequ_vel_num; i++)
    {
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( 0.0, maxVel * maxVel ) ); 
        con_bdk.push_back(cb_ie);   
    }

/*    cout<<"_inequ_acc_num: "<<_inequ_acc_num<<endl;
    cout<<"con_bdk size 1: "<<con_bdk.size()<<endl;*/
    for(int i = 0; i < _inequ_acc_num; i++)
    {   
        double lo_bnd;
        double up_bnd;

        if(i < 3 )
        {
            lo_bnd = - maxJer;
            up_bnd = + maxJer;            
        }
        else if( i >= _inequ_acc_num - 3 )
        {
            lo_bnd = - maxJer;
            up_bnd = + maxJer;            
        }
        else
        {
            lo_bnd = - maxAcc;
            up_bnd = + maxAcc;
        }

        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( lo_bnd, up_bnd ) ); 
        con_bdk.push_back(cb_ie);   
    }
    /*cout<<"con_bdk size 1: "<<con_bdk.size()<<endl;

    cout<<"_inequ_acc_num: "<<_inequ_jer_num<<endl;*/
    for(int i = 0; i < _inequ_jer_num; i++)
    {
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxJer, + maxJer ) ); 
        con_bdk.push_back(cb_ie);   
    }
    //cout<<"con_bdk size 1: "<<con_bdk.size()<<endl;

   /* for(auto ptr:con_bdk)
        cout<<ptr.second.first<<" , "<<ptr.second.second<<endl;*/

    /*** ## Stacking bounds for all unknowns ## ***/ 
    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 
    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk_n; 

    // stack all boundary value for variables in one segment
    {
        // for a
        for(int i = 0; i < num_a; i ++ )
        {
            pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); 
            var_bdk.push_back(vb_x);
        }
        // for b
        for(int i = 0; i < num_b; i ++ )
        {
            pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_RA, make_pair( 0.0, + MSK_INFINITY ) ); 
            var_bdk.push_back(vb_x);
        }

        // for c + d
        for(int i = 0; i < num_c + num_d; i ++ )
        {
            pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); 
            var_bdk.push_back(vb_x);
        }
        // for t
        for(int i = 0; i < num_t; i ++ )
        {
            pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FX, make_pair( 1.41421, 1.41421 ) ); 
            var_bdk.push_back(vb_x);
        }

        // for e,f,g,h
        for(int i = 0; i < num_e + num_f + num_g + num_h; i ++ )
        {
            pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); 
            var_bdk.push_back(vb_x);
        }

        for(int k = 0; k < _seg_num; k++)
        {
            for(auto ptr:var_bdk)
                var_bdk_n.push_back(ptr);
        }

        // slacked objective t0
        {
            pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); 
            var_bdk_n.push_back(vb_x);
        }

        // slack variable 1.0 in the rotated objective induced cone
        {
            pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FX, make_pair( 1.0, 1.0 ) ); 
            var_bdk_n.push_back(vb_x);
        }
    }

    MSKint32t  j,i; 
    MSKenv_t   env; 
    MSKtask_t  task; 

    r = MSK_makeenv( &env, NULL ); 
    r = MSK_maketask(env, _con_num, _var_num, &task); 

// Parameters used in the optimizer
//######################################################################
    MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_CONIC );
    //MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 4);
    //MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-5);
/*    MSK_putdouparam (task, MSK_DPAR_INTPNT_CO_TOL_DFEAS,  1e-15);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_CO_TOL_PFEAS,  1e-15);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_CO_TOL_INFEAS, 1e-15);*/

//######################################################################
    
    r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
    
    if ( r == MSK_RES_OK ) 
    {
      r = MSK_appendcons(task, _con_num);  
      r = MSK_appendvars(task, _var_num); 
    }

    //ROS_WARN("[Time Optimizer] Stacking the variable bounds. ");   
    for(j = 0; j<_var_num && r == MSK_RES_OK; ++j){ 
      // Set the bounds on variable j : //  blx[j] <= x_j <= bux[j] 
        if (r == MSK_RES_OK) 
            r = MSK_putvarbound(task, 
                                j,                            // Index of variable. 
                                var_bdk_n[j].first,             // Bound key.
                                var_bdk_n[j].second.first,      // Numerical value of lower bound.
                                var_bdk_n[j].second.second );   // Numerical value of upper bound.      
    } 
    
    //ROS_WARN("[Time Optimizer] Stacking the constraints bounds. ");   
    //   for i=1, ...,con_num : blc[i] <= constraint i <= buc[i] 
    assert(_con_num == (int)con_bdk.size());

    for( i = 0; i < _con_num && r == MSK_RES_OK; i++ ) 
    {
        //cout<<con_bdk[i].second.first<<" , "<<con_bdk[i].second.second<<endl;
        
        r = MSK_putconbound(task, 
                            i,                            // Index of constraint. 
                            con_bdk[i].first,             // Bound key.
                            con_bdk[i].second.first,      // Numerical value of lower bound.
                            con_bdk[i].second.second );   // Numerical value of upper bound. 
    }

    //ROS_WARN("[Time Optimizer] Stacking the equality constraints. ");   
    int row_idx = 0;
    //ROS_WARN("[Time Optimizer] For equality constraints of mapping b_k to a_k");   
    // b_k+1 - b_K - 2 * a_k * d_s == 0.0;
    // ==> - 2 * d_s * a_k - b_K + b_k+1 == 0.0;
    for(int k = 0; k < _seg_num; k ++) // mapping b to a for each segment of the trajectory
    {   
        for(int i = 0; i < K; i++)
        {
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = -2.0 * d_s;
            aval[1] = -1.0;
            aval[2] =  1.0;

            asub[0] = k * num_x + i;
            asub[1] = k * num_x + num_a + i;
            asub[2] = k * num_x + num_a + i + 1;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }
    
    //ROS_WARN("[Time Optimizer] For equality constraints of continuity at velocity"); 
    // continuous in b_k between each consecutive curve 
    for(int k = 1; k < _seg_num; k ++) 
    {       
        int nzi = 2;
        MSKint32t asub[nzi];
        double aval[nzi];
        aval[0] =  1.0;
        aval[1] = -1.0;

        asub[0] = (k - 1) * num_x + num_a + num_b - 1;
        asub[1] =    k    * num_x + num_a;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
        row_idx ++;
    }

    //ROS_WARN("[Time Optimizer] For equality constraints of mapping c_k+1 + c_k to e_k"); 
    // e_k = c_k + c_k+1
    for(int k = 0; k < _seg_num; k ++) 
    {   
        for(int i = 0; i < K; i++ )
        {
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] =  1.0;
            aval[1] =  1.0;
            aval[2] = -1.0;

            asub[0] = k * num_x + num_a + num_b + i;
            asub[1] = k * num_x + num_a + num_b + i + 1;
            asub[2] = k * num_x + num_a + num_b + num_c + num_d + num_t + i;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    //ROS_WARN("[Time Optimizer] For equality constraints of mapping 2*c_k to f_k"); 
    // f_k = 2 * c_k
    for(int k = 0; k < _seg_num; k ++) 
    {   
        for(int i = 0; i < K + 1; i++ )
        {
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] =  2.0;
            aval[1] = -1.0;

            asub[0] = k * num_x + num_a + num_b + i;
            asub[1] = k * num_x + num_a + num_b + num_c + num_d + num_t + num_e + i;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    //ROS_WARN("[Time Optimizer] For equality constraints of mapping g_k to b_k - 1"); 
    // g_k = b_k - 1
    for(int k = 0; k < _seg_num; k ++) 
    {   
        for(int i = 0; i < K + 1; i++ )
        {
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = -1.0;
            aval[1] =  1.0;

            asub[0] = k * num_x + num_a + i;
            asub[1] = k * num_x + num_a + num_b + num_c + num_d + num_t + num_e + num_f + i;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    //ROS_WARN("[Time Optimizer] For equality constraints of mapping h_k to b_k + 1"); 
    // h_k = b_k + 1
    for(int k = 0; k < _seg_num; k ++) 
    {   
        for(int i = 0; i < K + 1; i++ )
        {
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = -1.0;
            aval[1] =  1.0;

            asub[0] = k * num_x + num_a + i;
            asub[1] = k * num_x + num_a + num_b + num_c + num_d + num_t + num_e + num_f + num_g + i;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }
  //  ROS_WARN("[Time Optimizer] For equality constraints of continuity at acceleration, by limiting piecewise deviation"); 
    // continuous in b_k between each consecutive curve 
    for(int k = 1; k < _seg_num; k ++) 
    {   
        double s0 = (s(0) + s(1)) / 2.0;
        double sf = (s(K) + s(K-1)) / 2.0;

/*        cout<<s(K)<< " , "<<s(K-1)<<endl;
        cout<<s0<< " , "<<sf<<endl;*/

        Vector3d f_0 = getVel(k, s0);
        Vector3d f_f = getVel(k-1, sf);

        Vector3d h_0 = getAcc(k, s0);
        Vector3d h_f = getAcc(k-1, sf);

        for(int i = 0; i < 3; i++) //for x, y, and z axis    
        {
            int nzi = 6;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = f_f(i);
            aval[1] = h_f(i) / 2.0;
            aval[2] = h_f(i) / 2.0;

            aval[3] = -f_0(i);
            aval[4] = -h_0(i) / 2.0;
            aval[5] = -h_0(i) / 2.0;

            asub[0] = (k - 1) * num_x + num_a - 1;
            asub[1] = (k - 1) * num_x + num_a + num_b - 2;
            asub[2] = (k - 1) * num_x + num_a + num_b - 1;

            asub[3] = k * num_x;
            asub[4] = k * num_x + num_a;
            asub[5] = k * num_x + num_a + 1;
          /*  cout<<"row_idx: "<<row_idx<<endl;
            cout<<"sub: "<<asub[0]<<" , "<<asub[1]<<" , "<<asub[2]<<" , "<<asub[3]<<" , "<<asub[4]<<" , "<<asub[5]<<endl;*/

            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    //ROS_WARN("[Time Optimizer] For inequality constraints of constraining velocity within limit"); 
    for(int k = 0; k < _seg_num; k ++) 
    {   
        for(int p = 0; p < K + 1; p ++)
        {
            Vector3d f = getVel(k, s(p));

            for(int i = 0; i < 3; i++) //for x, y, and z axis    
            {
                int nzi = 1;
                MSKint32t asub[nzi];
                double aval[nzi];

                aval[0] = pow(f(i), 2);
                
                if( fabs(aval[0]) < 0.0001 )
                    aval[0]  = 0.0001;

                asub[0] = k * num_x + num_a + p;
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
        }
    }
/*
    ROS_WARN("[Time Optimizer] For inequality constraints of constraining acceleration within limit"); 
    cout<<"row_idx: "<<row_idx<<endl;*/
    int acc_ine_num = 0;
    for(int k = 0; k < _seg_num; k ++) 
    {   
        for(int p = 0; p < K; p ++)
        {   
            double s_a = (s(p) + s(p+1)) / 2.0;
            Vector3d f = getVel(k, s_a);
            Vector3d h = getAcc(k, s_a);
            
            for(int i = 0; i < 3; i++) //for x, y and z axis    
            {
                int nzi = 3;
                MSKint32t asub[nzi];
                double aval[nzi];

/*                if( (k == (_seg_num - 1) && p == (K - 1) ) || (k == 0 && p == 0 ) )
                {   
                    if( fabs(f(i)) < 0.001 )
                    {
                        if( f(i) > 0)
                            f(i) = 0.001;
                        else
                            f(i) = - 0.001;
                    }

                    if( fabs(h(i)) < 0.002 )
                    {
                        if( h(i) > 0)
                            h(i) = 0.002;
                        else
                            h(i) = - 0.002;
                    }
                }*/

                aval[0] = f(i);
                aval[1] = h(i) / 2.0;
                aval[2] = h(i) / 2.0;

/*                if(k == (_seg_num - 1) && p == (K - 1) )
                {   
                    ROS_ERROR("Last Acc");
                    cout<<"row_idx: "<<row_idx<<endl;
                    cout<<aval[0]<<" , "<<aval[1]<<" , "<<aval[2]<<endl;
                }*/

                asub[0] = k * num_x + p;
                asub[1] = k * num_x + num_a + p;
                asub[2] = k * num_x + num_a + p + 1;
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
               /* cout<<"row_idx: "<<row_idx<<endl;
                cout<<"sub: "<<asub[0]<<" , "<<asub[1]<<" , "<<asub[2]<<endl;*/
                row_idx ++;
                acc_ine_num++;
            }
        }
    }
    //cout<<"row_idx: "<<row_idx<<endl;
    //ROS_ERROR("[SOLVER] acc limit constains num: %d", acc_ine_num);

#if 1
    //ROS_WARN("[Time Optimizer] For inequality constraints of constraining jerk within limit"); 
    for(int k = 0; k < _seg_num; k ++) 
    {   
        for(int p = 0; p < K - 1; p ++)
        {   
            double s_a_1 = (s(p)   + s(p+1)) / 2.0;
            double s_a_2 = (s(p+1) + s(p+2)) / 2.0;
            //cout<<"s_a_1: "<<s_a_1<<endl;

            Vector3d f_1 = getVel(k, s_a_1);
            Vector3d h_1 = getAcc(k, s_a_1);

            Vector3d f_2 = getVel(k, s_a_2);
            Vector3d h_2 = getAcc(k, s_a_2);

            for(int i = 0; i < 3; i++) //for x, y and z axis    
            {
                int nzi = 5; // 5
                MSKint32t asub[nzi];
                double aval[nzi];

                if( (k == _seg_num - 1 && p == K - 2) || (k == 0 && p == 0 ) )
                {
                    if( fabs(f_2(i)) < 0.0001 )
                    {
                        if(f_2(i) > 0)
                            f_2(i) = 0.0001;
                        else
                            f_2(i) = - 0.0001;
                    }

                    if( fabs(h_2(i)) < 0.0001 )
                    {
                        if(h_2(i) > 0)
                            h_2(i) = 0.0001;
                        else
                            h_2(i) = - 0.0001;
                    }
                }

                aval[0] = f_1(i);
                aval[1] = h_1(i) / 2.0;
                aval[2] = h_1(i) / 2.0 - h_2(i) / 2.0;

                aval[3] = - f_2(i);
                aval[4] = - h_2(i) / 2.0;
                /*aval[4] = - h_2(i) / 2.0;
                aval[5] = - h_2(i) / 2.0;*/

                asub[0] = k * num_x + p;
                asub[1] = k * num_x + num_a + p;
                asub[2] = k * num_x + num_a + p + 1;

                asub[3] = k * num_x + p + 1;
                asub[4] = k * num_x + num_a + p + 2;
                /*asub[4] = k * num_x + num_a + p + 1;
                asub[5] = k * num_x + num_a + p + 2;*/

                /*cout<<"row_idx: "<<row_idx<<endl;
                cout<<"sub: "<<asub[0]<<" , "<<asub[1]<<" , "<<asub[2]<<" , "<<asub[3]<<" , "<<asub[4]<<endl;*/
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
        }
    }
#endif   
    //cout<<"row_idx: "<<row_idx<<endl;

    //ROS_WARN("[Time Optimizer] Stacking all conic cones");
    {   
        for(int k = 0; k < _seg_num; k++)
        {
            // 2 * d_k * e_k >= t_k^2
            for(int i = 0; i < K; i++)
            {   
                int nzi = 3;
                MSKint32t csub[nzi];                
                
                csub[0] = k * num_x + num_a + num_b + num_c + i;
                csub[1] = k * num_x + num_a + num_b + num_c + num_d + num_t + i;

                csub[2] = k * num_x + num_a + num_b + num_c + num_d + i;
                
                r = MSK_appendcone(task, MSK_CT_RQUAD, 0.0, nzi, csub);
            }        

            // h_k >= norm(f_k, g_k)
            for(int i = 0; i < K + 1; i++)
            {   
                int nzi = 3;
                MSKint32t csub[nzi];                
                
                csub[0] = k * num_x + num_a + num_b + num_c + num_d + num_t + num_e + num_f + num_g + i;

                csub[1] = k * num_x + num_a + num_b + num_c + num_d + num_t + num_e + i;
                csub[2] = k * num_x + num_a + num_b + num_c + num_d + num_t + num_e + num_f + i;
                
                r = MSK_appendcone(task, MSK_CT_QUAD, 0.0, nzi, csub);
            }        
        }
    }

    //ROS_WARN("[Time Optimizer] Stacking the quadratic induced rotated cone");
    {           
        int nzi = 2 + num_a_n;
        MSKint32t csub[nzi];                
        
        csub[0] = num_x_n;
        csub[1] = num_x_n + 1;

        int idx = 0;
        for(int k = 0; k < _seg_num; k++)
        {
            for( int i = 0; i < num_a; i++ )
            {
                csub[idx + 2] = k * num_x + i;
                idx ++;
            }

        }

        r = MSK_appendcone(task, MSK_CT_RQUAD, 0.0, nzi, csub);
    }

    //ROS_WARN("[Time Optimizer] Start stacking the objective");
    int nzi = num_d_n + 1;
    MSKint32t asub[nzi];
    double aval[nzi];

    aval[0] = w_a * d_s;
    asub[0] = num_x_n;

    int idx = 1;
    for(int k = 0; k < _seg_num; k++)
    {   
        for(int i = 0; i < K; i++)
        {
            aval[idx] = 2.0 * d_s;
            asub[idx] = k * num_x + num_a + num_b + num_c + i;
            idx ++;
        }
    }
    
    r = MSK_putclist(task, nzi, asub, aval);

    if ( r==MSK_RES_OK ) 
         r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
    
    //ros::Time time_opt = ros::Time::now();
    bool solve_ok = false;
    if ( r==MSK_RES_OK ) 
      { 
        //ROS_WARN("Prepare to solve the problem ");   
        MSKrescodee trmcode; 
        r = MSK_optimizetrm(task,&trmcode); 
        MSK_solutionsummary (task,MSK_STREAM_LOG); 
          
        if ( r==MSK_RES_OK ) 
        { 
          MSKsolstae solsta; 
          MSK_getsolsta (task,MSK_SOL_ITR,&solsta); 
           
          switch(solsta) 
          { 
            case MSK_SOL_STA_OPTIMAL:    
            case MSK_SOL_STA_NEAR_OPTIMAL: 
              
            
            r = MSK_getxx(task, 
                          MSK_SOL_ITR,    // Request the interior solution.  
                          x_var); 

            r = MSK_getprimalobj(
                task,
                MSK_SOL_ITR,
                &primalobj);

            _objective = primalobj;

            solve_ok = true;
            
            break; 
            
            case MSK_SOL_STA_DUAL_INFEAS_CER: 
            case MSK_SOL_STA_PRIM_INFEAS_CER: 
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: 
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:   
              printf("Primal or dual infeasibility certificate found.\n"); 
              break; 
               
            case MSK_SOL_STA_UNKNOWN: 
              printf("The status of the solution could not be determined.\n"); 
              break; 
            default: 
              printf("Other solution status."); 
              break; 
          } 
        } 
        else 
        { 
          printf("Error while optimizing.\n"); 
        } 
      }
     
      if (r != MSK_RES_OK) 
      { 
        // In case of an error print error code and description. 
        char symname[MSK_MAX_STR_LEN]; 
        char desc[MSK_MAX_STR_LEN]; 
         
        printf("An error occurred while optimizing.\n");      
        MSK_getcodedesc (r, 
                         symname, 
                         desc); 
        printf("Error %s - '%s'\n",symname,desc); 
      } 
    
    //ROS_WARN("[Time Optimizer Solver] final time is %f", _objective);
    VectorXd sol(_var_num);
    VectorXd sol_k(num_x);
    VectorXd a_k(num_a);
    VectorXd b_k(num_b);
    VectorXd c_k(num_c);
    VectorXd d_k(num_d);

    VectorXd t_k(num_t);
    VectorXd e_k(num_e);
    VectorXd f_k(num_f);
    VectorXd g_k(num_g);
    VectorXd h_k(num_h);

    for(int i = 0; i < _var_num; i++)
        sol(i) = x_var[i];

    double T = 0.0;
    for(int k = 0; k < _seg_num; k++)
    {
        sol_k = sol.segment(k * num_x, num_x);
        b_k   = sol_k.segment(num_a, num_b);

        for(int i = 0; i < K; i++)
        {
            if( b_k(i) <= 0.0 || b_k(i+1) <= 0.0 )
                continue;

            T += 1.0 * 2 * d_s/(sqrt(b_k(i)) + sqrt(b_k(i+1)));
            /*cout<<"b_k: "<<b_k(i)<<" , b_k_1:"<<b_k(i+1)<<endl;
            cout<<(sqrt(b_k(i)) + sqrt(b_k(i+1)))<<endl;
            cout<<"T: "<<T<<endl;*/
        }
    }

    ROS_WARN("[Time Optimizer Solver] final time is %f", T);

/*    
    
    ROS_WARN("check final Minimum Time : %f", T);
*/ 
    //ROS_WARN("check velocity ");
    /*for(int k = 0; k < _seg_num; k++)
    {   
        sol_k = sol.segment(k * num_x, num_x);
        b_k   = sol_k.segment(num_a, num_b);
        for(int i = 0; i < K + 1; i++)
        {
            Vector3d vel_s = getVelPoly(k, s(i));
            Vector3d vel_t = sqrt(b_k(i)) * vel_s;
        
            cout<<vel_t(0)<<" , "<<vel_t(1)<<" , "<<vel_t(2)<<endl;
        }
    }*/

    // double _s_step, int _K, double _vel_m, double _acc_m, double _jer_m
    time_allocator = new Allocator(_seg_num, d_s, K, maxVel, maxAcc, maxJer);
    for(int k = 0; k < _seg_num; k++)
    {
        double T     = 0.0;
        //time_allocator->time(k, 0) = T;

        sol_k = sol.segment(k * num_x, num_x);
        
        a_k   = sol_k.segment(0, num_a);
        b_k   = sol_k.segment(num_a, num_b);
        c_k   = sol_k.segment(num_a + num_b, num_c);
        d_k   = sol_k.segment(num_a + num_b + num_c, num_d);

        t_k   = sol_k.segment(num_a + num_b + num_c + num_d, num_t);
        
        e_k   = sol_k.segment(num_a + num_b + num_c + num_d + num_t, num_e);
        f_k   = sol_k.segment(num_a + num_b + num_c + num_d + num_t + num_e, num_f);
        g_k   = sol_k.segment(num_a + num_b + num_c + num_d + num_t + num_e + num_f, num_g);
        h_k   = sol_k.segment(num_a + num_b + num_c + num_d + num_t + num_e + num_f + num_g, num_h);

        for(int i = 0; i < K + 1; i++)
        {
            if(i <  K)
            {
                time_allocator->a(k, i) = a_k(i);
                time_allocator->d(k, i) = d_k(i);

                if( b_k(i) <= 0.0 || b_k(i+1) <= 0.0 )
                    T += 0.0;
                else
                    T += 1.0 * 2 * d_s/(sqrt(b_k(i)) + sqrt(b_k(i+1)));
                
                time_allocator->time(k, i) = T;
                
                if(i == 0)
                    time_allocator->time_acc(k, i) = time_allocator->time(k, i) / 2.0;
                else
                    time_allocator->time_acc(k, i) = (time_allocator->time(k, i) + time_allocator->time(k, i - 1)) / 2.0;
            }
            
            time_allocator->b(k, i) = b_k(i);
            time_allocator->c(k, i) = c_k(i);
            time_allocator->s(i) = s(i);
        }

    }


    MSK_deletetask(&task); 
    MSK_deleteenv(&env); 
}

Vector3d MinimumTimeOptimizer::getVel(int k, double s)
{
    // if(_type == 0)
        return getVelPoly(k ,s);
    // else
        // return getVelBezier(k ,s);
}

Vector3d MinimumTimeOptimizer::getAcc(int k, double s)
{
    // if(_type == 0)
        return getAccPoly(k ,s);
    // else
        // return getAccBezier(k ,s);
}

Vector3d MinimumTimeOptimizer::getVelPoly(int k, double s)
{
    Vector3d ret;

    for ( int dim = 0; dim < 3; dim++ )
    {
        VectorXd coeff = (_P.row(k)).segment( dim * _poly_num1D, _poly_num1D );
        VectorXd t = VectorXd::Zero( _poly_num1D );
        
        for(int j = 0; j < _poly_num1D; j ++)
            if(j==0)
                t(j) = 0.0;
            else
                t(j) = j * pow(s, j-1);

        ret(dim) = coeff.dot(t);
    }

    return ret;
}

Vector3d MinimumTimeOptimizer::getAccPoly(int k, double s)
{
    Vector3d ret;

    for ( int dim = 0; dim < 3; dim++ )
    {
        VectorXd coeff = (_P.row(k)).segment( dim * _poly_num1D, _poly_num1D );
        VectorXd t = VectorXd::Zero( _poly_num1D );

        for(int j = 0; j < _poly_num1D; j ++)
            if( j==0 || j==1 )
                t(j) = 0.0;
            else
                t(j) = j * (j - 1) * pow(s, j-2);

        ret(dim) = coeff.dot(t);
    }

    return ret;
}