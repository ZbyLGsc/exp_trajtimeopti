#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include "pcd_trajectory/trajectory_generator_socp_lite.h"
#include "pcd_trajectory/mosek.h"
#include "pcd_trajectory/bezier_base.h"

using namespace std;    
using namespace Eigen;

#define inf 1>>30

typedef SparseMatrix<double> SMatrixXd;

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
} /* printstr */


TrajectoryGeneratorSOCP::TrajectoryGeneratorSOCP(){}

TrajectoryGeneratorSOCP::~TrajectoryGeneratorSOCP(){}

MatrixXd TrajectoryGeneratorSOCP::BezierPloyCoeffGenerationSOCP(
            const MatrixXd  &Path,
            const VectorXd  &Radius,
            const VectorXd  &Time,
            const vector<int>      &multi_poly_order,
            const vector<MatrixXd> &FMList,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int minimize_order )  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   

#define ENFORCE_VEL  1 // whether or not adding extra constraints for ensuring the velocity feasibility
#define ENFORCE_ACC  1 // whether or not adding extra constraints for ensuring the velocity feasibility
#define ENFORCE_Z_CONSTRAIN 0 // whether ot not adding explicitly constraints for z > 0

    //ROS_WARN("[Bezier Trajectory] In TrajectoryGeneratorSOCP");
    ros::Time time_start = ros::Time::now();
    MatrixXd PolyCoeff;

    _Time   = Time;
    _Radius = Radius;      
    _Scale = _Time;
    
    initScale = _Scale(0);
    lstScale  = _Scale(_Scale.size() - 1);

    assert(_Time.size() == _Radius.size() );
    _Path = Path.block(1, 0, Path.rows() - 2, 3 );
    
    int _segment_num   = _Time.size();

    Vector3d StartPt  = pos.row(0); 
    Vector3d EndPt    = pos.row(1); 

    Vector3d StartVel = vel.row(0); 
    Vector3d EndVel   = vel.row(1);

    Vector3d StartAcc = acc.row(0);
    Vector3d EndAcc   = acc.row(1);

/*    cout<<"StartVel: \n"<<StartVel<<endl;
    cout<<"StartAcc: \n"<<StartAcc<<endl;*/
//#####################################################################################################################
//Prepare for constaints and objective data stucture for Mosek solver .         

    int _s1d1CtrlP_num, _s1CtrlP_num;

    int equ_con_s_num = 3 * 3; // p, v, a in x, y, z axis at the start point
    int equ_con_e_num = 3 * 3; // p, v, a in x, y, z axis at the end point
    int equ_con_continuity_num = 3 * 3 * (_segment_num - 1);
    int _equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position

    int _inequ_con_num = 0;
    int _vel_con_num   = 0; 
    int _acc_con_num   = 0; 
    int _CtrlP_num     = 0;
    int _obj_nzero_num = 0;

    cout<<multi_poly_order.size()<<endl;
    for(int i = 0; i < _segment_num; i++)
    {   
        int traj_order       = multi_poly_order[i];
        _s1d1CtrlP_num       = traj_order + 1;
        _inequ_con_num      += _s1d1CtrlP_num;
        _CtrlP_num          += _s1d1CtrlP_num * 3;
        _obj_nzero_num      +=(_s1d1CtrlP_num - minimize_order) * 3; 
        
        //if( i > 0 ){ // not enforced at the 1st segment
            _vel_con_num +=  traj_order * 3;
            _acc_con_num += (traj_order - 1) * 3;
        //}
    }

    if( !ENFORCE_VEL )
        _vel_con_num = 0;

    if( !ENFORCE_ACC )
        _acc_con_num = 0;

    int _high_order_con_num = _vel_con_num + _acc_con_num; 

    int _vel_var_num = _vel_con_num;
    int _acc_var_num = _acc_con_num;
    int _high_order_var_num = _vel_var_num + _acc_var_num; 

    // additional equality constraints introduced by the reformulation of all quadratic constraints
    // each quadratic ==> 1 for u - c'x = w, 3 for y = Fx
    // the objective  ==> _obj_nzero_num for y = Fx
    //int _equ_con_extra_num = _obj_nzero_num; 
    int _equ_con_extra_num = _inequ_con_num + _inequ_con_num * 3 + _obj_nzero_num; 

    // additional functional variables introduced by the reformulation of all quadratic constraints
    // each quadratic ==> 1 for t, 1 for w, 3 for y = Fx ... objective ==> 1 w, several y
    // we should calculate non-zeros in FM, for introducing extra y variables
    int _var_extra_obj_num = 1              + 1              + _obj_nzero_num;
    int _var_extra_qua_num = _inequ_con_num + _inequ_con_num + _inequ_con_num * 3;
    int _var_extra_num     = _var_extra_qua_num + _var_extra_obj_num; 

    int _var_y_con = _inequ_con_num * 3;
    int _var_y_obj = _obj_nzero_num;
    int _var_y_num = _var_y_obj + _var_y_con;

    int _var_w_num = 1 + _inequ_con_num; // w for replacing all quadratic terms by linear terms
    int _var_t_num = 1 + _inequ_con_num; // t= 1, used in the conic cones

    assert( _var_extra_num == _var_t_num + _var_y_num + _var_w_num );
    assert( _inequ_con_num * 3 == _CtrlP_num );

    int _con_num = _equ_con_num + _equ_con_extra_num + _high_order_con_num;
    int _var_num = _CtrlP_num   + _var_extra_num     + _high_order_var_num; // another way is by introducing extra linear qualities and middle variables

/*    cout<<"_CtrlP_num: "<<_CtrlP_num<<endl;
    cout<<"_var_extra_num: "<<_var_extra_num<<endl;
    cout<<"_var_w_num: "<<_var_w_num<<endl;
    cout<<"_var_y_num: "<<_var_y_num<<endl;
    cout<<"_var_t_num: "<<_var_t_num<<endl;*/

    double x_var[_var_num];
    MSKrescodee  r; 
    double primalobj;
    /* ## define a container for constaints boundary and boundkey ## */ 
    /* ## dataType in the double pair is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    //ROS_WARN("[Bezier Trajectory] Start stacking the bounding value for constraints");
    /*** Here is most important stuff, we will elliminate all quadratical constraints and this would introduce many equality constraints and more additional variables***/
    /*************************************************************************************************/
    /** Now firstly we omit the reformulation of the quadratic part in the objective **/
    //ROS_WARN(" ... bounding value for original equality induced value");
    for(int i = 0; i < _equ_con_num; i ++ ){ 
        double beq_i;
        if(i < 3)                    beq_i = StartPt(i); 
        else if (i >= 3  && i < 6  ) beq_i = StartVel(i - 3); 
        else if (i >= 6  && i < 9  ) beq_i = StartAcc(i - 6);
        else if (i >= 9  && i < 12 ) beq_i = EndPt (i - 9 );
        else if (i >= 12 && i < 15 ) beq_i = EndVel(i - 12);
        else if (i >= 15 && i < 18 ) beq_i = EndAcc(i - 15);
        else beq_i = 0.0;

        pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
        con_bdk.push_back(cb_eq);
    }

    //ROS_WARN(" ... bounding value for corridor induced value");
    /***  Stack the bounding value for equality constraints induced by the corridor constraints  ***/
    for(int k = 0; k < _segment_num; k++ )
    {
        _s1d1CtrlP_num = multi_poly_order[k] + 1;
        for(int i = 0; i < _s1d1CtrlP_num; i++)
        {
            double bin_i = _Radius(k) * _Radius(k) - _Path(k, 0) * _Path(k, 0) - _Path(k, 1) * _Path(k, 1) - _Path(k, 2) * _Path(k, 2) ;            
            pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_FX, make_pair( bin_i, bin_i ) ); 
            con_bdk.push_back(cb_ie);   
        }
    }

    //ROS_WARN(" ... bounding value for mapping x to y induced value");
    /***  Stack the bounding value for mapping Fx to additional variables y  ***/
    for(int i = 0; i < _var_y_num; i++)
    {
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_FX, make_pair( 0.0, 0.0 ) ); // # cb_ie means: constriants boundary of inequality constrain      
        con_bdk.push_back(cb_ie);   
    }

if(ENFORCE_VEL)
{
    /***  Stack the bounding value for the linear inequality for the velocity constraints  ***/
    // for velocity constraints
    for(int i = 0; i < _vel_con_num; i++)
    {
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_FX, make_pair( 0.0, 0.0) );
        con_bdk.push_back(cb_ie);   
    }
}

if(ENFORCE_ACC)
{
    /***  Stack the bounding value for the linear inequality for the velocity constraints  ***/
    // for velocity constraints
    for(int i = 0; i < _acc_con_num; i++)
    {
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_FX, make_pair( 0.0, 0.0) );
        con_bdk.push_back(cb_ie);   
    }
}

    /*** ## Stacking bounds for all unknowns ## ***/ 
    /*** The sequence is control points + additional variables: x, w, y, t ***/
    //ROS_WARN("[Bezier Trajectory] Start stacking the bounding value for variables");
    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 
    for(int i = 0; i < _CtrlP_num; i ++ ){
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); 
        var_bdk.push_back(vb_x);
    } 

if(ENFORCE_Z_CONSTRAIN)
{
    for(int k = 0; k < _segment_num ; k ++ )
    {   
        int order = multi_poly_order[k];
        _s1d1CtrlP_num = order + 1;
        for(int i = 0; i < 3; i++)
        {   
            pair<MSKboundkeye, pair<double, double> > vb_x;
            
            if( i == 2)
                vb_x  = make_pair( MSK_BK_FR, make_pair( 0.0, + MSK_INFINITY ) ); 
            else
                vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); 
            
            for(int p = 0; p < _s1d1CtrlP_num; p++)
                var_bdk.push_back(vb_x);
        }
    }
}

    /* ## Variable bounds for addtional variables y and w ## */ 
    for(int i = 0; i < (_var_y_num + _var_w_num); i ++ ){
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) );
        var_bdk.push_back(vb_x);
    }

    /* ## Variable bounds for addtional variables t, all pre-set as 1 ## */ 
    for(int i = 0; i < _var_t_num; i ++ ){
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FX, make_pair( 1.0, 1.0 ) ); 
        var_bdk.push_back(vb_x);
    }

if(ENFORCE_VEL)
{
    for(int i = 0; i < _vel_var_num; i ++ ){
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_RA, make_pair( - maxVel, + maxVel ) ); 
        var_bdk.push_back(vb_x);
    } 
}

if(ENFORCE_ACC)
{   
    for(int i = 0; i < _acc_var_num; i ++ ){
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_RA, make_pair( - maxAcc, + maxAcc ) ); 
        var_bdk.push_back(vb_x);
    } 
}
    //ROS_WARN("all variables number is: %d", _var_num);
    //ROS_WARN("all constraints number is: %d", _con_num);

    MSKint32t  j,i; 
    MSKenv_t   env; 
    MSKtask_t  task; 
    // Create the mosek environment. 
    r = MSK_makeenv( &env, NULL ); 
  
    // Create the optimization task. 
    r = MSK_maketask(env,_con_num, _var_num, &task); 

// Parameters used in the optimizer
//######################################################################
    MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_CONIC );
    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);

    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-5);
    //MSK_putdouparam (task, MSK_DPAR_INTPNT_CO_TOL_MU_RED, 1e-10);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_CO_TOL_DFEAS,  1e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_CO_TOL_PFEAS,  1e-3);
    //MSK_putdouparam (task, MSK_DPAR_INTPNT_CO_TOL_REL_GAP,5e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_CO_TOL_INFEAS, 1e-3);

//######################################################################
    
    //ROS_WARN("Append all bound values and keys");
    r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
    // Append 'con_num' empty constraints. 
     //The constraints will initially have no bounds. 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendcons(task, _con_num);  

    // Append '_var_num' variables. The variables will initially be fixed at zero (x=0). 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendvars(task, _var_num); 

    for(j = 0; j<_var_num && r == MSK_RES_OK; ++j){ 
      // Set the bounds on variable j : //  blx[j] <= x_j <= bux[j] 
        if (r == MSK_RES_OK) 
            r = MSK_putvarbound(task, 
                                j,                            // Index of variable. 
                                var_bdk[j].first,             // Bound key.
                                var_bdk[j].second.first,      // Numerical value of lower bound.
                                var_bdk[j].second.second );   // Numerical value of upper bound.      
    } 
    
    // Set the bounds on constraints. 
    //   for i=1, ...,con_num : blc[i] <= constraint i <= buc[i] 
    for( i = 0; i < _con_num && r == MSK_RES_OK; i++ ) {
        //cout<<"bounding key: "<<con_bdk[i].first<<endl;
        r = MSK_putconbound(task, 
                            i,                            // Index of constraint. 
                            con_bdk[i].first,             // Bound key.
                            con_bdk[i].second.first,      // Numerical value of lower bound.
                            con_bdk[i].second.second );   // Numerical value of upper bound. 
    }

    //ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, equality part");   
    // #1   put the equality constraints in the start position
    // #1.1 positon constraints in the start point
    // #1.2 velocity constraints in the start point
    // #1.3 acceleration constraints in the start point
    // For position, velocity and acceleration, seperately
    
    int row_idx = 0;
    /*   Start position  */
    {
        // position :
        int order = multi_poly_order[0];
        _s1d1CtrlP_num = order + 1;
        for(int i = 0; i < 3; i++){  // loop for x, y, z       
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = i * _s1d1CtrlP_num;
            aval[0] = 1.0 * initScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++){  // loop for x, y, z       
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = i * _s1d1CtrlP_num;
            asub[1] = i * _s1d1CtrlP_num + 1;
            aval[0] = - 1.0 * order;
            aval[1] =   1.0 * order;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);   
            row_idx ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++){  // loop for x, y, z       
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = i * _s1d1CtrlP_num;
            asub[1] = i * _s1d1CtrlP_num + 1;
            asub[2] = i * _s1d1CtrlP_num + 2;
            aval[0] =   1.0 * order * (order - 1) / initScale;
            aval[1] = - 2.0 * order * (order - 1) / initScale;
            aval[2] =   1.0 * order * (order - 1) / initScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    // #2   put the equality constraints in the end position
    // #2.1 positon constraints in the end point
    // #2.2 velocity constraints in the end point
    // #2.3 acceleration constraints in the end point

    /*   End position  */
    {   
        int order = multi_poly_order[_segment_num - 1];
        _s1d1CtrlP_num = order + 1;
        // position :
        for(int i = 0; i < 3; i++){  // loop for x, y, z       
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
            aval[0] = 1.0 * lstScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++){  // loop for x, y, z       
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 1;
            asub[1] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
            aval[0] = - 1.0;
            aval[1] =   1.0;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++){  // loop for x, y, z       
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 2;
            asub[1] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 1;
            asub[2] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
            aval[0] =   1.0 / lstScale;
            aval[1] = - 2.0 / lstScale;
            aval[2] =   1.0 / lstScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    // #3   put the equality coinstraints in each joint positions 
    // #3.1 positon constraints in each joint positions
    // #3.2 velocity constraints in each joint positions
    // #3.3 acceleration constraints in each joint positions

    /*   joint points  */
    //ROS_WARN(" joint position");
    {
        int sub_shift = 0;
        double val0, val1;
        for(int k = 0; k < (_segment_num - 1); k ++ )
        {
            int order = multi_poly_order[k];
            _s1d1CtrlP_num = order + 1;
            _s1CtrlP_num = 3 * _s1d1CtrlP_num;

            int order_next = multi_poly_order[k + 1];
            int _s1d1CtrlP_next_num = order_next + 1;
            // position :
            val0 = _Scale(k);
            val1 = _Scale(k+1);
            for(int i = 0; i < 3; i++)
            {  // loop for x, y, z
                int nzi = 2;
                MSKint32t asub[nzi];
                double aval[nzi];

                // This segment's last control point
                aval[0] = 1.0 * val0;
                asub[0] = sub_shift + (i+1) * _s1d1CtrlP_num - 1;
                // Next segment's first control point
                aval[1] = -1.0 * val1;
                asub[1] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_next_num;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
            // velocity :
            val0 = order;
            val1 = order_next;
            for(int i = 0; i < 3; i++)
            {  // loop for x, y, z
                int nzi = 4;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] = -1.0 * val0;
                aval[1] =  1.0 * val0;
                asub[0] = sub_shift + (i+1) * _s1d1CtrlP_num - 2;    
                asub[1] = sub_shift + (i+1) * _s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[2] =  1.0 * val1;
                aval[3] = -1.0 * val1;
                asub[2] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_next_num;    
                asub[3] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_next_num + 1;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
            // acceleration :
            val0 = order * (order - 1) / _Scale(k);
            val1 = order_next * (order_next - 1) / _Scale(k+1);
            for(int i = 0; i < 3; i++)
            {  // loop for x, y, z
                int nzi = 6;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] =  1.0  * val0;
                aval[1] = -2.0  * val0;
                aval[2] =  1.0  * val0;
                asub[0] = sub_shift + (i+1) * _s1d1CtrlP_num - 3;    
                asub[1] = sub_shift + (i+1) * _s1d1CtrlP_num - 2;   
                asub[2] = sub_shift + (i+1) * _s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[3] =  -1.0  * val1;
                aval[4] =   2.0  * val1;
                aval[5] =  -1.0  * val1;
                asub[3] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_next_num;    
                asub[4] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_next_num + 1;
                asub[5] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_next_num + 2;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }

            sub_shift += _s1CtrlP_num;
        }
    }

    #if 1
    //ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A for additional linear equality constraints");
    {
        /*** Unknowns sequence: x, w, y, t ***/
        // #0  for all the c'x + w = u induced by corridor constraints
        //ROS_WARN("for linear conrridor constraints");
        int sub_idx   = 0;
        int sub_shift = 0;
        for(int k = 0; k < _segment_num ; k ++ )
        {   
            int order = multi_poly_order[k];
            _s1d1CtrlP_num = order + 1;
            _s1CtrlP_num = 3 * _s1d1CtrlP_num;

            for(int p = 0; p < _s1d1CtrlP_num; p++)
            {
                int nzi = 4;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                for(int i = 0; i < 3; i++)
                {   
                    // for x, y, z, no loop but in a row : x, y, z coupled
                    aval[i] = -2.0 * _Path( k, i ) * _Scale(k);
                    asub[i] = sub_shift + i * _s1d1CtrlP_num + p;    
                }

                aval[3] = 1.0;
                asub[3] = _CtrlP_num + sub_idx;    

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
                sub_idx ++;
            }
            sub_shift += _s1CtrlP_num;
        }

        // #1  for all the Fx = y mapping relationships induced by the corridor constraints
        //ROS_WARN("for variable mapping by corridor");
        sub_idx   = 0;
        sub_shift = 0;
        for(int k = 0; k < _segment_num ; k ++ )
        {   
            int order = multi_poly_order[k];
            _s1d1CtrlP_num = order + 1;
            _s1CtrlP_num = 3 * _s1d1CtrlP_num;
            for(int p = 0; p < _s1CtrlP_num; p++)
            {
                int nzi = 2;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                aval[0] = sqrt(2.0) *_Scale(k); //
                asub[0] = sub_shift + p;    

                aval[1] = -1.0;
                asub[1] = _CtrlP_num + _var_w_num + sub_idx;    

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
                sub_idx ++;
            }
            sub_shift += _s1CtrlP_num;
        }

        // #2  for all the Fx = y mapping relationships induced by minimum snap objective
        //ROS_WARN("for variable mapping by objective");
        sub_idx   = 0;
        sub_shift = 0;
        for(int k = 0; k < _segment_num; k ++)
        {   
            int order = multi_poly_order[k];
            _s1d1CtrlP_num = order + 1;
            _s1CtrlP_num = 3 * _s1d1CtrlP_num;
            for(int p = 0; p < 3; p ++ )
            {
                for( int i = 0; i < _s1d1CtrlP_num - minimize_order; i ++ )
                {   
                    int nzi = _s1d1CtrlP_num - i + 1;
                    MSKint32t asub[nzi];
                    double aval[nzi];

                    for(int j = 0; j < nzi - 1; j ++)
                    {
                        aval[j] = sqrt(2.0) * FMList[order](i, j) / pow(_Scale(k), (2 * minimize_order - 3) / 2.0); 
                        asub[j] = sub_shift + p * _s1d1CtrlP_num + j;    
                    }
                    
                    aval[nzi-1] = -1.0;
                    asub[nzi-1] = _CtrlP_num + _var_w_num + _var_y_con + sub_idx; 

                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                    sub_idx ++;
                }
            }  
            sub_shift += _s1CtrlP_num;
        }
    }
    #endif

    //ROS_WARN("Enforcing high-order constraints in velocity");
if(ENFORCE_VEL)
{
    // The velocity constraints
    int sub_shift   = 0;
    int sub_v_shift = 0;
    for(int k = 0; k < _segment_num ; k ++ )
    {   
        int order = multi_poly_order[k];
        _s1d1CtrlP_num = order + 1;
        _s1CtrlP_num   = 3 * _s1d1CtrlP_num;

        for(int i = 0; i < 3; i++)
        {  // for x, y, z loop
            for(int p = 0; p < order; p++)
            {
                int nzi = 3;
                MSKint32t asub[nzi];
                double aval[nzi];

                aval[0] = -1.0 * order;
                aval[1] =  1.0 * order;
                aval[2] = -1.0;

                asub[0] = sub_shift + i * _s1d1CtrlP_num + p;    
                asub[1] = sub_shift + i * _s1d1CtrlP_num + p + 1;    
                asub[2] = _CtrlP_num  + _var_extra_num + sub_v_shift + i * order + p;    

                //cout<<"segment num : "<<k << endl;
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
        }

        sub_shift   += _s1CtrlP_num;
        sub_v_shift += 3 * order;
    }
}

if(ENFORCE_ACC)
{
    // The velocity constraints
    int sub_shift   = 0;
    int sub_a_shift = 0;
    for(int k = 0; k < _segment_num ; k ++ )
    {   
        int order = multi_poly_order[k];
        _s1d1CtrlP_num = order + 1;
        _s1CtrlP_num   = 3 * _s1d1CtrlP_num;

        for(int i = 0; i < 3; i++)
        {  // for x, y, z loop
            for(int p = 0; p < order - 1; p++)
            {
                int nzi = 4;
                MSKint32t asub[nzi];
                double aval[nzi];

                aval[0] =  1.0 * order * (order - 1) / _Scale(k);
                aval[1] = -2.0 * order * (order - 1) / _Scale(k);
                aval[2] =  1.0 * order * (order - 1) / _Scale(k);
                aval[3] = -1.0;

                asub[0] = sub_shift + i * _s1d1CtrlP_num + p;    
                asub[1] = sub_shift + i * _s1d1CtrlP_num + p + 1;    
                asub[2] = sub_shift + i * _s1d1CtrlP_num + p + 2;    
                asub[3] = _CtrlP_num  + _var_extra_num + _vel_var_num + sub_a_shift + i * (order - 1) + p;    

                //cout<<"segment num : "<<k << endl;
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
        }

        sub_shift   += _s1CtrlP_num;
        sub_a_shift += 3 * (order - 1);
    }
}

    //ROS_WARN("check final row idx after stacking all equalities: %d", row_idx);
    /*** Unknowns sequence: x, w, y, t ***/
    /*
    w num = 1 + _inequ_con_num; 
    y num = (_s1d1CtrlP_num - 3) * 3 * _segment_num + _inequ_con_num * 3;
    t num = 1 + _inequ_con_num; 
    */
    //ROS_WARN("[Bezier Trajectory] Start stacking all conic cones");
    {   
        /** cones induced by quadratical constraints **/
        int sub_idx   = 0;
        int sub_shift = 0;
        for(int k = 0; k < _segment_num ; k ++ )
        {   
            int order = multi_poly_order[k];
            _s1d1CtrlP_num = order + 1;
            _s1CtrlP_num = 3 * _s1d1CtrlP_num;
            for(int p = 0; p < _s1d1CtrlP_num; p++)
            {
                int nzi = 5;
                MSKint32t csub[nzi];                
                
                // Append one rotated conic cone.
                csub[0] = _CtrlP_num + sub_idx;
                csub[1] = _CtrlP_num + _var_w_num + _var_y_num + sub_idx;

                for( int i = 0; i < 3; i++ )
                {   
                    //cout<<"idx for a y variable is : "<<_CtrlP_num + _inequ_con_num + k * _s1CtrlP_num + i * _s1d1CtrlP_num + p<<endl;
                    csub[i + 2] = _CtrlP_num + _var_w_num + sub_shift + i * _s1d1CtrlP_num + p;
                }

                r = MSK_appendcone(task, MSK_CT_RQUAD, 0.0, nzi, csub);
                sub_idx ++;
            }
            sub_shift += _s1CtrlP_num;
        }
        /** the cone induced by the quadratical objective **/
        int nzi = 2 + _obj_nzero_num;
        MSKint32t csub[nzi];                
        
        // Append one rotated conic cone.
        csub[0] = _CtrlP_num + sub_idx;
        csub[1] = _CtrlP_num + _var_w_num + _var_y_num + sub_idx;

        for( int i = 0; i < _obj_nzero_num; i++ )
            csub[i + 2] = _CtrlP_num + _var_w_num + _var_y_con + i;

        r = MSK_appendcone(task, MSK_CT_RQUAD, 0.0, nzi, csub);
    }

    ROS_WARN("[Bezier Trajectory] Start stacking the objective");
    /*** Now no quadratical objective exists, replaced it by a linear variable w ***/    
    /* Set the linear term c_j in the objective.*/  
    MSKint32t csub = _CtrlP_num + _var_w_num - 1; 
    r = MSK_putcj(task, csub, 1.0);

    ros::Time time_end1 = ros::Time::now();
    ROS_WARN("Time in variables stack is");
    cout<<time_end1 - time_start<<endl;

    /*if ( r== MSK_RES_OK )
         r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval); */
    
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
              //solve_ok = true; // debug
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
    
    MSK_deletetask(&task); 
    MSK_deleteenv(&env); 

    ros::Time time_end2 = ros::Time::now();
    ROS_WARN("time consume in optimize is :");
    cout<<time_end2 - time_end1<<endl;

    if(!solve_ok){
      MatrixXd poly_fail = MatrixXd::Identity(3,3);
      ROS_WARN("In solver, falied ");
      return poly_fail;
    }

    VectorXd d_var(_var_num);
    for(int i = 0; i < _var_num; i++)
      d_var(i) = x_var[i];

    int max_order = *max_element( begin( multi_poly_order ), end( multi_poly_order ) );    
    PolyCoeff = MatrixXd::Zero(_segment_num, 3 *(max_order + 1) );

    int var_shift = 0;
    for(int i = 0; i < _segment_num; i++ )
    {
        int order = multi_poly_order[i];
        int poly_num1d = order + 1;
        
        for(int j = 0; j < 3 * poly_num1d; j++)
            PolyCoeff(i , j) = d_var(j + var_shift);

        var_shift += 3 * poly_num1d;
    }   
    
    return PolyCoeff;
}