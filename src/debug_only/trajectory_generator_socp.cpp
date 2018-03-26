#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include "pcd_trajectory/trajectory_generator_socp.h"
#include "pcd_trajectory/mosek.h"
#include "pcd_trajectory/bezier_base.h"

// **************************** FIX ME: Finish before 2017.04.05
// ****************************        1 . Add extreme check and limit adding functions in the complete solver, delete useless sparse matirx functions
// ****************************        2 . Add high order constraints and extreme check and constraints functions
// ****************************        4 . Solver re-wrapper: standalone interface for calling the convex solver. Added loop for extreme iteratively adding

using namespace std;    
using namespace Eigen;

#define inf 1>>30

typedef Eigen::SparseMatrix<double> SMatrixXd;

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
} /* printstr */


TrajectoryGenerator::TrajectoryGenerator(){}

TrajectoryGenerator::~TrajectoryGenerator(){}

Eigen::MatrixXd TrajectoryGenerator::BezierPloyCoeffGenerationSOCP(
            const Eigen::MatrixXd &Path,
            const Eigen::VectorXd &Radius,
            const Eigen::VectorXd &Time,
            const Eigen::MatrixXd &FM,
            const Eigen::MatrixXd &pos,
            const Eigen::MatrixXd &vel,
            const Eigen::MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,       // define the order of the polynomial functions
            const int minimize_order )  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   

#define POLYORDER traj_order
#define MINORDER  minimize_order

    //ROS_WARN("[Bezier Trajectory] In TrajectoryGenerator");
    ros::Time time_start = ros::Time::now();
    Eigen::MatrixXd PolyCoeff;

    _Time   = Time;
    _Radius = Radius;      
    _Scale = _Time;

/*    cout<<"maxVel: "<<maxVel<<endl;
    cout<<"Time: \n"<<_Time<<endl;
    cout<<"Radius: \n"<<_Radius<<endl;*/
    
    //cout<<"Scale: \n"<<_Scale<<endl;
    initScale = _Scale(0);
    lstScale  = _Scale(_Scale.size() - 1);

    //cout<<"init scale: "<<initScale<<" , "<<"last scale: "<<lstScale<<endl;

    assert(_Time.size() == _Radius.size() );

    _Path = Path.block(1, 0, Path.rows() - 2, 3 );
    
    int _segment_num   = _Time.size();
    ROS_WARN("_segment_num : %d",_segment_num);

    Eigen::Vector3d StartPt  = pos.row(0); 
    Eigen::Vector3d EndPt    = pos.row(1); 

    Eigen::Vector3d StartVel = vel.row(0);
    Eigen::Vector3d EndVel   = vel.row(1);

    Eigen::Vector3d StartAcc = acc.row(0);
    Eigen::Vector3d EndAcc   = acc.row(1);

//#####################################################################################################################
//Prepare for constaints and objective data stucture for Mosek solver .         

    int _s1d1CtrlP_num    = POLYORDER + 1;                // The number of control points in ONE single segment and in ONE single Dimension
    int _s1CtrlP_num      = 3 * _s1d1CtrlP_num;           // The number of control points in ONE single segment, 3Dimension (X + Y + Z)
    int _CtrlP_num        = _segment_num * _s1CtrlP_num ; // The number of all control points to be calculated

    int equ_con_s_num = 3 * 3; // p, v, a in x, y, z axis at the start point
    int equ_con_e_num = 3 * 3; // p, v, a in x, y, z axis at the end point
    int equ_con_continuity_num = 3 * 3 * (_segment_num - 1);

    int _equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position
    int _inequ_con_num = _s1d1CtrlP_num * _segment_num; 

    int _obj_nzero_num = (_s1d1CtrlP_num - 3) * 3 * _segment_num;

    // additional equality constraints introduced by the reformulation of all quadratic constraints
    // each quadratic ==> 1 for u - c'x = w, 3 for y = Fx
    // the objective  ==> _obj_nzero_num for y = Fx
    //int _equ_con_extra_num = _obj_nzero_num; 

    int _equ_con_extra_num = _inequ_con_num + _inequ_con_num * 3 + _obj_nzero_num; 
    cout<<"_equ_con_extra_num: "<<_equ_con_extra_num<<endl;

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

    int _con_num = _equ_con_num + _equ_con_extra_num;
    int _var_num = _CtrlP_num + _var_extra_num;

    /*cout<<"_CtrlP_num: "<<_CtrlP_num<<endl;
    cout<<"_var_extra_num: "<<_var_extra_num<<endl;
    cout<<"_var_w_num: "<<_var_w_num<<endl;
    cout<<"_var_y_num: "<<_var_y_num<<endl;
    cout<<"_var_t_num: "<<_var_t_num<<endl;*/

    double x_var[_var_num];
    
    MSKrescodee  r; 
    /* ## define a container for constaints boundary and boundkey ## */ 
    /* ## dataType in the double pair is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    ROS_WARN("[Bezier Trajectory] Start stacking the bounding value for constraints");
    /*** Here is most important stuff, we will elliminate all quadratical constraints and this would introduce many equality constraints and more additional variables***/
    /*************************************************************************************************/
    /** Now firstly we omit the reformulation of the quadratic part in the objective **/
    ROS_WARN(" ... bounding value for original equality induced value");
    for(int i = 0; i < _equ_con_num; i ++ ){ 
        double beq_i;
        if(i < 3)                    beq_i = StartPt(i); 
        else if (i >= 3  && i < 6  ) beq_i = StartVel(i - 3); 
        else if (i >= 6  && i < 9  ) beq_i = StartAcc(i - 6);
        else if (i >= 9  && i < 12 ) beq_i = EndPt (i - 9 );
        else if (i >= 12 && i < 15 ) beq_i = EndVel(i - 12);
        else if (i >= 15 && i < 18 ) beq_i = EndAcc(i - 15);
        else beq_i = 0.0;

        //cout<<"beq_i: "<<beq_i<<endl;
        pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
        con_bdk.push_back(cb_eq);
    }

    #if 1
    ROS_WARN(" ... bounding value for corridor induced value");
    /***  Stack the bounding value for equality constraints induced by the corridor constraints  ***/
    for(int k =0; k < _segment_num; k++ ){
        for(int i = 0; i < _s1d1CtrlP_num; i++){
            double bin_i = _Radius(k) * _Radius(k) - _Path(k, 0) * _Path(k, 0) - _Path(k, 1) * _Path(k, 1) - _Path(k, 2) * _Path(k, 2) ;
            
            //cout<<"bin_i: "<<bin_i<<endl;
            pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_FX, make_pair( bin_i, bin_i ) ); // # cb_ie means: constriants boundary of inequality constrain      
            con_bdk.push_back(cb_ie);   
        }
    }

    ROS_WARN(" ... bounding value for mapping x to y induced value");
    /***  Stack the bounding value for mapping Fx to additional variables y  ***/
    for(int i = 0; i < _var_y_num; i++)
    {
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_FX, make_pair( 0.0, 0.0 ) ); // # cb_ie means: constriants boundary of inequality constrain      
        con_bdk.push_back(cb_ie);   
    }

    #endif

    /*ROS_WARN("con_bdk size is: %d", (int)con_bdk.size());
    cout<<"check constraints boundary"<<endl;
    for(auto ptr:con_bdk)
        cout<<ptr.second.first<<" , "<<ptr.second.second<<endl;*/

    /*** ## Stacking bounds for all unknowns ## ***/ 
    /*** The sequence is control points + additional variables: x, w, y, t ***/
    ROS_WARN("[Bezier Trajectory] Start stacking the bounding value for variables");
    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 
    for(int i = 0; i < _CtrlP_num; i ++ ){
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); 
        var_bdk.push_back(vb_x);
    } 
    
    #if 1
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
    #endif

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
    //MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_INTPNT );
    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);

    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-5);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-2);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-2);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-2 );

//######################################################################
    
    ROS_WARN("Append all bound values and keys");
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

    ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, equality part");   
    // #1   put the equality constraints in the start position
    // #1.1 positon constraints in the start point
    // #1.2 velocity constraints in the start point
    // #1.3 acceleration constraints in the start point
    // For position, velocity and acceleration, seperately
    
    int row_idx = 0;
    /*   Start position  */
    // position :
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
        aval[0] = - 1.0 * POLYORDER;
        aval[1] =   1.0 * POLYORDER;
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
        aval[0] =   1.0 * POLYORDER * (POLYORDER - 1) / initScale;
        aval[1] = - 2.0 * POLYORDER * (POLYORDER - 1) / initScale;
        aval[2] =   1.0 * POLYORDER * (POLYORDER - 1) / initScale;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
        row_idx ++;
    }

    // #2   put the equality constraints in the end position
    // #2.1 positon constraints in the end point
    // #2.2 velocity constraints in the end point
    // #2.3 acceleration constraints in the end point

    /*   End position  */
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

    // #3   put the equality coinstraints in each joint positions 
    // #3.1 positon constraints in each joint positions
    // #3.2 velocity constraints in each joint positions
    // #3.3 acceleration constraints in each joint positions

    /*   joint points  */
    double val0, val1;
    for(int k = 0; k < (_segment_num - 1); k ++ ){
        // position :
        val0 = _Scale(k);
        val1 = _Scale(k+1);
        //cout<<"val0 : "<<val0<<" , val1: "<<val1<<endl;
        for(int i = 0; i < 3; i++){  // loop for x, y, z
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];

            // This segment's last control point
            aval[0] = 1.0 * val0;
            asub[0] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 1;
            // Next segment's first control point
            aval[1] = -1.0 * val1;
            asub[1] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num;

            //cout<<"sub 0: "<<asub[0]<<" , sub 1: "<<asub[1]<<endl;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // velocity :
        val0 = val1 = 1.0;
        for(int i = 0; i < 3; i++){  // loop for x, y, z
            int nzi = 4;
            MSKint32t asub[nzi];
            double aval[nzi];
            
            // This segment's last velocity control point
            aval[0] = -1.0 * val0;
            aval[1] =  1.0 * val0;
            asub[0] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 2;    
            asub[1] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 1;   
            // Next segment's first velocity control point
            aval[2] =  1.0 * val1;
            aval[3] = -1.0 * val1;
            asub[2] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num;    
            asub[3] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num + 1;

            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // acceleration :
        val0 = 1.0/ _Scale(k);
        val1 = 1.0/ _Scale(k+1);
        for(int i = 0; i < 3; i++){  // loop for x, y, z
            int nzi = 6;
            MSKint32t asub[nzi];
            double aval[nzi];
            
            // This segment's last velocity control point
            aval[0] =  1.0  * val0;
            aval[1] = -2.0  * val0;
            aval[2] =  1.0  * val0;
            asub[0] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 3;    
            asub[1] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 2;   
            asub[2] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 1;   
            // Next segment's first velocity control point
            aval[3] =  -1.0  * val1;
            aval[4] =   2.0  * val1;
            aval[5] =  -1.0  * val1;
            asub[3] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num;    
            asub[4] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num + 1;
            asub[5] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num + 2;

            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    #if 1
    ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A for additional linear equality constraints");
    {
        /*** Unknowns sequence: x, w, y, t ***/
        // #0  for all the c'x + w = u induced by corridor constraints
        ROS_WARN("for linear conrridor constraints");
        int sub_idx = 0;
        for(int k = 0; k < _segment_num ; k ++ )
        {
            for(int p = 0; p < _s1d1CtrlP_num; p++)
            {
                int nzi = 4;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                for(int i = 0; i < 3; i++)
                {   
                    // for x, y, z, no loop but in a row : x, y, z coupled
                    aval[i] = -2.0 * _Path( k, i ) * _Scale(k);
                    asub[i] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p;    
                }

                aval[3] = 1.0;
                asub[3] = _CtrlP_num + sub_idx;    

                /*cout<<"No. "<<row_idx<<" constraints"<<endl;
                for(int i = 0; i < nzi; i++)
                {
                    //cout<<"val: "<<aval[i]<<endl;
                    cout<<"sub: "<<asub[i]<<" , ";
                }
                cout<<"\n"<<endl;*/

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
                sub_idx ++;
            }
        }

        // #1  for all the Fx = y mapping relationships induced by the corridor constraints
        ROS_WARN("for variable mapping by corridor");
        sub_idx = 0;
        for(int k = 0; k < _segment_num ; k ++ )
        {
            for(int p = 0; p < _s1CtrlP_num; p++)
            {
                int nzi = 2;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                aval[0] = sqrt(2.0) *_Scale(k); //
                asub[0] = k * _s1CtrlP_num + p;    

                aval[1] = -1.0;
                asub[1] = _CtrlP_num + _var_w_num + sub_idx;    

/*                cout<<"No. "<<row_idx<<" constraints"<<endl;
                for(int i = 0; i < nzi; i++)
                {
                    //cout<<"val: "<<aval[i]<<endl;
                    cout<<"sub: "<<asub[i]<<" , ";
                }
                cout<<"\n"<<endl;*/

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
                sub_idx ++;
            }
        }

        // #2  for all the Fx = y mapping relationships induced by minimum snap objective
        ROS_WARN("for variable mapping by objective");
        sub_idx = 0;
        for(int k = 0; k < _segment_num; k ++)
        {
            for(int p = 0; p < 3; p ++ )
            {
                for( int i = 0; i < _s1d1CtrlP_num - 3; i ++ )
                {   
                    int nzi = _s1d1CtrlP_num - i + 1;
                    MSKint32t asub[nzi];
                    double aval[nzi];

  /*                  cout<<"No :"<<row_idx<<" equality "<<endl;
                    cout<<"nzi: "<<nzi<<endl;
*/
                    for(int j = 0; j < nzi - 1; j ++)
                    {
                        aval[j] = FM(i, j) / pow(_Scale(k), 1.5); 
                        asub[j] = k * _s1CtrlP_num + p * _s1d1CtrlP_num + j;    
                        
                        //cout<<"val: "<<aval[j]<<endl;
                        //cout<<"sub: "<<asub[j]<<endl;
                    }
                    
                    //cout<<"j+1"<<j+1<<endl;
                    aval[nzi-1] = -1.0;
                    asub[nzi-1] = _CtrlP_num + _var_w_num + _var_y_con + sub_idx; 
                    
                    /*cout<<"No. "<<row_idx<<" constraints"<<endl;
                    for(int i = 0; i < nzi; i++)
                    {
                        //cout<<"val: "<<aval[i]<<endl;
                        cout<<"sub: "<<asub[i]<<" , ";
                    }
                    cout<<"\n"<<endl;*/

                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                    sub_idx ++;
                }
            }  
        }

    /*    int NUMQ_blk = (POLYORDER + 1); // default minimize the jerk and MINORDER = 3
    int NUMQNZ   = _segment_num * 3 * NUMQ_blk * (NUMQ_blk + 1) / 2;
    MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];
    double     qval[NUMQNZ];

    int idx = 0;
    for(int k = 0; k < _segment_num; k ++){
        for(int p = 0; p < 3; p ++ ){
            for( int i = 0; i < _s1d1CtrlP_num; i ++ ){
                for( int j = 0; j < _s1d1CtrlP_num; j ++ ){
                    if( i >= j ){
                        qsubi[idx] = k * _s1CtrlP_num + p * _s1d1CtrlP_num + i;   
                        qsubj[idx] = k * _s1CtrlP_num + p * _s1d1CtrlP_num + j;  
                        qval[idx]  = MQM(i, j) / pow(_Scale(k), 3);                
                        idx ++ ;
                    }
                }
            }
        }  
    } */

    }
    #endif

    ROS_WARN("check final row idx after stacking all equalities: %d", row_idx);
    /*** Unknowns sequence: x, w, y, t ***/
    /*
    w num = 1 + _inequ_con_num; 
    y num = (_s1d1CtrlP_num - 3) * 3 * _segment_num + _inequ_con_num * 3;
    t num = 1 + _inequ_con_num; 
    */

    ROS_WARN("[Bezier Trajectory] Start stacking all conic cones");
    {   
        /** cones induced by quadratical constraints **/
        int sub_idx = 0;
        for(int k = 0; k < _segment_num ; k ++ )
        {
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
                    csub[i + 2] = _CtrlP_num + _var_w_num + k * _s1CtrlP_num + i * _s1d1CtrlP_num + p;
                }

                /*cout<<"No. "<<sub_idx<<" cone "<<endl;
                for(int i = 0; i < nzi; i++)
                {
                    //cout<<"val: "<<aval[i]<<endl;
                    cout<<"sub: "<<csub[i]<<" , ";
                }
                cout<<"\n"<<endl;*/

                r = MSK_appendcone(task, MSK_CT_RQUAD, 0.0, nzi, csub);
                sub_idx ++;
            }
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

    PolyCoeff.resize(_segment_num, _s1CtrlP_num);

    if(!solve_ok){
      MatrixXd poly_fail = MatrixXd::Identity(3,3);
      ROS_WARN("In solver, falied ");
      return poly_fail;
    }

    VectorXd d_var(_var_num);
    for(int i = 0; i < _var_num; i++)
      d_var(i) = x_var[i];

    VectorXd p_var = d_var;
    //cout<<"solution: "<<p_var<<endl;

    for(int i = 0; i < _CtrlP_num; i++)
        PolyCoeff(i / _s1CtrlP_num, i % _s1CtrlP_num) = p_var(i);
    
    return PolyCoeff;
}

Eigen::MatrixXd TrajectoryGenerator::BezierPloyCoeffGenerationQCQP(
            const Eigen::MatrixXd &Path,
            const Eigen::VectorXd &Radius,
            const Eigen::VectorXd &Time,
            const Eigen::MatrixXd &MQM,
            const Eigen::MatrixXd &pos,
            const Eigen::MatrixXd &vel,
            const Eigen::MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,       // define the order of the polynomial functions
            const int minimize_order )  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   

#define POLYORDER traj_order
#define MINORDER  minimize_order

    //ROS_WARN("[Bezier Trajectory] In TrajectoryGenerator");
    ros::Time time_start = ros::Time::now();
    Eigen::MatrixXd PolyCoeff;

    _Time   = Time;
    _Radius = Radius;      
    _Scale = _Time;

/*    for(int i = 0; i < _Time.size(); i++)
        cout<<"ratio "<<i<<" : "<<2 * _Radius(i) / _Time(i)<<endl;*/
    
    //cout<<"Scale: \n"<<_Scale<<endl;
    initScale = _Scale(0);
    lstScale  = _Scale(_Scale.size() - 1);

    //cout<<"init scale: "<<initScale<<" , "<<"last scale: "<<lstScale<<endl;

    assert(_Time.size() == _Radius.size() );

    _Path = Path.block(1, 0, Path.rows() - 2, 3 );
    
    int _segment_num   = _Time.size();


    Eigen::Vector3d StartPt  = pos.row(0); 
    Eigen::Vector3d EndPt    = pos.row(1); 

    Eigen::Vector3d StartVel = vel.row(0);
    Eigen::Vector3d EndVel   = vel.row(1);

    Eigen::Vector3d StartAcc = acc.row(0);
    Eigen::Vector3d EndAcc   = acc.row(1);

//#####################################################################################################################
//Prepare for constaints and objective data stucture for Mosek solver .         

    int _s1d1CtrlP_num    = POLYORDER + 1;                // The number of control points in ONE single segment and in ONE single Dimension
    int _s1CtrlP_num      = 3 * _s1d1CtrlP_num;           // The number of control points in ONE single segment, 3Dimension (X + Y + Z)
    int _CtrlP_num        = _segment_num * _s1CtrlP_num ; // The number of all control points to be calculated

    int equ_con_s_num = 3 * 3; // p, v, a in x, y, z axis at the start point
    int equ_con_e_num = 3 * 3; // p, v, a in x, y, z axis at the end point
    int equ_con_continuity_num = 3 * 3 * (_segment_num - 1);
    int _equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position

    int _inequ_con_num = _s1d1CtrlP_num * _segment_num ; 
    //int _high_order_con_num =  (_s1d1CtrlP_num - 1 + _s1d1CtrlP_num - 2) * 3 * _segment_num; //Now we do bound the velocity and the acceleration
    int _high_order_con_num = 0;// POLYORDER * 3 * _segment_num; //Now we only bound the velocity
    //int _high_order_con_num   =  0;  // Now we do not bound the velocity and the acceleration
    int con_num = _equ_con_num + _inequ_con_num + _high_order_con_num;

    double x_var[ _CtrlP_num];
    
    MSKrescodee  r; 
    /* ## define a container for constaints boundary and boundkey ## */ 
    /* ## dataType in the double pair is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    
    //ROS_WARN("[Bezier Trajectory] Start stacking the bounding value");
    /***  Stack the bounding value for the quadratic inequality for the corridor constraints  ***/
    for(int k =0; k < _segment_num; k++ ){
        for(int i = 0; i < _s1d1CtrlP_num; i++){
            double bin_i = _Radius(k) * _Radius(k) - _Path(k, 0) * _Path(k, 0) - _Path(k, 1) * _Path(k, 1) - _Path(k, 2) * _Path(k, 2) ;
            
            pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_UP, make_pair( - MSK_INFINITY, bin_i ) ); // # cb_ie means: constriants boundary of inequality constrain      
            con_bdk.push_back(cb_ie);   
        }
    }

    /***  Stack the bounding value for the linear inequality for the high order dynamic constraints  ***/
    // for velocity constraints
/*    for(int k =0; k < _segment_num; k++ ){
        for(int p = 0; p < 3; p++) // loop for x, y, z
            for(int i = 0; i < POLYORDER; i++){
                pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxVel,  maxVel) );
                con_bdk.push_back(cb_ie);   
            }
    }*/

    // for acceleration constraints
/*    for(int k =0; k < _segment_num; k++ ){
        for(int p = 0; p < 3; p++) // loop for x, y, z
            for(int i = 0; i < _s1d1CtrlP_num - 2; i++){
                pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxAcc,  maxAcc) ); 
                con_bdk.push_back(cb_ie);   
            }
    }*/

 //   ROS_WARN("[Bezier Trajectory] equality bound %d", _equ_con_num);
    for(int i = 0; i < _equ_con_num; i ++ ){ 
        double beq_i;
        if(i < 3)                    beq_i = StartPt(i); 
        else if (i >= 3  && i < 6  ) beq_i = StartVel(i - 3); 
        else if (i >= 6  && i < 9  ) beq_i = StartAcc(i - 6);
        else if (i >= 9  && i < 12 ) beq_i = EndPt (i - 9 );
        else if (i >= 12 && i < 15 ) beq_i = EndVel(i - 12);
        else if (i >= 15 && i < 18 ) beq_i = EndAcc(i - 15);
        else beq_i = 0.0;

        //cout<<"beq_i: "<<beq_i<<endl;
        pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
        con_bdk.push_back(cb_eq);
    }

    /* ## define a container for control points' boundary and boundkey ## */ 
    /* ## dataType in one tuple is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 
    for(int i = 0; i < _CtrlP_num; i ++ ){
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)
        var_bdk.push_back(vb_x);
    } 
  
    MSKint32t  j,i; 
    MSKenv_t   env; 
    MSKtask_t  task; 
    // Create the mosek environment. 
    r = MSK_makeenv( &env, NULL ); 
  
    // Create the optimization task. 
    r = MSK_maketask(env,con_num, _CtrlP_num, &task); 

// Parameters used in the optimizer
//######################################################################
    MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_INTPNT );
    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);
    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-5);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-2);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-2);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-2 );

//######################################################################
    
    r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
    // Append 'con_num' empty constraints. 
     //The constraints will initially have no bounds. 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendcons(task,con_num);  

    // Append '_CtrlP_num' variables. The variables will initially be fixed at zero (x=0). 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendvars(task,_CtrlP_num); 

    for(j = 0; j<_CtrlP_num && r == MSK_RES_OK; ++j){ 
      // Set the bounds on variable j : //  blx[j] <= x_j <= bux[j] 
        if (r == MSK_RES_OK) 
            r = MSK_putvarbound(task, 
                                j,                            // Index of variable. 
                                var_bdk[j].first,             // Bound key.
                                var_bdk[j].second.first,      // Numerical value of lower bound.
                                var_bdk[j].second.second );   // Numerical value of upper bound.      
    } 
    
   /* ROS_WARN("con_num: %d", con_num);
    ROS_WARN("con_bdk size: %d", int(con_bdk.size() ));*/

    // Set the bounds on constraints. 
    //   for i=1, ...,con_num : blc[i] <= constraint i <= buc[i] 
    for( i = 0; i < con_num && r == MSK_RES_OK; i++ ) {
        //cout<<"bounding key: "<<con_bdk[i].first<<endl;
        r = MSK_putconbound(task, 
                            i,                            // Index of constraint. 
                            con_bdk[i].first,             // Bound key.
                            con_bdk[i].second.first,      // Numerical value of lower bound.
                            con_bdk[i].second.second );   // Numerical value of upper bound. 
    }

  //  ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, inequality part");
    // Put variables in A by row one by one
    int row_idx_shift = 0;
    int row_idx;
    // #0   put the in-equality coinstraints in each control points
    // #0.1 positon constraints of each control point in the corresponding sphere
    // #0.2 velocity constraints of each 1st order control points within the limit
    // #0.3 acceleration constraints of each 2nd order control points within the limit
    for(int k = 0; k < _segment_num ; k ++ )
    {
        for(int p = 0; p < _s1d1CtrlP_num; p++)
        {
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            
            row_idx = row_idx_shift + _s1d1CtrlP_num * k + p;
            for(int i = 0; i < 3; i++)
            {   // for x, y, z, no loop but in a row : x, y, z coupled
                aval[i] = -2.0 * _Path( k, i ) * _Scale(k);
                asub[i] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p;    
            }
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
        }
    }
    row_idx_shift += _s1d1CtrlP_num * _segment_num;

    // The velocity constraints
/*    row_idx = row_idx_shift;
    for(int k = 0; k < _segment_num ; k ++ ){
        for(int i = 0; i < 3; i++){  // for x, y, z loop
            for(int p = 0; p < POLYORDER; p++){
                
                int nzi = 2;
                MSKint32t asub[nzi];
                double aval[nzi];

                aval[0] = -1.0 * POLYORDER;
                aval[1] =  1.0 * POLYORDER;
                asub[0] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p;    
                asub[1] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p + 1;    
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
        }
    }
    row_idx_shift += 3 * (_s1d1CtrlP_num - 1) * _segment_num;*/

    // The acceleration constraints
/*    row_idx = row_idx_shift;
    for(int k = 0; k < _segment_num ; k ++ ){
        for(int i = 0; i < 3; i++){  // for x, y, z loop
            for(int p = 0; p < _s1d1CtrlP_num - 2; p++){
                
                int nzi = 3;
                MSKint32t asub[nzi];
                double aval[nzi];
                aval[0] =  1.0 * POLYORDER * (POLYORDER - 1) / _Scale(k);
                aval[1] = -2.0 * POLYORDER * (POLYORDER - 1) / _Scale(k);
                aval[2] =  1.0 * POLYORDER * (POLYORDER - 1) / _Scale(k);
                asub[0] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p;    
                asub[1] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p + 1;    
                asub[2] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p + 2;    
                
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
        }
    }
    row_idx_shift += 3 * (_s1d1CtrlP_num - 2) * _segment_num;*/

    // ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, equality part");   
    // #1   put the equality constraints in the start position
    // #1.1 positon constraints in the start point
    // #1.2 velocity constraints in the start point
    // #1.3 acceleration constraints in the start point
    // For position, velocity and acceleration, seperately
     
    /*   Start position  */
    // position :
    for(int i = 0; i < 3; i++){  // loop for x, y, z       
        int nzi = 1;
        MSKint32t asub[nzi];
        double aval[nzi];
        row_idx = row_idx_shift + i;
        asub[0] = i * _s1d1CtrlP_num;
        aval[0] = 1.0 * initScale;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
    }
    // velocity :
    for(int i = 0; i < 3; i++){  // loop for x, y, z       
        int nzi = 2;
        MSKint32t asub[nzi];
        double aval[nzi];
        row_idx = row_idx_shift + 3 + i;
        asub[0] = i * _s1d1CtrlP_num;
        asub[1] = i * _s1d1CtrlP_num + 1;
        aval[0] = - 1.0 * POLYORDER;
        aval[1] =   1.0 * POLYORDER;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);   
    }
    // acceleration : 
    for(int i = 0; i < 3; i++){  // loop for x, y, z       
        int nzi = 3;
        MSKint32t asub[nzi];
        double aval[nzi];
        row_idx = row_idx_shift + 6 + i;
        asub[0] = i * _s1d1CtrlP_num;
        asub[1] = i * _s1d1CtrlP_num + 1;
        asub[2] = i * _s1d1CtrlP_num + 2;
        aval[0] =   1.0 * POLYORDER * (POLYORDER - 1) / initScale;
        aval[1] = - 2.0 * POLYORDER * (POLYORDER - 1) / initScale;
        aval[2] =   1.0 * POLYORDER * (POLYORDER - 1) / initScale;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
    }
    row_idx_shift += 9;

    // #2   put the equality constraints in the end position
    // #2.1 positon constraints in the end point
    // #2.2 velocity constraints in the end point
    // #2.3 acceleration constraints in the end point

    /*   End position  */
    // position :
    for(int i = 0; i < 3; i++){  // loop for x, y, z       
        int nzi = 1;
        MSKint32t asub[nzi];
        double aval[nzi];
        row_idx = row_idx_shift + i;
        asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
        aval[0] = 1.0 * lstScale;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
    }
    // velocity :
    for(int i = 0; i < 3; i++){  // loop for x, y, z       
        int nzi = 2;
        MSKint32t asub[nzi];
        double aval[nzi];
        row_idx = row_idx_shift + 3 + i;
        asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 1;
        asub[1] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
        aval[0] = - 1.0;
        aval[1] =   1.0;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
    }
    // acceleration : 
    for(int i = 0; i < 3; i++){  // loop for x, y, z       
        int nzi = 3;
        MSKint32t asub[nzi];
        double aval[nzi];
        row_idx = row_idx_shift + 6 + i;
        asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 2;
        asub[1] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 1;
        asub[2] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
        aval[0] =   1.0 / lstScale;
        aval[1] = - 2.0 / lstScale;
        aval[2] =   1.0 / lstScale;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
    }
    row_idx_shift += 9;

    // #3   put the equality coinstraints in each joint positions 
    // #3.1 positon constraints in each joint positions
    // #3.2 velocity constraints in each joint positions
    // #3.3 acceleration constraints in each joint positions

    /*   joint points  */
    double val0, val1;
    for(int k = 0; k < (_segment_num - 1); k ++ ){
        // position :
        val0 = _Scale(k);
        val1 = _Scale(k+1);
        for(int i = 0; i < 3; i++){  // loop for x, y, z
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];

            row_idx = row_idx_shift + k * 9 + i;
            // This segment's last control point
            aval[0] = 1.0 * val0;
            asub[0] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 1;
            // Next segment's first control point
            aval[1] = -1.0 * val1;
            asub[1] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num;

            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
        }
        // velocity :
        val0 = val1 = 1.0;
        for(int i = 0; i < 3; i++){  // loop for x, y, z
            int nzi = 4;
            MSKint32t asub[nzi];
            double aval[nzi];
            
            row_idx = row_idx_shift + k * 9 + 3 + i;
            // This segment's last velocity control point
            aval[0] = -1.0 * val0;
            aval[1] =  1.0 * val0;
            asub[0] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 2;    
            asub[1] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 1;   
            // Next segment's first velocity control point
            aval[2] =  1.0 * val1;
            aval[3] = -1.0 * val1;
            asub[2] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num;    
            asub[3] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num + 1;

            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
        }
        // acceleration :
        val0 = 1.0/ _Scale(k);
        val1 = 1.0/ _Scale(k+1);
        for(int i = 0; i < 3; i++){  // loop for x, y, z
            int nzi = 6;
            MSKint32t asub[nzi];
            double aval[nzi];
            
            row_idx = row_idx_shift + k * 9 + 6 + i;
            // This segment's last velocity control point
            aval[0] =  1.0  * val0;
            aval[1] = -2.0  * val0;
            aval[2] =  1.0  * val0;
            asub[0] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 3;    
            asub[1] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 2;   
            asub[2] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 1;   
            // Next segment's first velocity control point
            aval[3] =  -1.0  * val1;
            aval[4] =   2.0  * val1;
            aval[5] =  -1.0  * val1;
            asub[3] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num;    
            asub[4] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num + 1;
            asub[5] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num + 2;

            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
        }
    }

    ROS_WARN("[Bezier Trajectory] Start stacking the objective");
    //cout<<"MQM: \n"<<MQM<<endl;
    // ROS_WARN("[Solver] prepare for stacking the objective");
    int NUMQ_blk = (POLYORDER + 1); // default minimize the jerk and MINORDER = 3
    int NUMQNZ   = _segment_num * 3 * NUMQ_blk * (NUMQ_blk + 1) / 2;
    MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];
    double     qval[NUMQNZ];

    int idx = 0;
    for(int k = 0; k < _segment_num; k ++){
        for(int p = 0; p < 3; p ++ ){
            for( int i = 0; i < _s1d1CtrlP_num; i ++ ){
                for( int j = 0; j < _s1d1CtrlP_num; j ++ ){
                    if( i >= j ){
                        qsubi[idx] = k * _s1CtrlP_num + p * _s1d1CtrlP_num + i;   
                        qsubj[idx] = k * _s1CtrlP_num + p * _s1d1CtrlP_num + j;  
                        qval[idx]  = MQM(i, j) / pow(_Scale(k), 2 * minimize_order - 3);                
                        idx ++ ;
                    }
                }
            }
        }  
    } 
    
    ros::Time time_end1 = ros::Time::now();
    ROS_WARN("Time in variables stack is");
    cout<<time_end1 - time_start<<endl;

    if ( r== MSK_RES_OK )
         r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval); 
    
    if ( r==MSK_RES_OK ) 
         r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
    
    //ROS_WARN("[Bezier Trajectory] Start stacking the Quadratic Matrix");
    int quad_idx;
    for( int k = 0; k < _segment_num; k ++ ){
        double val = 2.0 * _Scale(k) * _Scale(k);
        for(int p = 0; p < _s1d1CtrlP_num; p++)
        {
            int nzi = 3;
            MSKint32t   qisubi[nzi], qisubj[nzi];
            double      qival[nzi];

            qisubi[0] = k * _s1CtrlP_num + p;   
            qisubj[0] = k * _s1CtrlP_num + p;  
            qival [0] = val; 

            qisubi[1] = k * _s1CtrlP_num + _s1d1CtrlP_num + p;   
            qisubj[1] = k * _s1CtrlP_num + _s1d1CtrlP_num + p;  
            qival [1] = val;
            
            qisubi[2] = k * _s1CtrlP_num + 2 * _s1d1CtrlP_num + p;   
            qisubj[2] = k * _s1CtrlP_num + 2 * _s1d1CtrlP_num + p;  
            qival [2] = val;

            quad_idx = _s1d1CtrlP_num * k + p;
            r = MSK_putqconk(task, quad_idx, nzi, qisubi, qisubj, qival);  
        }
    }

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

    PolyCoeff.resize(_segment_num, _s1CtrlP_num);

    if(!solve_ok){
      MatrixXd poly_fail = MatrixXd::Identity(3,3);
      ROS_WARN("In solver, falied ");
      return poly_fail;
    }

    VectorXd d_var(_CtrlP_num);
    for(int i = 0; i < _CtrlP_num; i++)
      d_var(i) = x_var[i];

    VectorXd p_var = d_var;
    //cout<<"solution: "<<p_var<<endl;

    for(int i = 0; i < _CtrlP_num; i++)
        PolyCoeff(i / _s1CtrlP_num, i % _s1CtrlP_num) = p_var(i);
    
    return PolyCoeff;
}