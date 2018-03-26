#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include "pcd_trajectory/trajectory_generator_discrete.h"
#include "pcd_trajectory/mosek.h"

using namespace std;    
using namespace Eigen;

#define inf 1>>30

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
} /* printstr */


TrajectoryGeneratorDiscrete::TrajectoryGeneratorDiscrete(){}

TrajectoryGeneratorDiscrete::~TrajectoryGeneratorDiscrete(){}

MatrixXd TrajectoryGeneratorDiscrete::AccDiscreteGeneration(
            const vector<Vector3d> &Path,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double & max_v,
            const double & max_a,
            const double & h,
            const double & l)
{
    ros::Time time_start = ros::Time::now();

    _Path = Path;

    int var_num_1d  = Path.size();
    int var_num = 3 * var_num_1d;

    MatrixXd AccValue(3, var_num_1d);

    Vector3d P0  = Path.front(); 
    Vector3d Pg  = Path.back(); 

/*    cout<<"h: "<<h<<endl;
    cout<<"l: "<<l<<endl;

    ROS_WARN("check start position: ");
    cout<<P0<<endl;

    ROS_WARN("check target position: ");
    cout<<Pg<<endl;*/

    Vector3d V0 = vel.row(0);
    Vector3d Vg = vel.row(1);

    Vector3d A0 = acc.row(0);
    Vector3d Ag = acc.row(1);

    MatrixXd W(var_num_1d, var_num_1d);
    MatrixXd H(var_num_1d, var_num_1d);

    for(int i = 0; i < var_num_1d; i++ )
        for(int j = 0; j < var_num_1d; j++ )
        {
            if( i == j)
                W(i, j) = - 1.0 / h;
            else if ( i == j - 1)
                W(i, j) = 1.0 / h;
            else
                W(i, j) = 0.0;
        }

    H = W.transpose() * W;

    //ROS_WARN("[Discrete Solver] finish stack H matrix");
/*    cout<<"H: \n"<<H<<endl;*/
//#####################################################################################################################
    int _equ_con_num  = 3 * 2;    
    int _inequ_con_num = 2 * var_num - 6;

    int con_num = _inequ_con_num + _equ_con_num;

    double x_var[ var_num];
    MSKrescodee  r; 
    /* ## define a container for constaints boundary and boundkey ## */ 
    /* ## dataType in the double pair is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    pair<MSKboundkeye, pair<double, double> > cb_ie;
    
    // for position Pk constraints
    for(int p = 0; p < 3; p++)
    {   
        double lo_bound, up_bound;
        for(int k = 1; k < var_num_1d - 1; k++)
        {   
            lo_bound = Path[k](p) - l - P0(p) - h * k * V0(p);
            up_bound = Path[k](p) + l - P0(p) - h * k * V0(p);

            cb_ie = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); 
            con_bdk.push_back(cb_ie);   
        }
    }

    // for velocity Vk constraints
    for(int p = 0; p < 3; p++)
    {   
        double lo_bound, up_bound;
        for(int k = 1; k < var_num_1d - 1; k++)
        {   
            lo_bound = - max_v - V0(p);
            up_bound = + max_v - V0(p);

            cb_ie = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); 
            con_bdk.push_back(cb_ie);   
        }
    }

    for(int i = 0; i < _equ_con_num; i ++ )
    { 
        double beq_i;
        pair<MSKboundkeye, pair<double, double> > cb_eq;

        if(i < 3)                    
        {
            beq_i = Pg(i) - P0(i) - h * (var_num_1d - 1) * V0(i); 
            cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
        }
        else 
        {
            beq_i = Vg(i - 3) - V0(i - 3); 
            cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
        }

        con_bdk.push_back(cb_eq);
    }

/*    ROS_WARN("[Discrete Solver] finish stack constraint bounds");
    ROS_WARN("[Discrete Solver] constraints num is %d", (int)con_bdk.size());
    for(auto ptr:con_bdk)
        cout<<ptr.second.first<<" , "<<ptr.second.second<<endl;*/

    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk;
    pair<MSKboundkeye, pair<double, double> > var_ie;

    for(int p = 0; p < 3; p++)
    {   
        double lo_bound, up_bound;
        
        // the initial acceleraion in each axis
        lo_bound = A0(p);
        up_bound = A0(p);
        var_ie = make_pair( MSK_BK_FX, make_pair( lo_bound, up_bound ) ); 
        var_bdk.push_back(var_ie);   

        // the accelerations a0 ~ ak-1
        for(int k = 1; k < var_num_1d - 1; k++)
        {   
            lo_bound = - max_a;
            up_bound = + max_a;

            var_ie = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); 
            var_bdk.push_back(var_ie);   
        }
        
        // the last acceleration ak
        
        /*lo_bound = Ag(p);
        up_bound = Ag(p);*/
        
        lo_bound = - max_a;
        up_bound = + max_a;

        var_ie = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); 
        var_bdk.push_back(var_ie);   
    }
    
/*    ROS_WARN("[Discrete Solver] finish stack variable bounds");
    ROS_WARN("[Discrete Solver] variable num is %d", (int)var_bdk.size());
    for(auto ptr:var_bdk)
        cout<<ptr.second.first<<" , "<<ptr.second.second<<endl;*/

    // variable sequence: ax, ay, az
    MSKint32t  j,i; 
    MSKenv_t   env; 
    MSKtask_t  task; 
    // Create the mosek environment. 
    r = MSK_makeenv( &env, NULL ); 
  
    // Create the optimization task. 
    r = MSK_maketask(env,con_num, var_num, &task); 

// Parameters used in the optimizer
//######################################################################
    MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_INTPNT );
    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);

    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-5);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-3);
    //MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-3 );

//######################################################################
    
    r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
    if ( r == MSK_RES_OK ) 
    {
      r = MSK_appendcons(task,con_num);  
      r = MSK_appendvars(task,var_num); 
    }
 
    for(j = 0; j < var_num && r == MSK_RES_OK; ++j){ 
        if (r == MSK_RES_OK) 
            r = MSK_putvarbound(task, 
                                j,                            // Index of variable. 
                                var_bdk[j].first,             // Bound key.
                                var_bdk[j].second.first,      // Numerical value of lower bound.
                                var_bdk[j].second.second );   // Numerical value of upper bound.      
    } 
    
    for( i = 0; i < con_num && r == MSK_RES_OK; i++ ) {
        //cout<<"bounding key: "<<con_bdk[i].first<<endl;
        r = MSK_putconbound(task, 
                            i,                            // Index of constraint. 
                            con_bdk[i].first,             // Bound key.
                            con_bdk[i].second.first,      // Numerical value of lower bound.
                            con_bdk[i].second.second );   // Numerical value of upper bound. 
    }

    //ROS_WARN("[Discrete Solver] Start stacking the Linear Matrix A, inequality part");
    // in-equality constraints for each positon :
    int row_idx = 0;
    for(int p = 0; p < 3; p++)
    {
        for(int k = 1; k < var_num_1d - 1; k ++)
        {   
            // p_k = p_0 + h * k * v_0 + h^2/2 + sum(2(k-1) + 3)a_i|0, k-1
            int nzi = k;
            MSKint32t asub[nzi];
            double aval[nzi];

            for(int i = 0; i < nzi; i++)
            { 
                aval[i] = h * h * (2 * ( k - i ) + 3) / 2.0;
                asub[i] = p * var_num_1d + i;    
            }

            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    // in-equality constraints for each velocity :
    for(int p = 0; p < 3; p++)
    {
        for(int k = 1; k < var_num_1d - 1; k ++)
        {   
            // p_k = p_0 + h * k * v_0 + h^2/2 + sum(2(k-1) + 3)a_i|0, k-1
            int nzi = k + 1;
            MSKint32t asub[nzi];
            double aval[nzi];

            for(int i = 0; i < nzi; i++)
            { 
                aval[i] = h;
                asub[i] = p * var_num_1d + i;    
            }

            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }
    
    // #2.1   put the equality constraints in the end position
    for(int p = 0; p < 3; p++)
    {
        int nzi = var_num_1d;
        MSKint32t asub[nzi];
        double aval[nzi];

        for(int i = 0; i < nzi; i++)
        { 
            aval[i] = h * h / 2.0 * (2 * ( var_num_1d - i ) + 3);
            asub[i] = p * var_num_1d + i;    
        }

        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
        row_idx ++;
    }

    // #2.2   put the equality constraints in the end velocity
    for(int p = 0; p < 3; p++)
    {
        int nzi = var_num_1d;
        MSKint32t asub[nzi];
        double aval[nzi];

        for(int i = 0; i < nzi; i++)
        { 
            aval[i] = h;
            asub[i] = p * var_num_1d + i;    
        }

        r = MSK_putarow(task, row_idx, nzi, asub, aval);    
        row_idx ++;
    }

    //ROS_WARN("[Discrete Solver] Start stacking the objective");

    int NUMQNZ = 3 * ( var_num_1d + var_num_1d - 1);
    MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];
    double     qval[NUMQNZ];
    
    int idx = 0;
    for(int p = 0; p < 3; p++)
        for( int i = 0; i < var_num_1d; i ++ )
            for( int j = 0; j < var_num_1d; j ++ )
                if( i == j || i == (j + 1))
                {
                    qsubi[idx] = p * var_num_1d + i;   
                    qsubj[idx] = p * var_num_1d + j;  
                    qval[idx]  = H(i , j);    
                    idx ++;            
                }   

/*    int NUMQNZ = 3 *  var_num_1d;
    MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];
    double     qval[NUMQNZ];
    
    int idx = 0;
    for(int p = 0; p < 3; p++)
        for( int i = 0; i < var_num_1d; i ++ )
            for( int j = 0; j < var_num_1d; j ++ )
                if( i == j )
                {
                    qsubi[idx] = p * var_num_1d + i;   
                    qsubj[idx] = p * var_num_1d + j;  
                    qval[idx]  = 1.0;    
                    idx ++;            
                }   */

    if ( r== MSK_RES_OK )
         r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval); 
    
    if ( r==MSK_RES_OK ) 
         r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);

    ros::Time time_end1 = ros::Time::now();
    ROS_WARN("Time in variables stack is");
    cout<<time_end1 - time_start<<endl;
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
    cout<<time_end2 - time_start<<endl;


    if(!solve_ok){
      ROS_WARN("In solver, falied ");
      ROS_BREAK();
    }

    for(int i = 0; i < var_num; i++)
    {   
        if(i < var_num_1d)
            AccValue(0, i) = x_var[i];
        else if( i < 2 * var_num_1d)
            AccValue(1, i - var_num_1d) = x_var[i];
        else
            AccValue(2, i - 2 * var_num_1d) = x_var[i];

    }
    
    //cout<<"solution: \n"<<AccValue<<endl;
    return AccValue;
}