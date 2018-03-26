#include "pcd_trajectory/trajectory_generator.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include "pcd_trajectory/mosek.h"

using namespace std;    
using namespace Eigen;

#define POLYORDER 7
#define inf 1>>30

int TotalPoly_num;
int Segment_num;
int D1Poly_num;
int Poly_num;
int Equ_con_num;
int Inequ_con_num;

typedef Eigen::SparseMatrix<double> SMatrixXd;

SMatrixXd MapP_D;
SMatrixXd invMapP_D;

static void MSKAPI printstr(void *handle,
                            MSKCONST char str[])
{
  printf("%s",str);
} /* printstr */

static void 
          insertBlock(SMatrixXd &m, int lr, int lc, const SMatrixXd &n);

static pair<SMatrixXd, VectorXd> 
          combineRowsPr(const pair<SMatrixXd, VectorXd> &x, const pair<SMatrixXd, VectorXd> &y);

static SMatrixXd 
          combineRowsSM(const SMatrixXd &a, const SMatrixXd &b);

static SMatrixXd 
          combineRowsSM(const SMatrixXd &a, const SMatrixXd &b, const SMatrixXd &c);

static void 
          printSM(SMatrixXd &x);

static void 
          printSMtoFile(SMatrixXd x, string s);

static void 
          printSMListtoFile(vector<SMatrixXd> list, string s);

static void 
          printDMtoFile( MatrixXd &x, string s);

static void 
          printDVtoFile( VectorXd &x, string s);

static MatrixXd 
          getDM(SMatrixXd &x);

static void 
          printPr(pair<SMatrixXd, VectorXd> &pr);

static pair<SMatrixXd, VectorXd> 
          combineRowsPr(const pair<SMatrixXd, VectorXd> &x, const pair<SMatrixXd, VectorXd> &y, const pair<SMatrixXd, VectorXd> &z);

static pair<SMatrixXd, SMatrixXd>
          MapPtoD( VectorXd Time, int seg_num );

vector<double> findRoots( double start_time, double end_time, Eigen::VectorXd poly_coef );

int checkPosLimit( double root, double radius, Vector3d center, VectorXd poly );

void addExConstrain( vector<SMatrixXd> & ConQua, vector<SMatrixXd> & ConLin, VectorXd & bin,
                     vector<int> segment_no, vector<double> t_ex );

Eigen::MatrixXd solver( const int & Equ_con_num,  const int & Inequ_con_num,
                        const vector<SMatrixXd> & ConQua, const vector<SMatrixXd> & ConLin, const vector<SMatrixXd> & ConLin_z,
                        const SMatrixXd & Qobj, const SMatrixXd & Aeq, const VectorXd & beq, const VectorXd & bin, double z_bound );

vector<double>  
          TrajectoryGenerator::getCost() { return qp_cost; };

TrajectoryGenerator::TrajectoryGenerator(){  ex_ex_num = 0;}

TrajectoryGenerator::~TrajectoryGenerator(){}

Eigen::MatrixXd TrajectoryGenerator::PloyCoeffGeneration(
            const Eigen::MatrixXd &PathCorridor,
            const Eigen::MatrixXd &PathConnect,
            const Eigen::VectorXd &Radius,
            const Eigen::VectorXd &Path_Radius,
            const Eigen::VectorXd &Time,
            const Eigen::MatrixXd &vel,
            const Eigen::MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const double z_bound )
{
        _Time   = Time;
        _Path_Radius = Path_Radius;       
        ros::Time time_start = ros::Time::now();
        Eigen::MatrixXd PolyCoeff;
        Eigen::Vector3d StartPt  = PathConnect.row( 0 ); 
        Eigen::Vector3d StartVel = vel.row(0);
        Eigen::Vector3d EndPt    = PathConnect.row( PathConnect.rows() - 1 ); 

        Eigen::MatrixXd Path_new = PathConnect.block(1, 0, PathConnect.rows() - 2, PathConnect.cols() );
        _Path = Path_new;
        
        _PathCorridor = PathCorridor.block(1, 0, PathCorridor.rows() - 2, PathCorridor.cols() );
        /*ROS_WARN("check _PathCorridor");
        cout<<_PathCorridor<<endl;*/
        
/*        //Path.resize(Path_new.rows(), Path_new.cols());
        Path = Path_new;*/
//#####################################################################################################################
//Prepare for constaints and objective data stucture for Mosek solver .         

        Segment_num   = Time.size();
        
        TotalPoly_num = 3 * Segment_num * ( POLYORDER + 1 ) ; // The number of all polynomial coeff to be calculated
        D1Poly_num    = POLYORDER + 1;  // The number of unknowns in ONE single segment(in ONE single Dimension)
        Poly_num      = 3 * D1Poly_num; // The number of unknowns in ONE single segment(X + Y + Z Dimension)

/* ###### Objective Matirx ###### */

        SMatrixXd Qobj(TotalPoly_num, TotalPoly_num);
        /*SMatrixXd Qo_k(Poly_num, Poly_num);*/
        SMatrixXd Qo_k_1d(D1Poly_num, D1Poly_num);

        vector<SMatrixXd> ConQua; // List all the quadratic part in the inequality constraints
        vector<SMatrixXd> ConLin; // List all the linear part in the inequality constraints
        vector<SMatrixXd> ConLin_z; // List all the inequality constraints in z axis, to contrain the quadrotor fly above the ground

        SMatrixXd Aeq; // total equality matirx in the left hand side
        
        VectorXd beq; // total equality value in the right hand side
        VectorXd bin; // total Quadratic inequality value in the right hand side

/*        ros::Time time_begin = ros::Time::now();*/
/*        Qobj.setZero();*/

        for(int k = 0; k < Segment_num; k ++){
/*          Qo_k.setZero();*/
          Qo_k_1d.setZero();
            
          for( int i = 3; i < D1Poly_num; i ++ ){
            for( int j = 3; j < D1Poly_num; j ++ ){
              Qo_k_1d.insert( i, j ) = i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) * pow(Time(k), i + j - 5) / (i + j - 5);
              //cout<<"non-zero number in a 7 x 7 small block in Objective: "<<Qo_k_1d.nonZeros()<<endl;
            }
          }
          double coeff_dim = 1.0;
          for(int p = 0; p < 3; p ++ ){
            if(p == 2) coeff_dim = 1.0; //1.1
            insertBlock( Qobj, k*Poly_num + p * D1Poly_num, k*Poly_num + p * D1Poly_num, coeff_dim * Qo_k_1d);
          }  
      } 

/* ######  Linear Equality  ###### */
    
    /* ### p,v,a Constraints at Start point ### */
      // Aeq : The total equality constrain matrix, will be resize later 

      SMatrixXd Aeq_start  (9, TotalPoly_num);
      SMatrixXd Aeq_start_p(3, TotalPoly_num); /* ## Pos Constraints at start point ## */ 
      SMatrixXd Aeq_start_v(3, TotalPoly_num); /* ## Vel Constraints at start point ## */ 
      SMatrixXd Aeq_start_a(3, TotalPoly_num); /* ## Acc Constraints at start point ## */  
      
      Aeq_start  .setZero();
      Aeq_start_p.setZero();
      Aeq_start_v.setZero();
      Aeq_start_a.setZero();

      for(int i = 0; i < 3; i ++){
          Aeq_start_p.insert(i, 0 + i * D1Poly_num) = 1;
          Aeq_start_v.insert(i, 1 + i * D1Poly_num) = 1;
          Aeq_start_a.insert(i, 2 + i * D1Poly_num) = 2;
      }

      insertBlock(Aeq_start, 0, 0, Aeq_start_p);
      insertBlock(Aeq_start, 3, 0, Aeq_start_v);
      insertBlock(Aeq_start, 6, 0, Aeq_start_a);

    /* ### p,v,a Constraints at End point ### */
      SMatrixXd Aeq_end  (9, TotalPoly_num);
      SMatrixXd Aeq_end_p(3, TotalPoly_num); /* ## Pos Constraints at end point ## */ 
      SMatrixXd Aeq_end_v(3, TotalPoly_num); /* ## Vel Constraints at end point ## */ 
      SMatrixXd Aeq_end_a(3, TotalPoly_num); /* ## Acc Constraints at end point ## */                  

      Aeq_end  .setZero();
      Aeq_end_p.setZero();
      Aeq_end_v.setZero();
      Aeq_end_a.setZero();

      for(int i = 0; i < 3; i ++){
        for(int j = 0; j < D1Poly_num; j ++ ){
          Aeq_end_p.insert(i, TotalPoly_num - (3 - i) * D1Poly_num + j ) = pow(Time(Segment_num - 1), j);
          if( j > 0 )
            Aeq_end_v.insert(i, TotalPoly_num - (3 - i) * D1Poly_num + j ) = j * pow(Time(Segment_num - 1), j - 1);
          if( j > 1 )
            Aeq_end_a.insert(i, TotalPoly_num - (3 - i) * D1Poly_num + j ) = j *(j - 1) * pow(Time(Segment_num - 1), j - 2);
        }
      }

      insertBlock(Aeq_end, 0, 0, Aeq_end_p);
      insertBlock(Aeq_end, 3, 0, Aeq_end_v);
      insertBlock(Aeq_end, 6, 0, Aeq_end_a);

      //ROS_WARN("start and end constrain is OK ");
    /* ### Continuous Constraints at junction point ### */
      SMatrixXd Aeq_con  ( 9 * (Segment_num - 1), TotalPoly_num );
      SMatrixXd Aeq_con_p( 3 * (Segment_num - 1), TotalPoly_num );
      SMatrixXd Aeq_con_v( 3 * (Segment_num - 1), TotalPoly_num );
      SMatrixXd Aeq_con_a( 3 * (Segment_num - 1), TotalPoly_num );
      
      Aeq_con  .setZero();
      Aeq_con_p.setZero();
      Aeq_con_v.setZero();
      Aeq_con_a.setZero();

      for(int k = 0; k < Segment_num - 1; k ++){
        
        SMatrixXd Aeq_seg_con_p(3, TotalPoly_num);
        SMatrixXd Aeq_seg_con_v(3, TotalPoly_num);
        SMatrixXd Aeq_seg_con_a(3, TotalPoly_num);
        
        Aeq_seg_con_p.setZero();
        Aeq_seg_con_v.setZero();
        Aeq_seg_con_a.setZero();

        for(int i = 0; i < 3; i ++){
          for(int j = 0; j < D1Poly_num; j ++){
            Aeq_seg_con_p.insert( i, i * D1Poly_num + k * Poly_num + j ) = pow( Time(k), j );
            if(j > 0)
              Aeq_seg_con_v.insert( i, i * D1Poly_num + k * Poly_num + j ) = j * pow( Time(k), j -1 );
            if(j > 1)
              Aeq_seg_con_a.insert( i, i * D1Poly_num + k * Poly_num + j ) = j * (j - 1) * pow( Time(k), j -2 );
          }
          /*for(int j = 1; j < D1Poly_num; j ++){
          }
          for(int j = 2; j < D1Poly_num; j ++){
          }*/
          
          Aeq_seg_con_p.insert( i, i * D1Poly_num + k * Poly_num + 3 * D1Poly_num + 0 ) = -1;
          Aeq_seg_con_v.insert( i, i * D1Poly_num + k * Poly_num + 3 * D1Poly_num + 1 ) = -1;
          Aeq_seg_con_a.insert( i, i * D1Poly_num + k * Poly_num + 3 * D1Poly_num + 2 ) = -2;
        }

        insertBlock(Aeq_con_p, 3 * k, 0, Aeq_seg_con_p);
        insertBlock(Aeq_con_v, 3 * k, 0, Aeq_seg_con_v);
        insertBlock(Aeq_con_a, 3 * k, 0, Aeq_seg_con_a);
      }

      insertBlock(Aeq_con, 0, 0, Aeq_con_p);
      insertBlock(Aeq_con, 3 * (Segment_num - 1), 0, Aeq_con_v);
      insertBlock(Aeq_con, 6 * (Segment_num - 1), 0, Aeq_con_a);

      Aeq.resize(2 * 9 + 9 * (Segment_num - 1), TotalPoly_num);
      Aeq.setZero();

      assert(Aeq.rows() == (Aeq_start.rows() + Aeq_end.rows() + Aeq_con.rows() ));

      insertBlock(Aeq, 0,  0, Aeq_start);
      insertBlock(Aeq, 9,  0, Aeq_end);
      insertBlock(Aeq, 18, 0, Aeq_con);

/*      ROS_WARN("ARe you OK ?");*/
      
/*      MatrixXd Aeq_dense = getDM(Aeq);*/
      /*ROS_WARN("ARe you OK2 ?");
      cout<<Aeq_dense.rows()<<","<<Aeq_dense.cols()<<endl;*/

      SMatrixXd Aeq_tmp = Aeq;
      Aeq.setZero();
      Aeq = Aeq_tmp;// * invMapP_D;

      /*ros::Time time_aft_linear = ros::Time::now();*/
/*      ROS_WARN("check time in linear equality");
      cout<< time_aft_linear - time_aft_mapping<<endl;*/
      //printSMtoFile(Aeq, "EqualityConstraints.txt");
      
      /*ROS_WARN("Linear Equality constrain is OK ");*/

  /* ######  Quadratic Inequality  ###### */

      for( int k = 0; k < Segment_num - 1; k ++ ){
        SMatrixXd Q_qc_k  ( TotalPoly_num, TotalPoly_num );
        SMatrixXd Q_qc_k1D( D1Poly_num, D1Poly_num);
        SMatrixXd F_qc_k  ( 1, TotalPoly_num);
        SMatrixXd F_qc_k1D( 1, D1Poly_num);
        SMatrixXd F_z_k   ( 1, TotalPoly_num);

      /* ### Quadratic Part ### */
          for(int i = 0; i < D1Poly_num; i ++){
            for(int j = 0; j < D1Poly_num; j ++){
              Q_qc_k1D.insert(i, j) = pow( Time( k ), i + j );
            }
          }
        
        for(int i = 0; i < 3; i ++)
          insertBlock(Q_qc_k, i * D1Poly_num + 3 * k * D1Poly_num, i * D1Poly_num + 3 * k * D1Poly_num, 
            Q_qc_k1D);

      SMatrixXd Q_qc_k_tmp = Q_qc_k;
      
      ConQua.push_back(Q_qc_k_tmp);

      /*ROS_WARN("OK here ... ");*/

      /* ### Linear Part ### */
        for(int i = 0; i < D1Poly_num; i ++ )
            F_qc_k1D.insert(0, i) = pow( Time( k ), i );

        insertBlock(F_qc_k, 0, k * Poly_num,                  -2.0 * Path_new( k, 0 ) * F_qc_k1D ); 
        insertBlock(F_qc_k, 0, k * Poly_num + D1Poly_num,     -2.0 * Path_new( k, 1 ) * F_qc_k1D );  
        insertBlock(F_qc_k, 0, k * Poly_num + 2 * D1Poly_num, -2.0 * Path_new( k, 2 ) * F_qc_k1D ); 

        SMatrixXd F_qc_k_tmp = F_qc_k; // * invMapP_D;
        ConLin.push_back(F_qc_k_tmp);

        //F_z_k.insert(0, k * Poly_num + 2 * D1Poly_num) = 1.0;
/*        insertBlock(F_z_k, 0, k * Poly_num + 2 * D1Poly_num, -1.0 * F_qc_k1D ); 
        SMatrixXd F_z_k_tmp = F_z_k; // * invMapP_D;
        ConLin_z.push_back(F_z_k_tmp);*/
      }

/*      ros::Time time_aft_qc_con = ros::Time::now();
      ROS_WARN("check time in Quadratic constriants");
      cout<< time_aft_qc_con - time_aft_linear<<endl;*/

/*        ROS_WARN("ConQua size is : ");
        cout<<ConQua.size()<<endl;

        ROS_WARN("ConLin size is : ");
        cout<<ConLin.size()<<endl;

        ROS_WARN("ConLin_z size is : ");
        cout<<ConLin_z.size()<<endl;

      ofstream quad_quadFile ("Quad_Quad_constraints.txt");
      if (quad_quadFile.is_open() ){
        for(auto ptr:ConQua)
          quad_quadFile << " quadratic part in inequality constraints is: \n " << ptr << "\n"<<endl;
        quad_quadFile.close(); } else cout << "Unable to open file";

      ofstream quad_linFile ("Quad_lin_constraints.txt");
      if (quad_linFile.is_open() ){
        for(auto ptr:ConLin)
          quad_linFile << " linear part in inequality constraints is: \n " << ptr <<  "\n"<< endl;
        quad_linFile.close(); } else cout << "Unable to open file";

      ofstream quad_lin_zFile ("Quad_lin_z_constraints.txt");
      if (quad_lin_zFile.is_open() ){
        for(auto ptr:ConLin_z)
          quad_lin_zFile << " linear part in inequality constraints is: \n " << ptr <<  "\n"<< endl;
        quad_lin_zFile.close(); } else cout << "Unable to open file";*/

/*      printSMListtoFile(ConQua,   "Quad_Quad_constraints.txt" );
      printSMListtoFile(ConLin,   "Quad_lin_constraints.txt"  );
      printSMListtoFile(ConLin_z, "Quad_lin_z_constraints.txt");*/

/*      ROS_WARN("Quadratic Inequality constrain is OK ");*/

/* ######  Value in Linear Equality  ###### */      
        
        beq.resize( Aeq.rows(), 1 );
        beq.setZero();
        
        for(int i = 0 ; i < 3; i ++ ){
          beq( i )     = StartPt(i); 
          beq( i + 3)  = StartVel(i);
          beq( i + 9 ) = EndPt(i);
        }

      //printDVtoFile(beq, "beq.txt");

/*      ROS_WARN("Right hand side Value in Equality constrain is OK ");*/
/* ######  Value in Quadratic Inequality  ###### */      
        
        bin.resize( Segment_num - 1, 1 );
        bin.setZero();
        
        for(int i = 0; i < Segment_num - 1; i ++ )
          bin( i ) = Radius(i) * Radius(i) - Path_new(i, 0) * Path_new(i, 0) - Path_new(i, 1) * Path_new(i, 1) - Path_new(i, 2) * Path_new(i, 2) ;
      //printDVtoFile(bin, "bin.txt");       
      /*ROS_WARN("Right hand side Value in Quadratic Inequality constrain is OK ");*/

      Equ_con_num   = beq.size();
      /*int Inequ_con_num = 2 * ( bin.size() - 1 );*/
      Inequ_con_num = bin.size();//2 * bin.size();

      ROS_WARN("[TrajectoryGenerator] prepare for the convex solver");
      PolyCoeff = solver( Equ_con_num, Inequ_con_num, ConQua, ConLin, ConLin_z, Qobj, Aeq, beq, bin, z_bound );   
      if(PolyCoeff.rows() == 3 && PolyCoeff.cols() == 3){
            ROS_WARN("Solver Wrong, Exit");
      }

      ROS_WARN("[TrajectoryGenerator] Solver finish");
/*      cout<<"non-zero number in Objective: "<<Qobj.nonZeros()<<endl;

      cout<<"beq size :"<<beq.size()<<endl;
      cout<<"bin size :"<<bin.size()<<endl;
      cout<<"Inequ_con_num: "<<Inequ_con_num<<endl;*/

      //ROS_WARN("prepare work is OK");

/*      int cycle_ct = 0;
      for( int i = 0; i < 5; i ++ ){
          ROS_WARN("Time bef");
          cout<<ros::Time::now()<<endl;
          PolyCoeff = solver( Equ_con_num, Inequ_con_num, ConQua, ConLin, ConLin_z, Qobj, Aeq, beq, bin, z_bound );
          if(PolyCoeff.rows() == 3 && PolyCoeff.cols() == 3){
            ROS_WARN("Solver Wrong, Exit");
            break;
          }

          pair< vector<int>, vector<double> > check_result = checkEx(PolyCoeff);
          vector<int>     seg_no  = check_result.first;
          vector<double>  ex_time = check_result.second;

          ex_ex_num = seg_no.size();
          ROS_WARN("total number of extreme out of limit is : ");
          cout<<ex_ex_num<<endl;

          if(ex_ex_num == 0)
            break;

          else{
            addExConstrain( ConQua, ConLin, bin, seg_no, ex_time );
            Inequ_con_num += ex_time.size();
            //Inequ_con_num += 0;
          }
          cycle_ct ++;
      }

    ROS_WARN("check cycle num");
    cout<<cycle_ct<<endl;

    if(ex_ex_num != 0){
      ROS_WARN("can't compress extreme within time, Exit");
      PolyCoeff = MatrixXd::Identity(3,3);
    } */

    ros::Time time_end = ros::Time::now();
    ROS_WARN("time consume is :");
    cout<<time_end - time_start<<endl;

    return PolyCoeff;
}

Eigen::MatrixXd solver( const int & Equ_con_num,  const int & Inequ_con_num,
                        const vector<SMatrixXd> & ConQua, const vector<SMatrixXd> & ConLin, const vector<SMatrixXd> & ConLin_z,
                        const SMatrixXd & Qobj, const SMatrixXd & Aeq, const VectorXd & beq, const VectorXd & bin, double z_bound )
{
/*            ROS_WARN("Time aft");
          cout<<ros::Time::now()<<endl;*/
    //int ieqn = Inequ_con_num / 2;
    int con_num = Equ_con_num + Inequ_con_num;// / 2;

    #define NUMVAR TotalPoly_num // Number of unknowns X.
    #define NUMCON con_num       // Number of all constraints

    MatrixXd PolyCoeff;

    ros::Time time_opt = ros::Time::now();
    double x_var[ NUMVAR];
    
    MSKrescodee  r; 
    /* ## define a container for constaints boundary and boundkey ## */ 
    /* ## dataType in the double pair is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    
    for(int i =0; i < Inequ_con_num; i++ ){
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_UP, make_pair( -MSK_INFINITY, bin(i) ) ); // # cb_ie means: constriants boundary of inequality constrain      
        con_bdk.push_back(cb_ie); 
    }

/*    // add a group of constraints: the height of the quadrotor( in z axis) must > 0.0 m
    for(int i =0; i < Inequ_con_num / 2; i++ ){
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_UP, make_pair( -MSK_INFINITY, -z_bound ) ); // # cb_ie means: constriants boundary of inequality constrain      
        con_bdk.push_back(cb_ie); 
    }*/

    for(int i = 0; i < Equ_con_num; i ++ ){ 
        pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair(    beq(i),     beq(i) ) ); // # cb_eq means: constriants boundary of equality constrain
        con_bdk.push_back(cb_eq);
    }

    cout<<"Equ_con_num: "<<Equ_con_num<<endl;
    cout<<"Inequ_con_num: "<<Inequ_con_num<<endl;
    cout<<"con_num: "<<con_num<<endl;
    cout<<"size "<<con_bdk.size()<<endl;

    /* ## define a container for unknowns boundary and boundkey ## */ 
    /* ## dataType in one tuple is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 
    
    for(int i = 0; i < TotalPoly_num; i ++ ){
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair(    -MSK_INFINITY,      +MSK_INFINITY ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)
        var_bdk.push_back(vb_x);
    } 
  
    MSKint32t  j,i; 
    MSKenv_t   env; 
    MSKtask_t  task; 
 
    // Create the mosek environment. 
    r = MSK_makeenv( &env, NULL ); 
  
    if ( r != MSK_RES_OK )
      ROS_BREAK();
    
    // Create the optimization task. 
    r = MSK_maketask(env,NUMCON,NUMVAR,&task); 

//######################################################################
    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);
    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-5);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 1e-2 );
    //MSK_putdouparam (task, MSK_DPAR_PRESOLVE_TOL_X, 1e-10 );
    
    /*MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 1e-2);
    MSK_putdouparam (task, MSK_DPAR_QCQO_REFORMULATE_REL_DROP_TOL, 1e-5);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_STEP_SIZE, 1e-3);*/
//######################################################################

      r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
 
    // Append 'NUMCON' empty constraints. 
     //The constraints will initially have no bounds. 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendcons(task,NUMCON);  

    // Append 'NUMVAR' variables. 
    // The variables will initially be fixed at zero (x=0). 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendvars(task,NUMVAR); 

    for(j=0; j<NUMVAR && r == MSK_RES_OK; ++j) 
    { 
      // Set the bounds on variable j. 
      //  blx[j] <= x_j <= bux[j] 
      if(r == MSK_RES_OK) 
        r = MSK_putvarbound(task, 
                            j,                            // Index of variable. 
                            var_bdk[j].first,             // Bound key.
                            var_bdk[j].second.first,      // Numerical value of lower bound.
                            var_bdk[j].second.second );   // Numerical value of upper bound.      
    } 
    
    //ROS_WARN("Start the Linear Matrix A Part ");
    // Input Linear matrix A 
    vector< pair< pair<int,int>, double > > tmp_value;
    
    for( size_t j = 0; j < ConLin.size(); j++ ){
      for ( int k = 0; k < ConLin[j].outerSize(); ++k )
        for ( SMatrixXd::InnerIterator it(ConLin[j],k); it; ++it){
          tmp_value.push_back( make_pair( make_pair(it.row() + j, it.col() ), 1.0 * it.value() ) );   // linear part in quadratic inequality
      }
    }

/*    for( size_t j = 0; j < ConLin_z.size(); j++ ){
      for ( int k = 0; k < ConLin_z[j].outerSize(); ++k )
        for ( SMatrixXd::InnerIterator it(ConLin_z[j],k); it; ++it){
          tmp_value.push_back( make_pair( make_pair(it.row() + ConLin.size() + j, it.col() ), it.value() ) );   // linear part in z axis
      }
    }
*/
    for ( int k = 0; k < Aeq.outerSize(); ++k )
      for ( SMatrixXd::InnerIterator it(Aeq,k); it; ++it){
        //tmp_value.push_back( make_pair( make_pair(it.row() + ConLin_z.size() + ConLin.size(), it.col() ), it.value() ) );   // linear equality 
        tmp_value.push_back( make_pair( make_pair(it.row() + ConLin.size(), it.col() ), it.value() ) );   // linear equality 
    }

    int NUMANZ = tmp_value.size();
    MSKint32t   asubi[NUMANZ], asubj[NUMANZ];
    double      aval[NUMANZ];

    for(int i =0; i < NUMANZ; i++ ){
       asubi[i] = tmp_value[i].first.first;   
       asubj[i] = tmp_value[i].first.second;  
       aval[i]  = tmp_value[i].second;
       //cout<<"col: "<<asubj[i]<<" row: "<<asubi[i]<<" Value in A : "<<aval[i]<<endl;
    }

    if(r == MSK_RES_OK)
        r = MSK_putaijlist( task, NUMANZ, asubi, asubj, aval ); 

    // Set the bounds on constraints. 
    //   for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] 
    for( i = 0; i < NUMCON && r == MSK_RES_OK; i++ ) {
    /*cout<<"NO."<<i<<" Bound key: "<<con_bdk[i].first<<" ,lower bound: "<< con_bdk[i].second.first << " ,upper bound: "<<con_bdk[i].second.second<<endl;*/
        r = MSK_putconbound(task, 
                            i,                            // Index of constraint. 
                            con_bdk[i].first,             // Bound key.
                            con_bdk[i].second.first,      // Numerical value of lower bound.
                            con_bdk[i].second.second );   // Numerical value of upper bound. 
    }
    
    //ROS_WARN("Start the Objective part ");
    if ( r== MSK_RES_OK ) { 
    //  The lower triangular part of the Q^o 
    //  matrix in the objective is specified. 
    tmp_value.clear();
    for ( int k = 0; k < Qobj.outerSize(); ++k ){
      for ( SMatrixXd::InnerIterator it(Qobj, k); it; ++it){
        if( (it.row() > it.col()) || (it.row() == it.col()) )
          tmp_value.push_back( make_pair( make_pair(it.row(), it.col() ), 2.0 * it.value() ) );   // col index 
      }
    }

    int NUMQNZ = tmp_value.size();
    MSKint32t   qsubi[NUMQNZ], qsubj[NUMQNZ];
    double      qval[NUMQNZ];
    
    for(int i =0; i < NUMQNZ; i++ ){
       qsubi[i] = tmp_value[i].first.first;   
       qsubj[i] = tmp_value[i].first.second;  
       qval[i]  = tmp_value[i].second;
       //cout<<"col: "<<qsubj[i]<<" row: "<<qsubi[i]<<" Value in Qobj: "<<qval[i]<<endl;
    }
        // Input the Q^o for the objective. 
        r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval); 
    } 
    
    if ( r==MSK_RES_OK ) 
        r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);       
    //ROS_WARN("Finish the Objective part ");

/*    r = MSK_toconic(task);    
    if(r != MSK_RES_OK)
      ROS_WARN("toconic fail ");
    cout<<r<<endl;*/

    //ROS_WARN("Start the Quadratic part ");
    // The lower triangular part of the Q^0 matrix in the first constraint is specified.     
    /*cout<<"Inequ_con_num: "<<Inequ_con_num<<endl;*/

    double coeff = 2.0;
    for(int j = 0; j < Inequ_con_num; j ++ ){
      
      tmp_value.clear();
      int idx = j;

      for ( int k = 0; k < ConQua[idx].outerSize(); ++k ){
        for ( SMatrixXd::InnerIterator it(ConQua[idx], k); it; ++it){
          if( (it.row() > it.col()) || (it.row() == it.col()) )
            tmp_value.push_back( make_pair( make_pair(it.row(), it.col() ), coeff * it.value() ) );   // col index 
        }
      }
    int numQiNZ = tmp_value.size();

    /*cout<<"numQiNZ: "<<numQiNZ<<endl;*/

    MSKint32t   qisubi[numQiNZ], qisubj[numQiNZ];
    double      qival[numQiNZ];
    
    for(int i =0; i < numQiNZ; i++ ){
       qisubi[i] = tmp_value[i].first.first;   
       qisubj[i] = tmp_value[i].first.second;  
       qival[i]  = tmp_value[i].second; 
      /*cout<<"col: "<<qisubj[i]<<" row: "<<qisubi[i]<<" Value in Q"<<j<<" : "<<qival[i]<<endl;*/
    } 

    // Put Qj in constraint with index j.  
       /* cout<<" constrain NO. :"<<j<<endl;*/
        r = MSK_putqconk(task, j, numQiNZ, qisubi, qisubj, qival);  
    }
    
/*    r = MSK_toconic( task ); */
    
    /*if ( r!=MSK_RES_OK ) 
        ROS_BREAK();*/
    
    bool solve_ok = false;
      if ( r==MSK_RES_OK ) 
      { 
        ROS_WARN("Prepare to solve the problem ");   

        MSKrescodee trmcode; 

        // Run optimizer 
        r = MSK_optimizetrm(task,&trmcode); 
 
        // Print a summary containing information 
        //   about the solution for debugging purposes
        
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
/*              printf("Optimal primal solution\n"); 
              
              for(j=0; j<NUMVAR; ++j) 
                printf("x[%d]: %e\n", j, x_var[j]); */
            
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
    
    MSK_deletetask(&task); 
    MSK_deleteenv(&env); 

    ros::Time time_stall = ros::Time::now();
    ROS_WARN("time consume in optimize is :");
    cout<<time_stall - time_opt<<endl;

    PolyCoeff.resize(Segment_num, Poly_num);

    if(!solve_ok){
      MatrixXd poly_fail = MatrixXd::Identity(3,3);
      ROS_WARN("In solver, falied ");
      return poly_fail;
    }

    VectorXd d_var(TotalPoly_num);
    for(int i = 0; i < TotalPoly_num; i++)
      d_var(i) = x_var[i];

    VectorXd p_var = d_var;
    //invMapP_D * d_var;
    //invMapP_D * d_var; 

    for(int i = 0; i < TotalPoly_num; i++)
        PolyCoeff(i / Poly_num, i % Poly_num) = p_var(i);

    return PolyCoeff;
}

pair< vector<int>, vector<double> > TrajectoryGenerator::checkEx( Eigen::MatrixXd Poly )
{   
    int pcn = POLYORDER + 1; // number of trajectory polynomial in one segment, 1 dimension (X, Y, Z seperately) pcn: poly coeff number
    int frn = POLYORDER * 2 + 1; // number of polynomial used to check roots in one segment, 3 dimension ( norm 2), frn: finding roots number

    vector<int>    seg_no;
    vector<double> ex_time;

    VectorXd poly_seg; // for polynomial coefficients in one single segment
    VectorXd polyX(pcn), polyY(pcn), polyZ(pcn);
    VectorXd polyCheck(frn);
    for(int i = 0; i < Poly.rows(); i++ )
    {   
        double x0 = _PathCorridor(i, 0), y0 = _PathCorridor(i, 1), z0 = _PathCorridor(i, 2);

        poly_seg = Poly.row(i);
        polyX    = poly_seg.segment(0 * pcn, pcn);
        polyY    = poly_seg.segment(1 * pcn, pcn);
        polyZ    = poly_seg.segment(2 * pcn, pcn);
        
        polyCheck(0)   =   polyX(0) * polyX(0) + polyY(0) * polyY(0) + polyZ(0) * polyZ(0) 
                         - 2 * ( x0 * polyX(0) + y0 * polyY(0) + z0 * polyZ(0) ) ;

        polyCheck(2)   =   polyX(1) * polyX(1) + polyY(1) * polyY(1) + polyZ(1) * polyZ(1)
                         + 2 * ( polyX(0) * polyX(2) + polyY(0) * polyY(2) + polyZ(0) * polyZ(2) )
                         - 2 * ( x0 * polyX(2) + y0 * polyY(2) + z0 * polyZ(2) ) ;

        polyCheck(4)   =   polyX(2) * polyX(2) + polyY(2) * polyY(2) + polyZ(2) * polyZ(2)
                         + 2 * ( polyX(0) * polyX(4) + polyY(0) * polyY(4) + polyZ(0) * polyZ(4) )
                         + 2 * ( polyX(1) * polyX(3) + polyY(1) * polyY(3) + polyZ(1) * polyZ(3) )
                         - 2 * ( x0 * polyX(4) + y0 * polyY(4) + z0 * polyZ(4) ) ;

        polyCheck(6)   =   polyX(3) * polyX(3) + polyY(3) * polyY(3) + polyZ(3) * polyZ(3)
                         + 2 * ( polyX(1) * polyX(5) + polyY(1) * polyY(5) + polyZ(1) * polyZ(5) )
                         + 2 * ( polyX(2) * polyX(4) + polyY(2) * polyY(4) + polyZ(2) * polyZ(4) );

        polyCheck(8)   =   polyX(4) * polyX(4) + polyY(4) * polyY(4) + polyZ(4) * polyZ(4)
                         + 2 * ( polyX(3) * polyX(5) + polyY(3) * polyY(5) + polyZ(3) * polyZ(5) );

        polyCheck(10)  =   polyX(5) * polyX(5) + polyY(5) * polyY(5) + polyZ(5) * polyZ(5);
        
        polyCheck(1)   =   2 * ( polyX(0) * polyX(1) + polyY(0) * polyY(1) + polyZ(0) * polyZ(1) )
                         - 2 * ( x0 * polyX(1) + y0 * polyY(1) + z0 * polyZ(1) ) ;

        polyCheck(3)   =   2 * ( polyX(0) * polyX(3) + polyY(0) * polyY(3) + polyZ(0) * polyZ(3) )
                         + 2 * ( polyX(1) * polyX(2) + polyY(1) * polyY(2) + polyZ(1) * polyZ(2) )
                         - 2 * ( x0 * polyX(3) + y0 * polyY(3) + z0 * polyZ(3) ) ;

        polyCheck(5)   =   2 * ( polyX(0) * polyX(5) + polyY(0) * polyY(5) + polyZ(0) * polyZ(5) )
                         + 2 * ( polyX(1) * polyX(4) + polyY(1) * polyY(4) + polyZ(1) * polyZ(4) )
                         + 2 * ( polyX(2) * polyX(3) + polyY(2) * polyY(3) + polyZ(2) * polyZ(3) )
                         - 2 * ( x0 * polyX(5) + y0 * polyY(5) + z0 * polyZ(5) ) ;

        polyCheck(7)   =   2 * ( polyX(2) * polyX(5) + polyY(2) * polyY(5) + polyZ(2) * polyZ(5) )
                         + 2 * ( polyX(3) * polyX(4) + polyY(3) * polyY(4) + polyZ(3) * polyZ(4) );
        
        polyCheck(9)   =   2 * ( polyX(4) * polyX(5) + polyY(4) * polyY(5) + polyZ(4) * polyZ(5) );
        

/*        ROS_WARN("The segment of polynomial coeffs are ");
        cout<<poly_seg<<endl;
        ROS_WARN("The polynomial to be checked extreme");
        cout<<polyCheck<<endl;*/

        //cout<<"i: "<<i<<endl;
        vector<double> rts = findRoots(0.0, _Time(i), polyCheck);
        double   rad    = _Path_Radius(i);
        Vector3d center = _PathCorridor.row(i);

        for(int j = 0; j < int(rts.size()); j++ )
            if( checkPosLimit( rts[j], rad, center, polyCheck ) ){
              seg_no.push_back(i);
              ex_time.push_back(rts[j]);
            }
        
    }

    return make_pair(seg_no, ex_time);
    // use findRoots function to find extreme point( where derivative equal to 0 ) at pos, vel and acc
}

double eps = 1e-5;
vector<double> findRoots( double start_time, double end_time, Eigen::VectorXd poly_coef )
{   
    // Step 0: formulate the derivative of input polynomial function
    int der_size = poly_coef.size() - 1;
    VectorXd poly_der(der_size);
    
    for(int i = 0; i < der_size; i++)
      poly_der(i) = (i + 1) * poly_coef(i + 1);

/*    ROS_WARN("check the derivative of Polynomial");
    cout<<poly_der<<endl;
*/    
    MatrixXd Der_comp(der_size, der_size);
    Der_comp.setZero();

    for(int i = 0; i < der_size; i++){
        
        if(i > 0)
          Der_comp(i, i - 1)        = 1.0;
        
        Der_comp(i, der_size - 1) =  - poly_der(i);
      }

/*      ROS_WARN("check companian matrix of poly derivative");
      cout<<Der_comp<<endl;*/

              // the complex roots
      VectorXcd eig  = Der_comp.eigenvalues();

      /*ROS_WARN("check eigenvalues number is : %lu", eig.size());*/

      vector<double> rts;

      for (int i = 0; i < eig.size(); i++){
            if (abs(eig(i).imag()) < eps){   
                double rt_real = eig(i).real();
                if( rt_real >= start_time && rt_real <= end_time )
                  rts.push_back( rt_real );
            }
      }

/*      ROS_WARN("check extreme point within interval is : %lu", rts.size());
      for(auto ptr:rts)
        cout<<"root : "<<ptr<<endl;*/

      return rts;

    // Step 1: formulate the companian matrix of polynomial's derivative in one single segment. ouptut: a list of all resulting extreme points ( t_ex )
    // Step 2: check whether t_ex in the extreme list is within this segment of time ( start_time -- end_time )
    // Step 3: return the roots in the time interval
}

int checkPosLimit( double root, double radius, Vector3d center, VectorXd poly ) // return signal: 1 for extreme exceeds limits,  1 indicates exceed, 0 indicates within constrain
{    
    // Calculate extreme value( pos, vel, acc ) in each extreme time per segment, check whether its pos is within the ball path, and whether its vel and acc is below dynamic maxium constrain.
/*    ROS_WARN("check center and radius");
    cout<<center<<endl;
    cout<<radius<<endl;*/

    double lvalue = center(0) * center(0) + center(1) * center(1) + center(2) * center(2);
    double t_ex = root;
    /*cout<<"lvalue is : "<<lvalue<<endl;*/
    
    
    for(int i = 0; i < poly.size(); i++)
      lvalue += poly(i) * pow( t_ex, i );

/*    cout<<"rvalue is : "<<radius * radius<<endl;
    cout<<"lvalue is : "<<lvalue<<endl;*/
    
    if( lvalue > radius * radius )
      return 1;
    else
      return 0;

}

void TrajectoryGenerator::addExConstrain( vector<SMatrixXd> & ConQua, vector<SMatrixXd> & ConLin, VectorXd & bin,
                                          vector<int> segment_no, vector<double> t_ex )
{
    // add extreme constraint at pos(t_ex) extreme where out of the ball path ( Quadratic constraint )
    // add dynamic constraint at vel(t_ex) and acc(t_ex) where exceeds the maximum vel and acc ( Linear constraint )
    ROS_WARN("In function addExConstrain");
    cout<<t_ex.size()<<endl;
    cout<<segment_no.size()<<endl;

   /* ROS_WARN("check extreme segment");
    for(auto ptr:segment_no)
      cout<<ptr<<endl;

    ROS_WARN("check extreme time");
    for(auto ptr:t_ex)
      cout<<ptr<<endl;*/

    assert(t_ex.size() == segment_no.size());
    int addConsNum = t_ex.size();

    VectorXd bin_tmp = bin;
    bin.resize( bin_tmp.size() + t_ex.size(), 1 );
    
    VectorXd bin_add(t_ex.size());

    for(int i = 0; i < addConsNum; i ++ ){
          int idx = segment_no[i];
           bin_add( i ) =  ( _Path_Radius(idx) / 1.1 ) * (_Path_Radius(idx) / 1.1) 
                          -  _PathCorridor(idx, 0     ) *  _PathCorridor(idx, 0    ) 
                          -  _PathCorridor(idx, 1     ) *  _PathCorridor(idx, 1    ) 
                          -  _PathCorridor(idx, 2     ) *  _PathCorridor(idx, 2    );

          /*ROS_WARN("check bin_add .");
          cout<<bin_add<<endl;*/
    }

    bin.segment(0, bin_tmp.size()) = bin_tmp;
    bin.segment(bin_tmp.size(), bin_add.size()) =bin_add;


    for( int k = 0; k < addConsNum; k ++ ){
        int idx = segment_no[k];
        SMatrixXd Q_qc_add  ( TotalPoly_num, TotalPoly_num );
        SMatrixXd Q_qc_add1D( D1Poly_num, D1Poly_num);
        SMatrixXd F_qc_add  ( 1, TotalPoly_num);
        SMatrixXd F_qc_add1D( 1, D1Poly_num);

      /* ### Quadratic Part ### */
          for(int i = 0; i < D1Poly_num; i ++){
            for(int j = 0; j < D1Poly_num; j ++){
              Q_qc_add1D.insert(i, j) = pow( t_ex[k], i + j );
            }
          }
        
        for(int i = 0; i < 3; i ++)
          insertBlock(Q_qc_add, i * D1Poly_num + 3 * idx * D1Poly_num, i * D1Poly_num + 3 * idx * D1Poly_num, Q_qc_add1D);

    SMatrixXd Q_qc_add_tmp = invMapP_D.transpose() * Q_qc_add * invMapP_D;
    //invMapP_D.transpose() * Q_qc_add * invMapP_D;
      
    ConQua.push_back(Q_qc_add_tmp);

    for(int i = 0; i < D1Poly_num; i ++ )
      F_qc_add1D.insert(0, i) = pow( t_ex[k], i );

    insertBlock(F_qc_add, 0, idx * Poly_num, -2.0 * _PathCorridor( idx, 0 ) * F_qc_add1D ); 
    insertBlock(F_qc_add, 0, idx * Poly_num + D1Poly_num, -2.0 * _PathCorridor( idx, 1 ) * F_qc_add1D );  
    insertBlock(F_qc_add, 0, idx * Poly_num + 2 * D1Poly_num, -2.0 * _PathCorridor( idx, 2 ) * F_qc_add1D ); 

    SMatrixXd F_qc_add_tmp = F_qc_add * invMapP_D;
    ConLin.push_back(F_qc_add_tmp);

    }

}

static void insertBlock( SMatrixXd &m, int lr, int lc, const SMatrixXd &n )
{
    for (int k = 0; k < n.outerSize(); k++)
        for (SMatrixXd::InnerIterator it(n,k); it; ++it)
            m.insert(lr + it.row(), lc + it.col()) = it.value();
}

static pair<SMatrixXd, VectorXd> combineRowsPr( const pair<SMatrixXd, VectorXd> &x, const pair<SMatrixXd, VectorXd> &y)
{
        size_t rows_x = x.first.rows(), rows_y = y.first.rows();
        size_t cols = x.first.cols();

        assert( (unsigned)y.first.cols() == cols );

        SMatrixXd m(rows_x + rows_y, cols);
        m.reserve(x.first.nonZeros() + y.first.nonZeros());

        VectorXd b = VectorXd::Zero( rows_x + rows_y);
        
        insertBlock(m, 0, 0, x.first);
        insertBlock(m, rows_x, 0, y.first);

        b << x.second , y.second;
        
        return make_pair(m,b);
}

    static SMatrixXd combineRowsSM(const SMatrixXd &a, const SMatrixXd &b)
    {
        size_t cols = a.cols(), rows_a = a.rows(), rows_b = b.rows();
        assert((unsigned)b.cols() == cols);

        SMatrixXd m(rows_a + rows_b, cols);
        m.reserve(a.nonZeros() + b.nonZeros());
        insertBlock(m, 0, 0, a);
        insertBlock(m, rows_a, 0, b);

        return m;
    }

    static SMatrixXd combineRowsSM(const SMatrixXd &a, const SMatrixXd &b, const SMatrixXd &c)
    {
        return combineRowsSM(combineRowsSM(a, b), c);
    }

    static void printSM(SMatrixXd &x)
    {
        MatrixXd m(x.rows(), x.cols());
        for (int i = 0; i < x.rows(); i++)
            for (int j = 0; j < x.cols(); j++)
                m(i, j) = x.coeffRef(i, j);
        cout << m <<endl;
    }

    static void printSMtoFile(SMatrixXd x, string s)
    {
        MatrixXd m(x.rows(), x.cols());
        for (int i = 0; i < x.rows(); i++)
            for (int j = 0; j < x.cols(); j++)
                m(i, j) = x.coeffRef(i, j);

      ofstream File(s);
      if (File.is_open() ){
        File << "Sparse Matirx is: \n " << m << endl;
        File.close(); } else cout << "Unable to open file";

        ROS_WARN("Print sparse matrix is OK ");
    }

    static void printSMListtoFile(vector<SMatrixXd> list, string s)
    {
      ofstream File (s);
      if (File.is_open() ){
        for(auto ptr:list)
        {   
            SMatrixXd x = ptr;
            MatrixXd m(x.rows(), x.cols());
            for (int i = 0; i < x.rows(); i++)
              for (int j = 0; j < x.cols(); j++)
                  m(i, j) = x.coeffRef(i, j);
          File << " quadratic part in inequality constraints is: \n " << m << "\n"<<endl;
        }
      
      File.close(); } else cout << "Unable to open file";

      ROS_WARN("Print sparse matrix list is OK ");
    
    }

    static void printDMtoFile( MatrixXd &x, string s)
    {
        ofstream objFile (s);
        if (objFile.is_open() ){
          objFile << s <<" is: \n " << x << endl;
          objFile.close(); } else cout << "Unable to open file";
      
    }

    static void printDVtoFile( VectorXd &x, string s)
    {
        ofstream objFile (s);
        if (objFile.is_open() ){
          objFile << s <<" is: \n " << x << endl;
          objFile.close(); } else cout << "Unable to open file";
      
    }

    static MatrixXd getDM(SMatrixXd &x)
    {
        MatrixXd m = MatrixXd::Zero(x.rows(), x.cols());

        for (int k = 0; k < x.outerSize(); k++)
            for (SMatrixXd::InnerIterator it(x,k); it; ++it)
                m(it.row(), it.col()) = it.value();
        return m;
    }

    static void printPr(pair<SMatrixXd, VectorXd> &pr)
    {
        MatrixXd m(pr.first.rows(), pr.first.cols() +1);
        for (int i = 0; i < pr.first.rows(); i++)
        {
            for (int j = 0; j < pr.first.cols(); j++)
                m(i, j) = pr.first.coeffRef(i, j);
            m(i, pr.first.cols()) = pr.second(i);
        }
        clog<<m<<endl;
    }

    static pair<SMatrixXd, VectorXd> combineRowsPr(
            const pair<SMatrixXd, VectorXd> &x, 
            const pair<SMatrixXd, VectorXd> &y,
            const pair<SMatrixXd, VectorXd> &z)
    {
        return combineRowsPr(x, combineRowsPr(y, z));
    }

static pair<SMatrixXd, SMatrixXd> MapPtoD( VectorXd Time, int seg_num )
{     
      int TotalPoly_num = 3 * seg_num * ( POLYORDER + 1 );
      int D1Poly_num    = POLYORDER + 1;
      int Poly_num      = 3 * D1Poly_num;

      /*ROS_WARN("Generate P2DMapping ");*/
      SMatrixXd P2DMapping(TotalPoly_num, TotalPoly_num);
      SMatrixXd invP2DMapping(TotalPoly_num, TotalPoly_num);

/*      SMatrixXd Map_k (Poly_num, Poly_num);
      SMatrixXd invMap_k (Poly_num, Poly_num);
      */
      SMatrixXd Map_1d(D1Poly_num, D1Poly_num);
      SMatrixXd invMap_1d(D1Poly_num, D1Poly_num);

/*      P2DMapping.setZero();
      invP2DMapping.setZero();*/

      for(int k = 0; k < seg_num; k ++){
/*            Map_k.setZero();
            invMap_k.setZero();*/
            Map_1d.setZero();
            invMap_1d.setZero();
            
            for( int i = 0; i < D1Poly_num; i ++ ){
                for( int j = i; j < D1Poly_num; j ++ ){
                  int factorial_j = 1, factorial_ji = 1;
                  
                  for(int p = 1; p < j + 1; p ++)     factorial_j *= p;
                  for(int p = 1; p < j - i + 1; p ++) factorial_ji *= p;  
                  
                  Map_1d.insert( i, j ) = factorial_j / factorial_ji * pow(Time(k), j - i);

                }
            }
            MatrixXd Map_1d_dense = getDM(Map_1d);
            MatrixXd invMap_1d_dense = Map_1d_dense.inverse();

            for(int i = 0; i < invMap_1d_dense.rows(); i++)
                for(int j = 0; j < invMap_1d_dense.cols(); j++)
                    invMap_1d.insert( i, j ) =  invMap_1d_dense(i,j);   
            
/*            for(int p = 0; p < 3; p ++ ){
                insertBlock( Map_k, p * D1Poly_num, p * D1Poly_num, Map_1d );
                insertBlock( invMap_k, p * D1Poly_num, p * D1Poly_num, invMap_1d );
            }
           
           insertBlock( P2DMapping, k*Poly_num, k*Poly_num, Map_k);
           insertBlock( invP2DMapping, k*Poly_num, k*Poly_num, invMap_k);*/

            for(int p = 0; p < 3; p ++ ){
                insertBlock( P2DMapping, k*Poly_num + p * D1Poly_num, k*Poly_num + p * D1Poly_num, Map_1d );
                insertBlock( invP2DMapping, k*Poly_num + p * D1Poly_num, k*Poly_num + p * D1Poly_num, invMap_1d );
            }
      } 

      /*printSMtoFile(invP2DMapping, "invP2DMapping");*/
      return make_pair(P2DMapping, invP2DMapping);
}
