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

//#define POLYORDER 5
#define inf 1>>30

int TotalPoly_num, Segment_num, D1Poly_num, Poly_num;
int Equ_con_num, Inequ_con_num, Inequ_con_v_num, Inequ_con_a_num;

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
          MapPtoD( VectorXd Time, int seg_num, int POLYORDER );

vector<double> findRoots( double start_time, double end_time, Eigen::VectorXd poly_coef );

void addPosExConstrain( vector<SMatrixXd> & ConQua, vector<SMatrixXd> & ConLin, VectorXd & bin,
                     vector<int> segment_no, vector<double> t_ex );

Eigen::MatrixXd solver( const int & Equ_con_num,  const int & Inequ_con_num, const int & Inequ_con_v_num, const int & Inequ_con_a_num,
                        const vector<SMatrixXd> & ConQua, const vector<SMatrixXd> & ConLin,
                        const SMatrixXd & Qobj, const SMatrixXd & Aeq, const SMatrixXd & Alin_v, const SMatrixXd & Alin_a, 
                        const VectorXd & beq, const VectorXd & bin,
                        const int & type );

vector<double>  TrajectoryGeneratorPoly::getObjective() { return qp_cost; };

TrajectoryGeneratorPoly::TrajectoryGeneratorPoly()
{
    _delta_pos = 0.0; 
    _delta_vel = 0.0; 
    _delta_acc = 0.0; 

    _d_eps = 0.025;
}

TrajectoryGeneratorPoly::~TrajectoryGeneratorPoly(){}

std::vector<Eigen::MatrixXd> TrajectoryGeneratorPoly::PloyCoeffGeneration(
            const Eigen::MatrixXd &Path,
            const Eigen::VectorXd &Radius,
            const Eigen::VectorXd &Time,
            const Eigen::MatrixXd &vel,
            const Eigen::MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int POLYORDER,
            double scale, double w_l, double w_s )
{  
#define ENFORCE_VEL  0 // whether or not adding extra constraints for ensuring the velocity feasibility
#define ENFORCE_ACC  0 // whether or not adding extra constraints for ensuring the velocity feasibility

        _scale = scale;
        ros::Time time_start = ros::Time::now();
        Eigen::MatrixXd PolyCoeff;
        std::vector<Eigen::MatrixXd> polycoeffList;

        Eigen::Vector3d StartPt  = Path.row(0); 
        Eigen::Vector3d StartVel = vel.row(0);
        Eigen::Vector3d StartAcc = acc.row(0) * _scale;

        Eigen::Vector3d EndPt    = Path.row( Path.rows() - 1 ); 
        
        _Path   = Path.block(1, 0, Path.rows() - 2, Path.cols() );
        _Time   = Time;
        _Radius = Radius;     

        _max_v  = maxVel;
        _max_a  = maxAcc * _scale;
        /*cout<<"maximum velocity: "<<_max_v<<endl;
        cout<<"maximum acceleration: "<<_max_a<<endl;*/
//#####################################################################################################################
//Prepare for constaints and objective data stucture for Mosek .         
        Segment_num   = Radius.size();
        
        TotalPoly_num = 3 * Segment_num * ( POLYORDER + 1 ) ; // The number of all polynomial coeff to be calculated
        D1Poly_num    = POLYORDER + 1;  // The number of unknowns in ONE single segment(in ONE single Dimension)
        Poly_num      = 3 * D1Poly_num; // The number of unknowns in ONE single segment(X + Y + Z Dimension)

/* ###### Objective Matirx ###### */
        SMatrixXd Qobj(TotalPoly_num, TotalPoly_num);
        /*SMatrixXd Qo_k(Poly_num, Poly_num);*/
        SMatrixXd Qo_k_1d(D1Poly_num, D1Poly_num);
        SMatrixXd Qo_k_1d_smo(D1Poly_num, D1Poly_num);
        SMatrixXd Qo_k_1d_arc(D1Poly_num, D1Poly_num);

        vector<SMatrixXd> ConQua; // List all the quadratic part in the inequality constraints
        vector<SMatrixXd> ConLin; // List all the linear part in the inequality constraints
        
        SMatrixXd Aeq;    // total equality matirx in the left hand side
        SMatrixXd Alin_v; // total linear equality matirx in the left hand side for adding velocity constraints
        SMatrixXd Alin_a; // total linear equality matirx in the left hand side for adding acceleration constraints
        
        VectorXd beq; // total equality value in the right hand side
        VectorXd bin; // total Quadratic inequality value in the right hand side

        for(int k = 0; k < Segment_num; k ++)
        {
          Qo_k_1d.setZero();
          Qo_k_1d_smo.setZero();
          Qo_k_1d_arc.setZero();
            
          for( int i = 3; i < D1Poly_num; i ++ ){
            for( int j = 3; j < D1Poly_num; j ++ ){
              Qo_k_1d_smo.insert( i, j ) = w_s * (double)i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) * pow(Time(k), i + j - 5) / (double)(i + j - 5);
            }
          }

        /*for i=2:n_poly_perSeg
          for j=2:n_poly_perSeg
            Ql(i,j)=(i-1)*(j-1)/(i+j-3);
          end
        end*/
          
          for( int i = 1; i < D1Poly_num; i ++ ){
            for( int j = 1; j < D1Poly_num; j ++ ){
              Qo_k_1d_arc.insert( i, j ) = w_l * (double)i * j * pow(Time(k), i + j - 1) / (double)(i + j - 1);
            }
          }

          Qo_k_1d = Qo_k_1d_arc + Qo_k_1d_smo;

          for(int p = 0; p < 3; p ++ )
            insertBlock( Qobj, k*Poly_num + p * D1Poly_num, k*Poly_num + p * D1Poly_num, Qo_k_1d);
      } 

      SMatrixXd QobjD = Qobj;//invMapP_D.transpose() * Qobj * invMapP_D;
      Qobj.setZero();      
      Qobj = QobjD;
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

      for(int k = 0; k < Segment_num - 1; k ++)
      {  
          SMatrixXd Aeq_seg_con_p(3, TotalPoly_num);
          SMatrixXd Aeq_seg_con_v(3, TotalPoly_num);
          SMatrixXd Aeq_seg_con_a(3, TotalPoly_num);
          
          Aeq_seg_con_p.setZero();
          Aeq_seg_con_v.setZero();
          Aeq_seg_con_a.setZero();

          for(int i = 0; i < 3; i ++)
          {
              for(int j = 0; j < D1Poly_num; j ++)
              {
                  Aeq_seg_con_p.insert( i, i * D1Poly_num + k * Poly_num + j ) = pow( Time(k), j );
                  if(j > 0)
                    Aeq_seg_con_v.insert( i, i * D1Poly_num + k * Poly_num + j ) = j * pow( Time(k), j -1 );
                  if(j > 1)
                    Aeq_seg_con_a.insert( i, i * D1Poly_num + k * Poly_num + j ) = j * (j - 1) * pow( Time(k), j -2 );
              }
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

      Alin_v.resize( 3 * Segment_num, TotalPoly_num );
      Alin_a.resize( 3 * Segment_num, TotalPoly_num );

      Alin_v.setZero();
      Alin_a.setZero();

      for(int k = 0; k < Segment_num; k ++)
      {  
          SMatrixXd Alin_v_seg(3, TotalPoly_num);
          SMatrixXd Alin_a_seg(3, TotalPoly_num);
          
          Alin_v_seg.setZero();
          Alin_a_seg.setZero();
          
          for(int i = 0; i < 3; i ++)
          {
              Alin_v_seg.insert( i, i * D1Poly_num + k * Poly_num + 1 ) = 1.0;
              Alin_a_seg.insert( i, i * D1Poly_num + k * Poly_num + 2 ) = 2.0;
          }

          insertBlock(Alin_v, 3 * k, 0, Alin_v_seg);
          insertBlock(Alin_a, 3 * k, 0, Alin_a_seg);
      }

      //ROS_WARN("continuity constraint is OK ");
      /* ######  Quadratic Inequality  ###### */
      for( int k = 0; k < Segment_num - 1; k ++ )
      {
          for(int p = 0; p < 2; p++)
          {
              SMatrixXd Q_qc_k  ( TotalPoly_num, TotalPoly_num );
              SMatrixXd Q_qc_k1D( D1Poly_num, D1Poly_num);

              /* ### Quadratic Part ### */
              for(int i = 0; i < D1Poly_num; i ++){
                for(int j = 0; j < D1Poly_num; j ++){
                  Q_qc_k1D.insert(i, j) = pow( Time( k ), i + j );
                }
              }
              
              for(int i = 0; i < 3; i ++)
                insertBlock(Q_qc_k, i * D1Poly_num + 3 * k * D1Poly_num, i * D1Poly_num + 3 * k * D1Poly_num, Q_qc_k1D);

              SMatrixXd Q_qc_k_tmp = Q_qc_k;//invMapP_D.transpose() * Q_qc_k * invMapP_D;             

              ConQua.push_back(Q_qc_k_tmp);
                
              SMatrixXd F_qc_k  ( 1, TotalPoly_num);
              SMatrixXd F_qc_k1D( 1, D1Poly_num);
              /* ### Linear Part ### */
              for(int i = 0; i < D1Poly_num; i ++ )
                F_qc_k1D.insert(0, i) = pow( Time( k ), i );

              insertBlock(F_qc_k, 0, (k) * Poly_num,                  -2.0 * _Path( k + p, 0 ) * F_qc_k1D ); 
              insertBlock(F_qc_k, 0, (k) * Poly_num + D1Poly_num,     -2.0 * _Path( k + p, 1 ) * F_qc_k1D );  
              insertBlock(F_qc_k, 0, (k) * Poly_num + 2 * D1Poly_num, -2.0 * _Path( k + p, 2 ) * F_qc_k1D ); 

              SMatrixXd F_qc_k_tmp = F_qc_k;// * invMapP_D;

              ConLin.push_back(F_qc_k_tmp);
          }
      }

      /* ######  Value in Linear Equality  ###### */       
      beq.resize( Aeq.rows(), 1 );
      beq.setZero();
      
      for(int i = 0 ; i < 3; i ++ )
      {
        beq( i )     = StartPt(i); 
        beq( i + 3)  = StartVel(i);
        beq( i + 6)  = StartAcc(i);
        beq( i + 9 ) = EndPt(i);
      }

      /* ######  Value in Quadratic Inequality  ###### */      
      bin.resize( 2 * (Segment_num - 1), 1 );
      bin.setZero();

      for(int i = 0; i < Segment_num - 1; i ++ )
      {
        bin( 2*i )   = pow(_Radius(i)  , 2.0) - pow(_Path(i, 0),   2.0) - pow(_Path(i, 1),   2.0) - pow(_Path(i, 2), 2.0);
        bin( 2*i+1 ) = pow(_Radius(i+1), 2.0) - pow(_Path(i+1, 0), 2.0) - pow(_Path(i+1, 1), 2.0) - pow(_Path(i+1, 2), 2.0);
      }
      
      //cout<<"check bin list: \n"<<bin<<endl;
      Equ_con_num   = beq.size();
      Inequ_con_num = bin.size();
      Inequ_con_v_num = 3 * Segment_num;
      Inequ_con_a_num = 3 * Segment_num;

      int cycle_ct = 0;
      int loop_num = 10;
      bool is_no_extraExtremum = true;
      
      int pos_ex_num = 0;
      int vel_ex_num = 0;
      int acc_ex_num = 0;

      for( int l = 0; l < loop_num; l ++ )
      {    
          is_no_extraExtremum = true;
          _delta_pos += _d_eps;
          _delta_vel += _d_eps;
          _delta_acc += _d_eps;

          if(ENFORCE_VEL == 0 && ENFORCE_ACC == 0)
            PolyCoeff = solver( Equ_con_num, Inequ_con_num, Inequ_con_v_num, Inequ_con_a_num, vel_ex_num, acc_ex_num, ConQua, ConLin, Qobj, Aeq, Alin_v, Alin_a, beq, bin, 0 );
          else if(ENFORCE_VEL == 1 && ENFORCE_ACC == 0)
            PolyCoeff = solver( Equ_con_num, Inequ_con_num, Inequ_con_v_num, Inequ_con_a_num, vel_ex_num, acc_ex_num, ConQua, ConLin, Qobj, Aeq, Alin_v, Alin_a, beq, bin, 1 );
          else if(ENFORCE_VEL == 1 && ENFORCE_ACC == 1)
            PolyCoeff = solver( Equ_con_num, Inequ_con_num, Inequ_con_v_num, Inequ_con_a_num, vel_ex_num, acc_ex_num, ConQua, ConLin, Qobj, Aeq, Alin_v, Alin_a, beq, bin, 2 );
          else
          { 
            ROS_ERROR("ERROR in polynomial solver type selection");
            ROS_BREAK();
          }

          polycoeffList.push_back(PolyCoeff);
          if(PolyCoeff.rows() == 3 && PolyCoeff.cols() == 3)
          {
            ROS_WARN("Solver Wrong, Exit");
            break;
          }

          /*check extremum in position of each piece trajectory*/
          pair< vector<int>, vector<double> > check_p_result = checkPosEx(PolyCoeff, POLYORDER);
          vector<int>     p_seg_no  = check_p_result.first;
          vector<double>  p_ex_time = check_p_result.second;
          pos_ex_num = p_seg_no.size();
          
          if(pos_ex_num != 0)
          { 
              is_no_extraExtremum = false;
              addPosExConstrain( ConQua, ConLin, bin, p_seg_no, p_ex_time );
              Inequ_con_num += pos_ex_num;
          }

          if(ENFORCE_VEL)
          {   
              vel_ex_num = 0;
              for(int i = 0; i < 3; i++)
              {   
                  //ROS_WARN("enforcing vel check at axis %d", i);
                  int poly_num = POLYORDER + 1;
                  MatrixXd PolyCoeff_i = PolyCoeff.block(0, i * poly_num, PolyCoeff.rows(), poly_num);

                  pair< vector<int>, vector<double> > check_v_result = checkVelEx(PolyCoeff_i, POLYORDER);
                  vector<int>     v_seg_no  = check_v_result.first;
                  vector<double>  v_ex_time = check_v_result.second;
                  int vel_ex_num_i = v_seg_no.size();
                  vel_ex_num += vel_ex_num_i;

                  if(vel_ex_num_i != 0)
                  {   
                      is_no_extraExtremum = false;
                      addVelExConstrain(Alin_v, v_seg_no, v_ex_time, i );
                      Inequ_con_v_num += vel_ex_num_i;
                  }
              }
          } 

          if(ENFORCE_ACC)
          {   
              acc_ex_num = 0;
              for(int i = 0; i < 3; i++)
              {   
                  //ROS_WARN("enforcing acc check at axis %d", i);
                  int poly_num = POLYORDER + 1;
                  MatrixXd PolyCoeff_i = PolyCoeff.block(0, i * poly_num, PolyCoeff.rows(), poly_num);

                  pair< vector<int>, vector<double> > check_a_result = checkAccEx(PolyCoeff_i, POLYORDER);
                  vector<int>     a_seg_no  = check_a_result.first;
                  vector<double>  a_ex_time = check_a_result.second;
                  int acc_ex_num_i = a_seg_no.size();
                  acc_ex_num += acc_ex_num_i;
                  
                  if(acc_ex_num_i != 0)
                  {   
                      is_no_extraExtremum = false;
                      addAccExConstrain( Alin_a, a_seg_no, a_ex_time, i );
                      Inequ_con_a_num += acc_ex_num_i;
                  }
              }
          }

          if(is_no_extraExtremum == true)
            break;
          
          cycle_ct ++;
    }

    ROS_WARN("check cycle num");
    cout<<cycle_ct<<endl;

    if(pos_ex_num != 0 || vel_ex_num != 0 || acc_ex_num != 0)
    {
        ROS_WARN("can't compress extreme within time, Exit");
        PolyCoeff = MatrixXd::Identity(3,3);
        polycoeffList.push_back(PolyCoeff);
    }

    ros::Time time_end = ros::Time::now();
    ROS_WARN("time consume is :");
    cout<<time_end - time_start<<endl;

    return polycoeffList;
}

Eigen::MatrixXd TrajectoryGeneratorPoly::solver
                      ( const int & Equ_con_num,  const int & Inequ_con_num, const int & Inequ_con_v_num, const int & Inequ_con_a_num,
                        const int & vel_ex_num, const int & acc_ex_num, const vector<SMatrixXd> & ConQua, const vector<SMatrixXd> & ConLin,
                        const SMatrixXd & Qobj, const SMatrixXd & Aeq, const SMatrixXd & Alin_v, const SMatrixXd & Alin_a, const VectorXd & beq, const VectorXd & bin,
                        const int & type )
{
  int con_num = Equ_con_num + Inequ_con_num;                                
  MatrixXd PolyCoeff;

  ros::Time time_opt = ros::Time::now();
  
  MSKrescodee  r;
  double primalobj; 
  /* ## define a container for constaints boundary and boundkey ## */ 
  /* ## dataType in the double pair is : boundary type, lower bound, upper bound ## */
  vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
  
  for(int i = 0; i < Inequ_con_num; i++ ){
      //pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_UP, make_pair( -MSK_INFINITY, bin( (i + 1)/2 ) ) ); // # cb_ie means: constriants boundary of inequality constrain
      pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_UP, make_pair( -MSK_INFINITY, bin(i) ) ); // # cb_ie means: constriants boundary of inequality constrain
      
      con_bdk.push_back(cb_ie); 
  }

  for(int i = 0; i < Equ_con_num; i ++ ){ 
      pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair(    beq(i),     beq(i) ) ); // # cb_eq means: constriants boundary of equality constrain
      con_bdk.push_back(cb_eq);
  }

  /*cout<<"vel_ex_num: "<<vel_ex_num<<endl;
  cout<<"acc_ex_num: "<<acc_ex_num<<endl;*/

  if(type == 1 || type == 2)
  { 
      con_num += Inequ_con_v_num;
      for(int i = 0; i < Inequ_con_v_num; i++ )
      {   
          pair<MSKboundkeye, pair<double, double> > cb_ie;

          if( i < (Inequ_con_v_num - vel_ex_num) )
            cb_ie = make_pair( MSK_BK_RA, make_pair( - _max_v, + _max_v ) );
          else
            cb_ie = make_pair( MSK_BK_RA, make_pair( - (_max_v + _delta_vel), + (_max_v + _delta_vel) ) );

          con_bdk.push_back(cb_ie); 
      }
  }

  if(type == 2)
  { 
      con_num += Inequ_con_a_num;
      for(int i = 0; i < Inequ_con_a_num; i++ )
      {   
          pair<MSKboundkeye, pair<double, double> > cb_ie;

          if( i < (Inequ_con_a_num - acc_ex_num) )
              cb_ie = make_pair( MSK_BK_RA, make_pair( - _max_a, + _max_a ) );      
          else
              cb_ie = make_pair( MSK_BK_RA, make_pair( - (_max_a + _delta_acc), + (_max_a + _delta_acc) ) );      

          con_bdk.push_back(cb_ie); 
      }
  }

  #define NUMVAR TotalPoly_num // Number of unknowns X.
  #define NUMCON con_num       // Number of all constraints
  double x_var[ NUMVAR];

  /* ## define a container for unknowns boundary and boundkey ## */ 
  /* ## dataType in one tuple is : boundary type, lower bound, upper bound ## */
  vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 
  
  for(int i = 0; i < TotalPoly_num; i ++ ){
      pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( -MSK_INFINITY, +MSK_INFINITY ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)
      var_bdk.push_back(vb_x);
  } 
  
  MSKint32t    j,i; 
  MSKenv_t     env; 
  MSKtask_t    task; 
 
  // Create the mosek environment. 
  r = MSK_makeenv( &env, NULL ); 
  
  if ( r != MSK_RES_OK ) 
    ROS_BREAK();
    
    /*ROS_WARN("Are you OK?1");*/
    // Create the optimization task. 
    r = MSK_maketask(env,NUMCON,NUMVAR,&task); 

    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);
    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-5);
/*    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 1e-2);*/

//######################################################################
   /* MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 1e-4);
    MSK_putdouparam (task, MSK_DPAR_QCQO_REFORMULATE_REL_DROP_TOL, 1e-5);*/
//######################################################################


    /*MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_STEP_SIZE, 1e-3);*/

    // do not print
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
    
    // Input Linear matrix A 
    vector< pair< pair<int,int>, double > > tmp_value;
    
    for( size_t j = 0; j < ConLin.size(); j++ ){
      for ( int k = 0; k < ConLin[j].outerSize(); ++k )
        for ( SMatrixXd::InnerIterator it(ConLin[j],k); it; ++it){
          tmp_value.push_back( make_pair( make_pair(it.row() + j, it.col() ), it.value() ) );   // linear part in quadratic inequality
      }
    }

    for ( int k = 0; k < Aeq.outerSize(); ++k )
      for ( SMatrixXd::InnerIterator it(Aeq,k); it; ++it)
        tmp_value.push_back( make_pair( make_pair(it.row() + ConLin.size(), it.col() ), it.value() ) ); 

    if( type == 1 || type == 2 )
        for ( int k = 0; k < Alin_v.outerSize(); ++k )
          for ( SMatrixXd::InnerIterator it(Alin_v,k); it; ++it)
            tmp_value.push_back( make_pair( make_pair(it.row() + ConLin.size() + Equ_con_num, it.col() ), it.value() ) ); 

    if( type == 2 )
        for ( int k = 0; k < Alin_a.outerSize(); ++k )
          for ( SMatrixXd::InnerIterator it(Alin_a,k); it; ++it)
            tmp_value.push_back( make_pair( make_pair(it.row() + ConLin.size() + Equ_con_num + Inequ_con_v_num, it.col() ), it.value() ) ); 

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
    
    if ( r==MSK_RES_OK ) 
    { 
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
/*    ROS_WARN("Objective part ");*/

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

    // The lower triangular part of the Q^0 matrix in the first constraint is specified.     
    /*cout<<"Inequ_con_num: "<<Inequ_con_num<<endl;*/
    for(int j = 0; j < Inequ_con_num; j ++ ){
      
      tmp_value.clear();
      int idx = j;

      for ( int k = 0; k < ConQua[idx].outerSize(); ++k ){
        for ( SMatrixXd::InnerIterator it(ConQua[idx], k); it; ++it){
          if( (it.row() > it.col()) || (it.row() == it.col()) )
            tmp_value.push_back( make_pair( make_pair(it.row(), it.col() ), 2.0 * it.value() ) );   // col index 
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
    
    bool solve_ok = false;
    if ( r==MSK_RES_OK ) 
    { 
      //ROS_WARN("Prepare to solve the problem ");   
      /*ROS_WARN("toconic success ");*/
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
            
          r = MSK_getxx(task, MSK_SOL_ITR, x_var); 

          r = MSK_getprimalobj(task, MSK_SOL_ITR, &primalobj);

          qp_cost.push_back( primalobj / pow(_scale, 3) );

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
  
    MSK_deletetask(&task); 
    MSK_deleteenv(&env); 

    ros::Time time_stall = ros::Time::now();
    /*ROS_WARN("time consume in optimize is :");
    cout<<time_stall - time_opt<<endl;*/

    PolyCoeff.resize(Segment_num, Poly_num);

    if(!solve_ok){
      MatrixXd poly_fail = MatrixXd::Identity(3,3);
      ROS_WARN("In solver, falied ");
      return poly_fail;
    }

    VectorXd d_var(TotalPoly_num);
    for(int i = 0; i < TotalPoly_num; i++)
      d_var(i) = x_var[i];

    VectorXd p_var = d_var;//invMapP_D * d_var; //MapP_D_dense.inverse() * d_var;

    for(int i = 0; i < TotalPoly_num; i++)
        PolyCoeff(i / Poly_num, i % Poly_num) = p_var(i);

    return PolyCoeff;
}

pair< vector<int>, vector<double> > TrajectoryGeneratorPoly::checkPosEx( Eigen::MatrixXd Poly, const int POLYORDER )
{   
    int pcn = POLYORDER + 1; // number of trajectory polynomial in one segment, 1 dimension (X, Y, Z seperately) pcn: poly coeff number
    int frn = pcn * 2 - 1;   // number of polynomial used to check roots in one segment, 3 dimension ( norm 2), frn: finding roots number

    vector<int>    seg_no;
    vector<double> ex_time;

    VectorXd poly_seg; // for polynomial coefficients in one single segment
    VectorXd polyX(pcn), polyY(pcn), polyZ(pcn);
    VectorXd polyCheck(frn);

    for(int k = 0; k < Poly.rows(); k++ )
    {   
        polyCheck = VectorXd::Zero(frn);
        double x0 = _Path(k, 0), y0 = _Path(k, 1), z0 = _Path(k, 2);

        poly_seg = Poly.row(k);
        polyX    = poly_seg.segment(0 * pcn, pcn);
        polyY    = poly_seg.segment(1 * pcn, pcn);
        polyZ    = poly_seg.segment(2 * pcn, pcn);
        
        for(int i = 0; i < pcn; i++)
          for(int j = 0; j < pcn; j++)
              polyCheck(i+j) += polyX(i) * polyX(j) + polyY(i) * polyY(j) + polyZ(i) * polyZ(j);

        for(int i  = 0; i < pcn; i++)
          polyCheck(i) += -2 * (x0 * polyX(i) + y0 * polyY(i) + z0 * polyZ(i));

        vector<double> rts = findRoots(0.0, _Time(k), polyCheck);
        double   rad    = _Radius(k);
        Vector3d center = _Path.row(k);

        for(int j = 0; j < int(rts.size()); j++ )
            if( checkPosLimit( rts[j], rad, center, polyCheck ) ){
              seg_no.push_back(k);
              ex_time.push_back(rts[j]);
            }
    }

    return make_pair(seg_no, ex_time);
    // use findRoots function to find extreme point( where derivative equal to 0 ) at pos, vel and acc
}

pair< vector<int>, vector<double> > TrajectoryGeneratorPoly::checkVelEx( Eigen::MatrixXd Poly, const int POLYORDER )
{   
    int pcn = POLYORDER + 1; // number of trajectory's polynomial in one segment, 1 dimension (X, Y, Z seperately) pcn: poly coeff number
    int vcn = POLYORDER;     // number of trajectory's velocity polynomial in one segment, 1 dimension (X, Y, Z seperately) pcn: poly coeff number
    int frn = vcn;           // number of polynomial used to check roots in one segment, 3 dimension ( norm 2), frn: finding roots number

    vector<int>    seg_no;
    vector<double> ex_time;

    VectorXd poly(pcn);
    VectorXd velCheck(frn);
    
    for(int k = 0; k < Poly.rows(); k++ )
    {   
        poly = Poly.row(k);
        for(int i = 1; i < pcn; i++)
            velCheck(i-1) = i * poly(i);

        vector<double> rts = findRoots(0.0, _Time(k), velCheck);
        double v_m = _max_v;

        for(int j = 0; j < int(rts.size()); j++ )
        {   
            //cout<<" No. "<<k<<" segment"<<endl;
            if( checkVelLimit( rts[j], v_m, velCheck ) )
            {
              seg_no.push_back(k);
              ex_time.push_back(rts[j]);
            }
        }
    }

    return make_pair(seg_no, ex_time);
}

pair< vector<int>, vector<double> > TrajectoryGeneratorPoly::checkAccEx( Eigen::MatrixXd Poly, const int POLYORDER )
{   
    int pcn = POLYORDER + 1; // number of trajectory's polynomial in one segment, 1 dimension (X, Y, Z seperately) pcn: poly coeff number
    int acn = POLYORDER - 1;     // number of trajectory's velocity polynomial in one segment, 1 dimension (X, Y, Z seperately) pcn: poly coeff number
    int frn = acn;           // number of polynomial used to check roots in one segment, 3 dimension ( norm 2), frn: finding roots number

    vector<int>    seg_no;
    vector<double> ex_time;

    VectorXd poly(pcn);
    VectorXd accCheck(frn);
    
    for(int k = 0; k < Poly.rows(); k++ )
    {   
        poly = Poly.row(k);
        for(int i = 2; i < pcn; i++)
            accCheck(i-2) = i * (i - 1) * poly(i);
          
        vector<double> rts = findRoots(0.0, _Time(k), accCheck);
        double a_m = _max_a;

        for(int j = 0; j < int(rts.size()); j++ )
        {    
            //cout<<" No. "<<k<<" segment"<<endl;
            if( checkAccLimit( rts[j], a_m, accCheck ) )
            {
              seg_no.push_back(k);
              ex_time.push_back(rts[j]);
            }
        }
    }

    return make_pair(seg_no, ex_time);
}

double eps = 1e-5;
vector<double> findRoots( double start_time, double end_time, Eigen::VectorXd poly_coef )
{   
    // Step 0: formulate the derivative of input polynomial function
    int der_size = poly_coef.size() - 1;
    VectorXd poly_der(der_size);
    
    for(int i = 0; i < der_size; i++)
    {
      poly_der(i) = (i + 1) * poly_coef(i + 1);
    }

/*    ROS_WARN("check the derivative of Polynomial");
    cout<<poly_der<<endl;*/
    
    MatrixXd Der_comp(der_size - 1, der_size - 1);
    Der_comp.setZero();

    for(int i = 0; i < der_size - 1; i++){
        
        if(i > 0)
          Der_comp(i, i - 1)        = 1.0;
        
        Der_comp(0, i) =  - poly_der(der_size - 1 - i - 1) / poly_der(der_size - 1);
      }

/*      ROS_WARN("check companian matrix of poly derivative");
      cout<<Der_comp<<endl;*/

              // the complex roots
      VectorXcd eig  = Der_comp.eigenvalues();

/*      ROS_WARN("check eigenvalues number is : %lu", eig.size());
      cout<<"check eigen value: "<<eig<<endl;
*/
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

int TrajectoryGeneratorPoly::checkPosLimit( double root, double radius, Vector3d center, VectorXd poly ) // return signal: 1 for extreme exceeds limits,  1 indicates exceed, 0 indicates within constrain
{    
    // Calculate extreme value( pos, vel, acc ) in each extreme time per segment, check whether its pos is within the ball path, and whether its vel and acc is below dynamic maxium constrain.
    double lvalue = center(0) * center(0) + center(1) * center(1) + center(2) * center(2);
    double t_ex = root;
    
    for(int i = 0; i < poly.size(); i++)
      lvalue += poly(i) * pow( t_ex, i );
    
    if( lvalue > radius * radius )
      return 1;
    else
      return 0;

}

int TrajectoryGeneratorPoly::checkVelLimit( double root, double v_m, VectorXd poly ) 
{    
    double lvalue = 0.0;
    double t_ex = root;
    
    for(int i = 0; i < poly.size(); i++)
      lvalue += poly(i) * pow( t_ex, i );
    
    //cout<<"extreme velocity value : "<<lvalue<<" , occurs at time : "<<t_ex<<endl;

    if( fabs(lvalue) > v_m )
      return 1;
    else
      return 0;
}

int TrajectoryGeneratorPoly::checkAccLimit( double root, double a_m, VectorXd poly )
{    
    double lvalue = 0.0;
    double t_ex = root;
    
    for(int i = 0; i < poly.size(); i++)
      lvalue += poly(i) * pow( t_ex, i );
    
    /*cout<<"lvalue: "<<lvalue<<endl; 
    cout<<"_scale: "<<_scale<<endl;*/

    lvalue /= _scale;
    
    /*cout<<"_scale: "<<_scale<<endl;    
    cout<<"acceleration poly: \n"<<poly<<endl;
    cout<<"extreme acceleration value : "<<lvalue<<" , occurs at time : "<<t_ex<<endl;*/
    
    if( fabs(lvalue) > a_m )
      return 1;
    else
      return 0;
}   
    
void TrajectoryGeneratorPoly::addVelExConstrain( SMatrixXd & Alin_v, vector<int> segment_no, vector<double> t_ex, int idx )
{
    assert(t_ex.size() == segment_no.size());
    int addConsNum = t_ex.size();
    int velConsNum = Alin_v.rows();

    //ROS_WARN("adding extra constraints on velocity extremum");
    SMatrixXd Alin_v_aug;
    Alin_v_aug.resize(velConsNum + addConsNum, TotalPoly_num);
    insertBlock(Alin_v_aug, 0,  0, Alin_v);

    for(int k = 0; k < addConsNum; k ++)
    {  
        SMatrixXd Alin_v_seg(1, TotalPoly_num);
        Alin_v_seg.setZero();
        int    seg_no  = segment_no[k];
        double time_ex = t_ex[k];

        //ROS_WARN("extreme occurs in No %d segment, in %d axis, at time %f", seg_no, idx, time_ex);
        for(int j = 0; j < D1Poly_num; j ++)
            if( j > 0)
                Alin_v_seg.insert( 0, idx * D1Poly_num + seg_no * Poly_num + j ) = j * pow( time_ex, j -1 );

        insertBlock(Alin_v_aug, velConsNum + k, 0, Alin_v_seg);
    }

    Alin_v = Alin_v_aug;
}

void TrajectoryGeneratorPoly::addAccExConstrain( SMatrixXd & Alin_a, vector<int> segment_no, vector<double> t_ex, int idx )
{
    assert(t_ex.size() == segment_no.size());
    int addConsNum = t_ex.size();
    int accConsNum = Alin_a.rows();

    SMatrixXd Alin_a_aug;
    Alin_a_aug.resize(accConsNum + addConsNum, TotalPoly_num);
    insertBlock(Alin_a_aug, 0,  0, Alin_a);

    for(int k = 0; k < addConsNum; k ++)
    {  
        SMatrixXd Alin_a_seg(1, TotalPoly_num);
        Alin_a_seg.setZero();
        int    seg_no  = segment_no[k];
        double time_ex = t_ex[k];

        for(int j = 0; j < D1Poly_num; j ++)
            if( j > 1)
                Alin_a_seg.insert( 0, idx * D1Poly_num + seg_no * Poly_num + j ) = j * (j - 1) * pow( time_ex, j -2 );

        insertBlock(Alin_a_aug, accConsNum + k, 0, Alin_a_seg);
    }

    Alin_a = Alin_a_aug;
}

void TrajectoryGeneratorPoly::addPosExConstrain( vector<SMatrixXd> & ConQua, vector<SMatrixXd> & ConLin, VectorXd & bin,
                                          vector<int> segment_no, vector<double> t_ex )
{
    // add extreme constraint at pos(t_ex) extreme where out of the ball path ( Quadratic constraint )
    // add dynamic constraint at vel(t_ex) and acc(t_ex) where exceeds the maximum vel and acc ( Linear constraint )
    assert(t_ex.size() == segment_no.size());
    int addConsNum = t_ex.size();

    VectorXd bin_tmp = bin;
    bin.resize( bin_tmp.size() + t_ex.size(), 1 );
    
    VectorXd bin_add(t_ex.size());

    for(int i = 0; i < addConsNum; i ++ ){
        int idx = segment_no[i];
        bin_add( i ) =  _Radius(idx)  *  _Radius(idx)
                     -  _Path(idx, 0) *  _Path(idx, 0) 
                     -  _Path(idx, 1) *  _Path(idx, 1) 
                     -  _Path(idx, 2) *  _Path(idx, 2) - _delta_pos;
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

    SMatrixXd Q_qc_add_tmp = Q_qc_add;//invMapP_D.transpose() * Q_qc_add * invMapP_D;
      
    ConQua.push_back(Q_qc_add_tmp);

    for(int i = 0; i < D1Poly_num; i ++ )
      F_qc_add1D.insert(0, i) = pow( t_ex[k], i );

    insertBlock(F_qc_add, 0, idx * Poly_num, -2.0 * _Path( idx, 0 ) * F_qc_add1D ); 
    insertBlock(F_qc_add, 0, idx * Poly_num + D1Poly_num, -2.0 * _Path( idx, 1 ) * F_qc_add1D );  
    insertBlock(F_qc_add, 0, idx * Poly_num + 2 * D1Poly_num, -2.0 * _Path( idx, 2 ) * F_qc_add1D ); 

    SMatrixXd F_qc_add_tmp = F_qc_add;// * invMapP_D;
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

static pair<SMatrixXd, SMatrixXd> MapPtoD( VectorXd Time, int seg_num, const int POLYORDER )
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
