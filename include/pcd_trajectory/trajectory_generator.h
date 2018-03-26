#ifndef _TRAJECTORY_GENERATOR_H_
#define _TRAJECTORY_GENERATOR_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SMatrixXd;
class TrajectoryGeneratorPoly {
private:
		vector<double> qp_cost;
public:
        MatrixXd _Path;
        VectorXd _Time;
        VectorXd _Radius;
        double _max_v, _max_a, _scale;
        double _delta_pos, _delta_vel, _delta_acc;
        double _d_eps;
        
        TrajectoryGeneratorPoly();

        ~TrajectoryGeneratorPoly();

        vector<MatrixXd> PloyCoeffGeneration(
            const MatrixXd &Path,
            const VectorXd &Radius,
            const VectorXd &Time,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int POLYORDER,
            double scale, double w_l, double w_s );

        MatrixXd solver( const int & Equ_con_num,  const int & Inequ_con_num, const int & Inequ_con_v_num, const int & Inequ_con_a_num,
                        const int & vel_ex_num, const int & acc_ex_num, const vector<SMatrixXd> & ConQua, const vector<SMatrixXd> & ConLin,
                        const SMatrixXd & Qobj, const SMatrixXd & Aeq, const SMatrixXd & Alin_v, const SMatrixXd & Alin_a, const VectorXd & beq, const VectorXd & bin,
                        const int & type );

        vector<double> getObjective();

        pair< vector<int>, vector<double> > checkPosEx( MatrixXd Poly, const int POLYORDER );
        pair< vector<int>, vector<double> > checkVelEx( MatrixXd Poly, const int POLYORDER );
        pair< vector<int>, vector<double> > checkAccEx( MatrixXd Poly, const int POLYORDER );

        void addPosExConstrain( vector< SMatrixXd > & ConQua, vector< SMatrixXd > & ConLin, VectorXd & bin, vector<int> segment_no, vector<double> t_ex );
        void addVelExConstrain( SMatrixXd & Alin_v, vector<int> segment_no, vector<double> t_ex, int idx );
        void addAccExConstrain( SMatrixXd & Alin_a, vector<int> segment_no, vector<double> t_ex, int idx );

        int checkPosLimit( double root, double radius, Vector3d center, VectorXd poly );
        int checkVelLimit( double root, double v_m, VectorXd poly ); 
        int checkAccLimit( double root, double a_m, VectorXd poly );
};

#endif
