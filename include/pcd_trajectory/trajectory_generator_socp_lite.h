#ifndef _TRAJECTORY_GENERATOR_SOCP_LITE_H_
#define _TRAJECTORY_GENERATOR_SOCP_LITE_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SMatrixXd;

class TrajectoryGeneratorSOCP {
private:
        double _objective;

public:
        MatrixXd _Path;
        VectorXd _Radius;
        VectorXd _Time;
        VectorXd _Scale;

        double initScale, lstScale;

        TrajectoryGeneratorSOCP();

        ~TrajectoryGeneratorSOCP();

        /* Use Bezier curve for the trajectory */
        Eigen::MatrixXd BezierPloyCoeffGenerationSOCP(
            const MatrixXd &Path,
            const VectorXd &Radius,
            const VectorXd &Time,
            const vector<int>     &multi_poly_order,
            const vector<MatrixXd> &FMList,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int minimize_order );  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  

        double getObjective()
        { 
            return _objective;
        };

};

#endif
