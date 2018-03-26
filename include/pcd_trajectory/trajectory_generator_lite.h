#ifndef _TRAJECTORY_GENERATOR_BEZIER_H_
#define _TRAJECTORY_GENERATOR_BEZIER_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SMatrixXd;

class TrajectoryGeneratorBezier {
private:

public:
        MatrixXd _Path;
        VectorXd _Radius;
        VectorXd _Time;
        VectorXd _Scale;

        double initScale, lstScale;

        TrajectoryGeneratorBezier();

        ~TrajectoryGeneratorBezier();

        /* Use Bezier curve for the trajectory */
        MatrixXd BezierPloyCoeffGeneration(
            const MatrixXd &Path,
            const VectorXd &Radius,
            const VectorXd &Time,
            const vector<int>     &multi_poly_order,
            const vector<MatrixXd> &MQMList,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int minimize_order );  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
};

#endif
