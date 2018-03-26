#ifndef _TRAJECTORY_GENERATOR_H_
#define _TRAJECTORY_GENERATOR_H_

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
Eigen::MatrixXd BezierPloyCoeffGenerationSOCP(
            const Eigen::MatrixXd &Path,
            const Eigen::VectorXd &Radius,
            const Eigen::VectorXd &Time,
            const Eigen::MatrixXd &FM,
            const Eigen::MatrixXd &pos,
            const Eigen::MatrixXd &vel,
            const Eigen::MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,        // define the order of the polynomial functions
            const int minimize_order );  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  


Eigen::MatrixXd BezierPloyCoeffGenerationQCQP(
            const Eigen::MatrixXd &Path,
            const Eigen::VectorXd &Radius,
            const Eigen::VectorXd &Time,
            const Eigen::MatrixXd &MQM,
            const Eigen::MatrixXd &pos,
            const Eigen::MatrixXd &vel,
            const Eigen::MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,      
            const int minimize_order );
};

#endif
