#ifndef _TRAJECTORY_GENERATOR_DISCRETE_H_
#define _TRAJECTORY_GENERATOR_DISCRETE_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SMatrixXd;
class TrajectoryGeneratorDiscrete {
private:
		double qp_cost;
public:
        vector<Vector3d> _Path;
        vector<double>   _Time;
        
        TrajectoryGeneratorDiscrete();
        ~TrajectoryGeneratorDiscrete();

        MatrixXd AccDiscreteGeneration(
            const vector<Vector3d> &Path,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double & max_v,
            const double & max_a,
            const double & h,
            const double & l);

        double getObjective() {return qp_cost;};
};

#endif