#ifndef _TRAJECTORY_QP_GENERATOR_H_
#define _TRAJECTORY_QP_GENERATOR_H_

#include <Eigen/Eigen>
#include <vector>
#include "pcd_trajectory/gridNode.h"

class TrajectoryGenerator {
private:
        Eigen::VectorXd _Dx;
        Eigen::VectorXd _Dy;
        Eigen::VectorXd _Dz;
public:
        TrajectoryGenerator();

        ~TrajectoryGenerator();

        Eigen::MatrixXd PolyQPGeneration( const Eigen::MatrixXd &Path, const State &init_s, double time, double &qp_cost); 

        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > getInitialD(); // Initial Derivatives variable for the following optimization 
};

#endif
