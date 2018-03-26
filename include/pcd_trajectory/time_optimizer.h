#ifndef _TIME_OPTIMIZER_H_
#define _TIME_OPTIMIZER_H_

#include <Eigen/Dense>
#include <vector>
#include "timeAllocator.h"

using namespace std;
using namespace Eigen;

class MinimumTimeOptimizer {
private:
        double _T;    // the reulsted minimum time for a feasible travesal
        MatrixXd _P;  // recording the polynomial's coefficients for further evaluation
        int _seg_num; // segment number of the piecewise trajectory
        int _poly_num1D;

        double _objective;
        Allocator * time_allocator; // for return the final result to high-level planer

        Vector3d getVelPoly(int k, double s);
        Vector3d getAccPoly(int k, double s);

        /*Vector3d getVelBezier(int k, double s);
        Vector3d getAccBezier(int k, double s);*/
        Vector3d getVel(int k, double s);
        Vector3d getAcc(int k, double s);
        // Vector3d getVelBezier(int, double);
        // Vector3d getAccBezier(int, double);

public:
        MinimumTimeOptimizer();
        ~MinimumTimeOptimizer();

        void MinimumTimeGeneration( const MatrixXd & polyCoeff, // the polynomial coefficients in virtual domain form 0.0 to 1.0
                                    const double & maxVel,      // maximum phisical limitation of v, a, and j
                                    const double & maxAcc,
                                    const double & maxJer,
                                    const int    & K,
                                    const double & w_a );   // discretize size of a grid in s from 0.0 to 1.0

        Allocator * GetTimeAllcoation() {return time_allocator;}
};

#endif
