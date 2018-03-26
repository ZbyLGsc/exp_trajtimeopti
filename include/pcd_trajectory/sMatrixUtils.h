#include "pcd_trajectory/trajectory_generator.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

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