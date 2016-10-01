#pragma once
#include "Eigen/dense"
#include "Eigen/SVD"
#include "Utility.hpp"

class Procruste
{   
    public:
    // X*T ≈ Y
    // find the best T that transform X to fit Y
    // return the Frobenius norm between X and Y

    template<typename MatrixType>
    static double procruste( MatrixType & X
                            , MatrixType & Y
                            ,Eigen::Matrix3d & T);


    static double procruste(std::vector<POINT>& X, std::vector<POINT>& Y);


    // convert POINT vector into matrix
    static void Point_to_Matrix(std::vector<POINT> & vp
            , Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> > & mat);

    // convert PDB data int matrix
    static void PDB_to_Matrix(std::vector<std::string> & data,Eigen::Matrix<double,Eigen::Dynamic,3> & mat);

    static void test();


};
