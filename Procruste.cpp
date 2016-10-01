#include "Procruste.hpp"
#include<iostream>
using namespace std;

template<typename MatrixType>
double Procruste::procruste(MatrixType & X
                            , MatrixType & Y
                            , Eigen::Matrix3d & T)
{
    if(X.rows() != Y.rows())
        return -1;

    Eigen::MatrixXd temp,A;
    Eigen::Vector3d v;
    double traceY,traceX,traceSource;

    //centralize X and Y    
    v=X.colwise().mean();
    X.rowwise()-=v.transpose();
    v=Y.colwise().mean();
    Y.rowwise()-=v.transpose();


    traceX = (X*X.transpose()).trace();
    traceY = (Y*Y.transpose()).trace();

    A = X.transpose()*Y;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // we want a porper rotation without reflection
    T = svd.matrixV() * svd.matrixU().transpose();
    
    if(T.determinant() < 0)
    {
        cout<<"correct rotation\n";
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        I(2,2) = -1;
        T = svd.matrixV() * I * svd.matrixU().transpose();
    }

    //return 1+traceX/traceY-2*svd.singularValues().sum()/traceY;

    return sqrt((traceX + traceY - 2*svd.singularValues().sum())/X.rows());

}



double Procruste::procruste(vector<point> & X, vector<POINT> & Y)
{
    //convert point vector to matrix
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,3, Eigen::RowMajor> > mX(NULL,0,3);
    Point_to_Matrix(X,mX);
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,3, Eigen::RowMajor> > mY(NULL,0,3);
    Point_to_Matrix(Y,mY);
    
    Eigen::Matrix3d T;
    return procruste(mX,mY,T);
}



void Procruste::Point_to_Matrix(std::vector<POINT> & vp
                                , Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> > & mat)
{
    new (&mat) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> >((double*)vp.data(), vp.size(),3);
}


void Procruste:: test()
{
    vector<POINT> vp;
    POINT pt;
    for(int i =0; i < 10; i++)
    {
        pt.x = i;
        pt.y = i*i;
        pt.z = i*i*i;
        vp.push_back(pt);
    }

    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,3, Eigen::RowMajor> > m(NULL,0,3);
    Point_to_Matrix(vp,m);

    vector<POINT> vp2;
    for(int i =0; i < 10; i++)
    {
        pt.x = i;
        pt.y = i*i;
        pt.z = i*i*i;
        vp2.push_back(pt);
    }

    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,3, Eigen::RowMajor> > m2(NULL,0,3);
    Point_to_Matrix(vp2,m2);


    Eigen::Matrix3d m3,T;
    T << 0.5000  , -0.8660 ,        0
           , 0.8660  ,  0.5000   ,      0
            ,         0   ,      0 ,   1.0000;

    Eigen::Vector3d v;
    v=m.colwise().mean();
    m.rowwise()-=v.transpose();
    v=m2.colwise().mean();
    m2.rowwise()-=v.transpose();


    
    m2 = m2*T;
    cout<<m<<'\n'<<endl;
    cout<<m2<<'\n'<<endl;


    cout<<procruste(m,m2,m3)<<endl;
    cout<<m3<<'\n'<<endl;
    cout<<m3.determinant()<<'\n'<<endl;

    cout<<m*m3<<'\n'<<endl;
    cout<<m*m3.transpose()<<'\n'<<endl;



    cout<<procruste(vp,vp2)<<endl;
}
