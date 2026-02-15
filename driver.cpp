#include <iostream>
#include<Eigen/Dense>
#include "Controller.hpp"

using namespace Eigen;
using Eigen::placeholders::all;

int main()
{

	unsigned int f=20;// prediction horizon
	unsigned int v=18;// control horizon

	double m1=2  ; double m2=2   ; double k1=100  ; double k2=200 ; double d1=1  ; double d2=5; 

    Matrix <double,4,4> Ac {{0, 1, 0, 0},
                            {-(k1+k2)/m1 ,  -(d1+d2)/m1 , k2/m1 , d2/m1},
                            {0 , 0 ,  0 , 1},
                            {k2/m2,  d2/m2, -k2/m2, -d2/m2}};
    Matrix <double,4,1> Bc {{0},{0},{0},{1/m2}};
    Matrix <double,1,4> Cc {{1,0,0,0}};

    Matrix <double,4,1> x0 {{0},{0},{0},{0}}; 

	unsigned int n = Ac.rows();  unsigned int m = Bc.cols(); unsigned int r = Cc.rows();

	//# discretization constant
	double sampling=0.05;

	// # model discretization
	// identity matrix
	MatrixXd In;
	In= MatrixXd::Identity(n,n);
	MatrixXd A;
	MatrixXd B;
	MatrixXd C;
	A.resize(n,n);
	B.resize(n,m);
	C.resize(r,n);
	A=(In-sampling*Ac).inverse();
	B=A*sampling*Bc;
	C=Cc;

	// Weights
	double Q0=0.0000000011;
	double Qother=0.0001;
	double predWeight=10;
	//Tuples
	auto systemMatrices = std::make_tuple(A,B,C);
	auto horizons = std::make_tuple(v,f);
	auto weights = std::make_tuple(Q0, Qother, predWeight);

	unsigned int timeSteps=300;

	//# pulse trajectory

	MatrixXd desiredTrajectory;
	desiredTrajectory.resize(timeSteps,1);
	desiredTrajectory.setZero();

	MatrixXd tmp1;
	tmp1=MatrixXd::Ones(100,1);

	desiredTrajectory(seq(0,100-1),all)=tmp1;
	desiredTrajectory(seq(200,timeSteps-1),all)=tmp1;

	MPC  mpc(systemMatrices, horizons,weights,x0,desiredTrajectory);

	// Main control loop
	for (int index1=0; index1<timeSteps-f-1; index1++)
	{
	  mpc.computeControlInputs();    
	}
    
	mpc.saveData("data/trajectory.csv", "data/computedInputs.csv", 
								"data/states.csv", "data/outputs.csv","data/Omatrix.csv","data/Mmatrix.csv");

	std::cout<<"Simulation completed!"<<std::endl;
	return 0;

}

