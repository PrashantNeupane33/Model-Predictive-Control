#include<iostream>
#include<string>
#include<fstream>
#include <tuple>
#include<Eigen/Dense>
#include "Controller.hpp"

using namespace Eigen;
using namespace std;
using Eigen::placeholders::all;


MPC::MPC(std::tuple <MatrixXd, MatrixXd, MatrixXd> systemMatrices,
		 std::tuple<unsigned int, unsigned int> horizons,
		 std::tuple<double, double, double> weights,
		 MatrixXd _initialState, MatrixXd _desiredTrajectory):

	A(get<0>(systemMatrices)), B(get<1>(systemMatrices)), C(get<2>(systemMatrices)),
    f(get<1>(horizons)), v(get<0>(horizons)),
    x0(_initialState), 
    desiredInput(_desiredTrajectory),
    k(0)
{
    n = A.rows();   m = B.cols(); r = C.rows();
    unsigned int maxSimulationSamples = desiredInput.rows() - f;
    
	// States
    states.resize(n, maxSimulationSamples); 
    states.setZero();
    states.col(0)=x0;

	// Inputs
    inputs.resize(m, maxSimulationSamples-1); 
    inputs.setZero();
   
	// Outpus
    outputs.resize(r, maxSimulationSamples-1);
    outputs.setZero();

    setObservabilityMatrix();
	setToeplitzMatrix();
	setGainMatrix(weights);
}

void MPC::setObservabilityMatrix(){

    O.resize(f*r,n);
    O.setZero();

    MatrixXd _temp;
    _temp.resize(n,n);
    
    for (int i=0; i<f;i++)
    {
        if (i == 0)
            _temp=A;
        else
            _temp = _temp*A;

        O(seq(i*r,(i+1)*r-1),all)=C*_temp;
    }
}

void MPC::setToeplitzMatrix(){

    M.resize(f*r,v*m);
    M.setZero();

    MatrixXd _temp;
    MatrixXd sumLast;
    sumLast.resize(n,n);
    sumLast.setZero();

    for (int i=0; i<f;i++)
    {
        // until the control horizon
        if(i<v)
        {
            for (int j=0; j<i+1;j++)
            {
                if (j==0)
                    _temp=MatrixXd::Identity(n,n);
                else
                    _temp=_temp*A;

                M(seq(i*r,(i+1)*r-1),seq((i-j)*m,(i-j+1)*m-1))=C*_temp*B;
            }
        }
        // from the control horizon until the prediction horizon
        else
        {
            for(int j=0;j<v;j++)
            {
                if (j==0)
                {
                        sumLast.setZero();
                        for (int s=0;s<i+v+2;s++)
                        {
                            if (s == 0)
                                _temp = MatrixXd::Identity(n,n);
                            else
                                _temp=_temp*A;

                            sumLast=sumLast+_temp;
                        }
                        M(seq(i*r,(i+1)*r-1),seq((v-1)*m,(v)*m-1))=C*sumLast*B;
				}
                else
                {
                    _temp=_temp*A;
                    M(seq(i*r,(i+1)*r-1),seq((v-1-j)*m,(v-j)*m-1))=C*_temp*B;
                }
            }
        }
	}
}

std::tuple<MatrixXd,MatrixXd> MPC::getWeightMatrices(tuple<double, double, double> weights)
{

	MatrixXd wt1;
	wt1.resize(v*m,v*m);
	wt1.setZero();

	for (int i=0; i<v;i++)
	{
	  if (i==0)
		 wt1(seq(i*m,(i+1)*m-1),seq(i*m,(i+1)*m-1)) = MatrixXd::Identity(m,m);
	  else
	  {
		 wt1(seq(i*m,(i+1)*m-1),seq(i*m,(i+1)*m-1)) = MatrixXd::Identity(m,m);
		 wt1(seq(i*m,(i+1)*m-1),seq((i-1)*m,(i)*m-1)) = -MatrixXd::Identity(m,m);
	  }
	}
	MatrixXd wt2;
	wt2.resize(v*m,v*m);
	wt2.setZero();
	wt2(0,0) = get<0>(weights);
	for (int i=1; i<v; i++)
		wt2(i*m,i*m)=get<1>(weights);

	MatrixXd wt3;
	wt3 = (wt1.transpose())*wt2*wt1;

	MatrixXd wt4;
	wt4.resize(f*r,f*r);
	wt4.setZero();
	
	for(int i =0;i<f;i++)
		wt4(i*r,i*r) = get<2>(weights);
	return std::make_tuple(wt3,wt4);
}

void MPC::setGainMatrix(std::tuple<double, double, double> weights){
	 auto weightMatrices = MPC::getWeightMatrices(weights);
	 MatrixXd wt1 = get<0>(weightMatrices);
	 MatrixXd wt2 = get<1>(weightMatrices);

     MatrixXd _tmp;
     _tmp.resize(v*m,v*m);
     _tmp=M.transpose()*wt2*M+wt1;
     gainMatrix=(_tmp.inverse())*(M.transpose())*wt2;
}

// Propagets the dynamics and computes the control input
void MPC::computeControlInputs()
{
    MatrixXd desiredControlTrajectory;
    desiredControlTrajectory = desiredInput(seq(k,k+f-1),all);

    MatrixXd trackingError;
    trackingError=desiredControlTrajectory-O*states.col(k);

    //# compute the control sequence
    MatrixXd inputSequence;
    inputSequence = gainMatrix*trackingError;

    inputs.col(k)=inputSequence(seq(0,m-1),all); // apply only the first term as input

    // propagate the dynamics
    states.col(k+1)=A*states.col(k)+B*inputs.col(k);
    outputs.col(k)=C*states.col(k);
    k=k+1;
}

// void MPC::computeControlInputs()
// {
//     MatrixXd desiredControlTrajectory;
//     desiredControlTrajectory = desiredInput(seq(k,k+f-1),all);  // f x r
//
//     // Reshape to column vector using Map
//     Map<MatrixXd> desiredVector(desiredControlTrajectory.data(), f*r, 1);
//
//     MatrixXd trackingError;
//     trackingError = desiredVector - O*states.col(k);
//
//     MatrixXd inputSequence;
//     inputSequence = gainMatrix*trackingError;
//
//     inputs.col(k) = inputSequence(seq(0,m-1),all);
//
//     states.col(k+1) = A*states.col(k) + B*inputs.col(k);
//     outputs.col(k) = C*states.col(k);
//     k = k+1;
// }

void MPC::writeToCSV(const std::string& filename, const MatrixXd& matrix) const {
	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
        ofstream file(filename);
        if (file.is_open()) 
            file << matrix.format(CSVFormat);
         else 
            std::cerr << "Error: Could not open " << filename << std::endl;
}

void MPC::saveData(string desiredInput_File, string input_File, 
							string state_File, string output_File, string O_File, string M_File) const
{
		writeToCSV(desiredInput_File, desiredInput);
        writeToCSV(input_File, inputs);
        writeToCSV(state_File, states);
        writeToCSV(output_File, outputs);
        writeToCSV(O_File, O);
        writeToCSV(M_File, M);
}
