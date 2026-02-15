#ifndef MODELPREDICTIVECONTROLLER_H
#define MODELPREDICTIVECONTROLLER_H

#include<string>
#include <tuple>
#include<Eigen/Dense>
using namespace Eigen;

class MPC{

    private:
        unsigned int k; // current time step
        unsigned int m,n,r; // (m,n,r) -> (input, state, output) dimension

	    MatrixXd A,B,C,Q; // system matrices
        MatrixXd W3,W4;   // weighting matrices
	    MatrixXd x0;      // initial state
        MatrixXd desiredInput; // total desired trajectory
        unsigned int f,v; // f- prediction horizon, v - control horizon

        MatrixXd states;
        MatrixXd inputs;
        MatrixXd outputs;

		// Lifted state matrix. O -> observability matrix, M -> Toeplitz matrix
        MatrixXd O;
        MatrixXd M;

        // control gain matrix
        MatrixXd gainMatrix;

		void setObservabilityMatrix();
		void setToeplitzMatrix();
		void setGainMatrix(std::tuple<double, double, double> weights);
		std::tuple<MatrixXd,MatrixXd> getWeightMatrices(std::tuple<double, double, double> weights);
		void writeToCSV(const std::string& filename, const MatrixXd& matrix) const;

    public:
		MPC(std::tuple <MatrixXd, MatrixXd, MatrixXd> systemMatrices,
								  std::tuple<unsigned int, unsigned int> horizons,
								  std::tuple<double, double, double> weights,
								  MatrixXd _initialState, MatrixXd _desiredTrajectory);
        
        // function to propagate the dynamics and to compute the solution of the MPC problem
        void computeControlInputs();

        // function to save data
        void saveData(std::string _desiredInput, std::string inputFile, 
							std::string stateFile, std::string outputFile,std::string OFile, std::string MFile) const;
};
#endif
