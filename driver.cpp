#include <iostream>
#include <Eigen/Dense>
#include "Controller.hpp"
#include <cmath>

using namespace Eigen;
using Eigen::placeholders::all;

int main()
{
    unsigned int f = 40;   // prediction horizon
    unsigned int v = 35;   // control horizon

    double sampling = 0.05;

    MatrixXd C = MatrixXd::Identity(3, 3);

    VectorXd x0(3);
    x0 << 0.0, 0.0, 0.0;

    // Actuator limits
    VectorXd u_min(3), u_max(3);
    u_min << -3.0, -3.0, -5.0;
    u_max <<  3.0,  3.0,  5.0;

    // Weights
    double Q0       = 0.0;
    double Qother   = 0.1;  // control smoothness weight
    double predWeight = 2000.0; // tracking weight

    // Tuples
    auto horizons = std::make_tuple(v, f);
    auto weights  = std::make_tuple(Q0, Qother, predWeight);

    // Reference trajectory: circle
    double t_total   = 10.0;
    double R_circle  = 2.0;
    unsigned int timeSteps = static_cast<unsigned int>(t_total / sampling) + f + 10;

    MatrixXd desiredTrajectory;
    desiredTrajectory.resize(timeSteps, 3);

    double omega_ref = 2.0 * M_PI / t_total;

    for (int i = 0; i < timeSteps; i++)
    {
        double t = i * sampling;
        double X_ref     =  R_circle * cos(omega_ref * t - M_PI/2.0);
        double Y_ref     =  R_circle * sin(omega_ref * t - M_PI/2.0) + R_circle;
        double theta_ref =  atan2( cos(omega_ref * t - M_PI/2.0),
                                  -sin(omega_ref * t - M_PI/2.0));
        desiredTrajectory(i, 0) = X_ref;
        desiredTrajectory(i, 1) = Y_ref;
        desiredTrajectory(i, 2) = theta_ref;
    }

	for(int i = 0; i < 5; i++)
    std::cout << "ref[" << i << "] = " 
         << desiredTrajectory(i,0) << ", " 
         << desiredTrajectory(i,1) << std::endl;

    // Construct MPC
    MPC mpc(C, horizons, weights, x0, desiredTrajectory, sampling, u_min, u_max);

    // Main control loop
    for (int i = 0; i < timeSteps - f - 1; i++)
    {
        mpc.computeControlInputs();
    }

    mpc.saveData("data/trajectory.csv",     "data/computedInputs.csv",
                 "data/states.csv",         "data/outputs.csv",
                 "data/Omatrix.csv",        "data/Mmatrix.csv");

    std::cout << "Simulation completed!" << std::endl;
    return 0;
}
