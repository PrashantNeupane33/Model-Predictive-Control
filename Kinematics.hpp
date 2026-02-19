#ifndef KINEMATICS_H
#define KINEMATICS_H

#include<Eigen/Dense>

struct ControlInput {
    double vx_cmd    = 0.0;
    double vy_cmd    = 0.0;
    double omega_cmd = 0.0;
};

struct WheelSpeeds {
    double w1 = 0.0;
    double w2 = 0.0;
    double w3 = 0.0;
    double w4 = 0.0;
};

struct RobotState {
    double x     = 0.0;
    double y     = 0.0;
    double theta = 0.0;
    double vx    = 0.0;
    double vy    = 0.0;
    double omega = 0.0;
};

class OmniKinematics{
	private:
    double wheelRadius;
    double robotRadius;
	Eigen::Matrix<double, 4, 3> Jacobian;
    Eigen::Matrix<double, 3, 4> JacobianInv;

    void buildJacobian();
	void updateState();

	public:
		OmniKinematics(double wheelRadius, double robotRadius);
		WheelSpeeds  inverseKinematics(const ControlInput& cmd) const;
    	ControlInput forwardKinematics(const WheelSpeeds& ws)   const;

};
#endif 
