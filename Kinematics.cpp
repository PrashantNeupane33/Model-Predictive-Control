#include "Kinematics.hpp"
#include<Eigen/Dense>

OmniKinematics::OmniKinematics(double wheelRadius, double robotRadius)
    : wheelRadius(wheelRadius), robotRadius(robotRadius)
{
    buildJacobian();
}

// Initalize the Jacobian of four wheel omni drive
void OmniKinematics::buildJacobian()
{
    double R = robotRadius;
    double r = wheelRadius;

    Jacobian << -1,  1,  R,
           1,  1,  R,
          -1, -1,  R,
           1, -1,  R;

    Jacobian/= r;

    JacobianInv = Jacobian.completeOrthogonalDecomposition().pseudoInverse();
}

void OmniKinematics::updateState(RobotState currentState){
	// updateState here
}

WheelSpeeds OmniKinematics::inverseKinematics(const ControlInput& cmd) const
{
    Eigen::Vector3d twist(cmd.vx_cmd, cmd.vy_cmd, cmd.omega_cmd);
    Eigen::Vector4d speeds = Jacobian * twist;

    WheelSpeeds ws;
    ws.w1 = speeds(0);
    ws.w2 = speeds(1);
    ws.w3 = speeds(2);
    ws.w4 = speeds(3);

    return ws;
}

ControlInput OmniKinematics::forwardKinematics(const WheelSpeeds& ws) const
{
    Eigen::Vector4d speeds(ws.w1, ws.w2, ws.w3, ws.w4);
    Eigen::Vector3d twist = JacobianInv * speeds;

    ControlInput cmd;
    cmd.vx_cmd    = twist(0);
    cmd.vy_cmd    = twist(1);
    cmd.omega_cmd = twist(2);

    return cmd;
}
