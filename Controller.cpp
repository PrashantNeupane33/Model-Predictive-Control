#include<iostream>
#include<string>
#include<fstream>
#include<tuple>
#include<cmath>
#include<Eigen/Dense>
#include "Controller.hpp"

using namespace Eigen;
using namespace std;
using Eigen::placeholders::all;

MPC::MPC(MatrixXd _C,
         std::tuple<unsigned int, unsigned int> horizons,
         std::tuple<double, double, double> weights,
         VectorXd _initialState,
         MatrixXd _desiredTrajectory,
         double _sampling,
         VectorXd _u_min,
         VectorXd _u_max):
    C(_C),
    f(get<1>(horizons)), v(get<0>(horizons)),
    x0(_initialState),
    desiredInput(_desiredTrajectory),
    sampling(_sampling),
    u_min(_u_min),
    u_max(_u_max),
    k(0)
{
    n = 3;
    m = _u_min.size();
    r = C.rows();

    unsigned int maxSimulationSamples = desiredInput.rows() - f;

    states.resize(n, maxSimulationSamples);
    states.setZero();
    states.col(0) = x0;

    inputs.resize(m, maxSimulationSamples - 1);
    inputs.setZero();

    outputs.resize(r, maxSimulationSamples - 1);
    outputs.setZero();

    u_prev = VectorXd::Zero(m);

    auto wm = getWeightMatrices(weights);
    W3 = get<0>(wm);
    W4 = get<1>(wm);
}

void MPC::setObservabilityMatrix(const MatrixXd& A_k)
{
    O.resize(f*r, n);
    O.setZero();

    MatrixXd _temp = MatrixXd::Identity(n, n);

    for (int i = 0; i < f; i++)
    {
        _temp = _temp * A_k;
        O(seq(i*r, (i+1)*r-1), all) = C * _temp;
    }
}

void MPC::setToeplitzMatrix(const MatrixXd& A_k, const MatrixXd& B_k)
{
    M.resize(f*r, v*m);
    M.setZero();

    MatrixXd _temp;

    for (int i = 0; i < f; i++)
    {
        if (i < v)
        {
            for (int j = 0; j < i+1; j++)
            {
                if (j == 0)
                    _temp = MatrixXd::Identity(n, n);
                else
                    _temp = _temp * A_k;

                M(seq(i*r,(i+1)*r-1), seq((i-j)*m,(i-j+1)*m-1)) = C * _temp * B_k;
            }
        }
        else
        {
            for (int j = 0; j < v; j++)
            {
                if (j == 0)
                {
                    MatrixXd sumLast = MatrixXd::Zero(n, n);
                    _temp = MatrixXd::Identity(n, n);
                    for (int s = 0; s <= i; s++)
                    {
                        sumLast = sumLast + _temp;
                        _temp = _temp * A_k;
                    }
                    M(seq(i*r,(i+1)*r-1), seq((v-1)*m, v*m-1)) = C * sumLast * B_k;
                }
                else
                {
                    _temp = _temp * A_k;
                    M(seq(i*r,(i+1)*r-1), seq((v-1-j)*m,(v-j)*m-1)) = C * _temp * B_k;
                }
            }
        }
    }
}

std::tuple<MatrixXd,MatrixXd> MPC::getWeightMatrices(tuple<double, double, double> weights)
{
    // wt1: finite difference matrix for control smoothness
    MatrixXd wt1 = MatrixXd::Zero(v*m, v*m);
    for (int i = 0; i < v; i++)
    {
        wt1(seq(i*m,(i+1)*m-1), seq(i*m,(i+1)*m-1)) = MatrixXd::Identity(m, m);
        if (i > 0)
            wt1(seq(i*m,(i+1)*m-1), seq((i-1)*m,i*m-1)) = -MatrixXd::Identity(m, m);
    }

    // wt2: control rate weight — full m x m blocks on diagonal
    MatrixXd wt2 = MatrixXd::Zero(v*m, v*m);
    wt2(seq(0, m-1), seq(0, m-1)) = get<0>(weights) * MatrixXd::Identity(m, m);
    for (int i = 1; i < v; i++)
        wt2(seq(i*m,(i+1)*m-1), seq(i*m,(i+1)*m-1)) = get<1>(weights) * MatrixXd::Identity(m, m);

    MatrixXd W3 = wt1.transpose() * wt2 * wt1;

    // W4: tracking weight — full r x r blocks on diagonal
    MatrixXd W4 = MatrixXd::Zero(f*r, f*r);
    for (int i = 0; i < f; i++)
        W4(seq(i*r,(i+1)*r-1), seq(i*r,(i+1)*r-1)) = get<2>(weights) * MatrixXd::Identity(r, r);

    return std::make_tuple(W3, W4);
}

VectorXd MPC::nonlinearDynamics(const VectorXd& x, const VectorXd& u) const
{
    double theta = x(2);
    VectorXd x_dot(n);
    x_dot(0) = u(0)*cos(theta) - u(1)*sin(theta);
    x_dot(1) = u(0)*sin(theta) + u(1)*cos(theta);
    x_dot(2) = u(2);
    return x + sampling * x_dot;
}

std::tuple<MatrixXd,MatrixXd> MPC::linearizeModel(const VectorXd& x_bar, const VectorXd& u_bar) const
{
    double eps = 1e-5;
    MatrixXd A_lin(n, n), B_lin(n, m);
    VectorXd f0 = nonlinearDynamics(x_bar, u_bar);

    for (int i = 0; i < n; i++) {
        VectorXd x_pert = x_bar;
        x_pert(i) += eps;
        A_lin.col(i) = (nonlinearDynamics(x_pert, u_bar) - f0) / eps;
    }

    for (int i = 0; i < m; i++) {
        VectorXd u_pert = u_bar;
        u_pert(i) += eps;
        B_lin.col(i) = (nonlinearDynamics(x_bar, u_pert) - f0) / eps;
    }

    return {A_lin, B_lin};
}

void MPC::computeControlInputs()
{
    VectorXd x_k = states.col(k);
    // DO NOT wrap x_k(2) here — let theta grow unbounded
    // wrapping causes discontinuity in A_k*x_k vs f0

    auto [A_k, B_k] = linearizeModel(x_k, u_prev);

    setObservabilityMatrix(A_k);
    setToeplitzMatrix(A_k, B_k);

    MatrixXd refWindow = desiredInput(seq(k, k+f-1), all);

    VectorXd refVec(f * r);
    for (int i = 0; i < f; i++)
        refVec(seq(i*r, (i+1)*r-1)) = refWindow.row(i).transpose();

    VectorXd f0  = nonlinearDynamics(x_k, u_prev);
    VectorXd d_k = f0 - A_k*x_k - B_k*u_prev;

    VectorXd offset = VectorXd::Zero(f * r);
    VectorXd d_prop = d_k;
    for (int i = 0; i < f; i++) {
        offset(seq(i*r, (i+1)*r-1)) = C * d_prop;
        d_prop = A_k * d_prop;
    }

    VectorXd trackingError = refVec - O*x_k - offset;

    // Wrap only the theta component of tracking error
    for (int i = 0; i < f; i++)
        trackingError(i*r + 2) = atan2(sin(trackingError(i*r + 2)),
                                        cos(trackingError(i*r + 2)));

    MatrixXd _tmp = M.transpose()*W4*M + W3;
    VectorXd inputSequence = (_tmp.inverse() * M.transpose() * W4) * trackingError;

    if (inputSequence.hasNaN()) {
        cerr << "NaN at step " << k << " | det(_tmp)=" << _tmp.determinant() << endl;
        inputs.col(k)   = VectorXd::Zero(m);
        u_prev          = VectorXd::Zero(m);
        states.col(k+1) = nonlinearDynamics(x_k, VectorXd::Zero(m));
        // Wrap only for output
        VectorXd state_out = states.col(k);
        state_out(2) = atan2(sin(state_out(2)), cos(state_out(2)));
        outputs.col(k) = C * state_out;
        k++;
        return;
    }

    VectorXd u_apply = inputSequence(seq(0, m-1));
    u_apply = u_apply.cwiseMax(u_min).cwiseMin(u_max);

    inputs.col(k) = u_apply;
    u_prev        = u_apply;

    // Propagate — NO angle wrap on state, let theta grow continuously
    states.col(k+1) = nonlinearDynamics(x_k, u_apply);

    // Wrap only for output/plotting
    VectorXd state_out = states.col(k);
    state_out(2) = atan2(sin(state_out(2)), cos(state_out(2)));
    outputs.col(k) = C * state_out;

    k++;
}

void MPC::writeToCSV(const std::string& filename, const MatrixXd& matrix) const {
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    ofstream file(filename);
    if (file.is_open())
        file << matrix.format(CSVFormat);
    else
        cerr << "Error: Could not open " << filename << endl;
}

void MPC::saveData(string desiredInput_File, string input_File,
                   string state_File, string output_File,
                   string O_File, string M_File) const
{
    writeToCSV(desiredInput_File, desiredInput);
    writeToCSV(input_File,  inputs);
    writeToCSV(state_File,  states);
    writeToCSV(output_File, outputs);
    writeToCSV(O_File, O);
    writeToCSV(M_File, M);
}
