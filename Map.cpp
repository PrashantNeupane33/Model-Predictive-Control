#include "Map.hpp"
#include <stdexcept>
#include <cmath>

PathGenerator::PathGenerator(int polynomialOrder, double resolution, std::vector<Coordinate> controlPoints)
    : polynomialOrder(polynomialOrder), resolution(resolution), controlPoints(controlPoints)
{}

void PathGenerator::addControlPoint(double x, double y)
{
    controlPoints.push_back({x, y});
}

SplineSegment PathGenerator::fitSegment(const Coordinate& p0,
                                         const Coordinate& p1,
                                         const Coordinate& p2,
                                         double t0, double t1)
{
    int n = polynomialOrder + 1;
    double dt = t1 - t0;

    Eigen::MatrixXd A(n, n);
    Eigen::VectorXd bx(n), by(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A(i, j) = std::pow(t0 + i * dt / (n - 1), j);
        }
        bx(i) = p0.x + (p1.x - p0.x) * i / (n - 1);
        by(i) = p0.y + (p1.y - p0.y) * i / (n - 1);
    }

    SplineSegment seg;
    seg.cx = A.colPivHouseholderQr().solve(bx);
    seg.cy = A.colPivHouseholderQr().solve(by);
    seg.t0 = t0;
    seg.t1 = t1;

    return seg;
}

void PathGenerator::generateSpline()
{
    if (controlPoints.size() < 2)
        throw std::runtime_error("Need at least 2 control points");

    segments.clear();

    int n = controlPoints.size();

    for (int i = 0; i < n - 1; i++) {
        double t0 = static_cast<double>(i);
        double t1 = static_cast<double>(i + 1);

        const Coordinate& p0 = controlPoints[i];
        const Coordinate& p1 = controlPoints[i + 1];
        const Coordinate& p2 = (i + 2 < n) ? controlPoints[i + 2] : controlPoints[i + 1];

        segments.push_back(fitSegment(p0, p1, p2, t0, t1));
    }
}

void PathGenerator::discretize()
{
    path.clear();

    for (const auto& seg : segments) {
        int n = polynomialOrder + 1;
        double t = seg.t0;

        while (t <= seg.t1) {
            Eigen::VectorXd tvec(n);
            for (int j = 0; j < n; j++)
                tvec(j) = std::pow(t, j);

            Coordinate pt;
            pt.x = seg.cx.dot(tvec);
            pt.y = seg.cy.dot(tvec);
            path.push_back(pt);

            t += resolution;
        }
    }
}

std::vector<Coordinate> PathGenerator::getPath()
{
    generateSpline();
    discretize();
    return path;
}
