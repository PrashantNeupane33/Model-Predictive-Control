#include <vector>
#include <Eigen/Dense>

struct Coordinate {
    double x = 0.0;
    double y = 0.0;
};

struct SplineSegment {
    Eigen::VectorXd cx;
    Eigen::VectorXd cy;
    double t0;
    double t1;
};

class PathGenerator {
public:
    PathGenerator(int polynomialOrder, double resolution, std::vector<Coordinate> controlPoints);

    std::vector<Coordinate> getPath();

private:
    int polynomialOrder;
    double resolution;
    std::vector<Coordinate>     controlPoints;
    std::vector<SplineSegment>  segments;
    std::vector<Coordinate>     path;

    void addControlPoint(double x, double y);
    void generateSpline();
    void discretize();

    SplineSegment fitSegment(const Coordinate& p0,
                             const Coordinate& p1,
                             const Coordinate& p2,
                             double t0, double t1);
};
