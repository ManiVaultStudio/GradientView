#pragma once

#include <Eigen/Dense>
#include <vector>

using DataMatrix = Eigen::MatrixXf;

namespace compute
{
    float computeProjectionDiameter(const DataMatrix& projection, int xDim, int yDim);

    void findNeighbourhood(const DataMatrix& projection, int centerId, float radius, std::vector<int>& neighbourhood, int xDim, int yDim);

    void computeSpatialLocalDimensionality(DataMatrix& dataMatrix, DataMatrix& projMatrix, std::vector<float>& colors);
}
