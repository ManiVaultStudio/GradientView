#pragma once

#include <Eigen/Eigen>

#include <vector>

using DataMatrix = Eigen::MatrixXf;

namespace filters
{
    void spatialCircleFilter(int seedPoint, float projSize, const DataMatrix& dataMatrix, const DataMatrix& projMatrix, std::vector<int>& dimRanking);
}
