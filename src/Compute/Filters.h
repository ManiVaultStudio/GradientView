#pragma once

#include <Eigen/Eigen>

#include <vector>

using DataMatrix = Eigen::MatrixXf;

namespace filters
{
    void spatialCircleFilter(int seedPoint, float projSize, const DataMatrix& dataMatrix, const DataMatrix& projMatrix, std::vector<int>& dimRanking);

    void radiusPeakFilterHD(int seedPoint, const DataMatrix& dataMatrix, std::vector<std::vector<int>> floodPoints, std::vector<int>& dimRanking);
}
