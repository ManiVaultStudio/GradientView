#pragma once

#include <Eigen/Eigen>

#include <vector>
#include <QString>

using DataMatrix = Eigen::MatrixXf;

void writeDimensionRanking(const std::vector<std::vector<int>>& ranking, const std::vector<QString>& names);

namespace filters
{
    enum class FilterType
    {
        SPATIAL_PEAK,
        HD_PEAK
    };

    void spatialCircleFilter(int seedPoint, float projSize, const DataMatrix& dataMatrix, const DataMatrix& projMatrix, std::vector<int>& dimRanking);

    void radiusPeakFilterHD(int seedPoint, const DataMatrix& dataMatrix, std::vector<std::vector<int>> floodPoints, std::vector<int>& dimRanking);
}
