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

    class SpatialPeakFilter
    {
    public:
        SpatialPeakFilter();

        float getInnerFilterRadius() { return _innerFilterRadius; }
        float getOuterFilterRadius() { return _outerFilterRadius; }
        void setInnerFilterRadius(float size);
        void setOuterFilterRadius(float size);

        void computeDimensionRanking(int pointId, const DataMatrix& dataMatrix, const DataMatrix& projMatrix, float projSize, std::vector<int>& dimRanking);

    private:
        float _innerFilterRadius;
        float _outerFilterRadius;
    };

    class HDFloodPeakFilter
    {
    public:
        HDFloodPeakFilter();

        void setInnerFilterSize(int size);
        void setOuterFilterSize(int size);

        void computeDimensionRanking(int pointId, const DataMatrix& dataMatrix, const std::vector<std::vector<int>>& floodPoints, std::vector<int>& dimRanking);
    private:
        int _innerFilterSize;
        int _outerFilterSize;
    };
}
