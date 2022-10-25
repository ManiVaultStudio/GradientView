#include "Filters.h"

#include "graphics/Vector2f.h"
#include "graphics/Vector3f.h"

#include <numeric>

using namespace hdps;

void findPointsInRadius(Vector2f center, float radius, const DataMatrix& projMatrix, std::vector<int>& indices)
{
    float radiusSqr = radius * radius;
    for (int i = 0; i < projMatrix.rows(); i++)
    {
        Vector2f pos(projMatrix(i, 0), projMatrix(i, 1));

        Vector2f diff = center - pos;
        float sqrLen = diff.x * diff.x + diff.y * diff.y;

        if (sqrLen < radiusSqr)
        {
            indices.push_back(i);
        }
    }
}

void computeDimensionAverage(const DataMatrix& data, const std::vector<int>& indices, std::vector<float>& averages)
{
    int numDimensions = data.cols();
    averages.resize(numDimensions);
    for (int d = 0; d < numDimensions; d++)
    {
        for (const int& index : indices)
        {
            float v = data(index, d);
            averages[d] += v;
        }
        averages[d] /= indices.size();
    }
}

namespace filters
{
    void spatialCircleFilter(int seedPoint, float projSize, const DataMatrix& dataMatrix, const DataMatrix& projMatrix, std::vector<int>& dimRanking)
    {
        int numDimensions = dataMatrix.cols();

        // Small and large circle averages
        std::vector<std::vector<float>> averages(2);
        std::vector<std::vector<int>> circleIndices(2);

        Vector2f center = Vector2f(projMatrix(seedPoint, 0), projMatrix(seedPoint, 1));

        findPointsInRadius(center, 0.05f * projSize, projMatrix, circleIndices[0]);
        computeDimensionAverage(dataMatrix, circleIndices[0], averages[0]);
        findPointsInRadius(center, 0.1f * projSize, projMatrix, circleIndices[1]);
        computeDimensionAverage(dataMatrix, circleIndices[1], averages[1]);

        std::vector<float> diffAverages(numDimensions);
        for (int d = 0; d < numDimensions; d++)
        {
            diffAverages[d] = fabs(averages[0][d] - averages[1][d]);
        }

        // Sort averages from high to low
        dimRanking.resize(numDimensions);
        std::iota(dimRanking.begin(), dimRanking.end(), 0);

        std::stable_sort(dimRanking.begin(), dimRanking.end(), [&diffAverages](size_t i1, size_t i2) {return diffAverages[i1] > diffAverages[i2]; });
    }
}
