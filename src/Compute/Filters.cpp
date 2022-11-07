#include "Filters.h"

#include "graphics/Vector2f.h"
#include "graphics/Vector3f.h"

#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>

#include <QString>

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
    averages.resize(numDimensions, 0);
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

void writeDimensionRanking(const std::vector<std::vector<int>>& ranking, const std::vector<QString>& names)
{
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream fileName;
    fileName << "rankings";
    fileName << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    fileName << ".csv";

    std::ofstream myfile;
    myfile.open(fileName.str());
    for (int i = 0; i < ranking.size(); i++)
    {
        for (int d = 0; d < ranking[i].size(); d++)
        {
            if (d != 0) myfile << ',';
            myfile << names[ranking[i][d]].toStdString();
        }
        myfile << std::endl;
    }

    myfile.close();
    std::cout << "Rankings written to file" << std::endl;
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
            diffAverages[d] = averages[0][d] - averages[1][d];
        }

        // Sort averages from high to low
        dimRanking.resize(numDimensions);
        std::iota(dimRanking.begin(), dimRanking.end(), 0);

        std::stable_sort(dimRanking.begin(), dimRanking.end(), [&diffAverages](size_t i1, size_t i2) {return diffAverages[i1] > diffAverages[i2]; });
    }

    void radiusPeakFilterHD(int seedPoint, const DataMatrix& dataMatrix, std::vector<std::vector<int>> floodPoints, std::vector<int>& dimRanking)
    {
        int numDimensions = dataMatrix.cols();

        std::vector<int> nearIndices;
        std::vector<int> farIndices;

        for (int wave = 0; wave < floodPoints.size()/2; wave++)
        {
            nearIndices.insert(nearIndices.end(), floodPoints[wave].begin(), floodPoints[wave].end());
        }
        for (int wave = floodPoints.size() / 2; wave < floodPoints.size(); wave++)
        {
            farIndices.insert(farIndices.end(), floodPoints[wave].begin(), floodPoints[wave].end());
        }

        std::vector<float> nearAverages;
        std::vector<float> farAverages;

        computeDimensionAverage(dataMatrix, nearIndices, nearAverages);
        computeDimensionAverage(dataMatrix, farIndices, farAverages);

        std::vector<float> diffAverages(numDimensions);
        for (int d = 0; d < numDimensions; d++)
            diffAverages[d] = nearAverages[d] - farAverages[d];

        // Sort averages from high to low
        dimRanking.resize(numDimensions);
        std::iota(dimRanking.begin(), dimRanking.end(), 0);

        std::stable_sort(dimRanking.begin(), dimRanking.end(), [&diffAverages](size_t i1, size_t i2) {return diffAverages[i1] > diffAverages[i2]; });
    }
}
