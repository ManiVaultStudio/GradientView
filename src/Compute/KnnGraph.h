#pragma once

#include <Eigen/Eigen>
#include <knncpp.h>

typedef knncpp::Matrixi Matrixi;
using DataMatrix = Eigen::MatrixXf;

using KdTree = knncpp::KDTreeMinkowskiX<float, knncpp::ManhattenDistance<float>>;

class KnnGraph
{
public:
    void build(const DataMatrix& data, KdTree* kdTree, int k);

    std::vector<std::vector<int>> neighbours;
    int numNeighbours;
};
