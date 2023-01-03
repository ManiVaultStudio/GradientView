#pragma once

#include <Eigen/Eigen>
#include <knncpp.h>

#pragma warning(push, 0)
#define ANNOYLIB_MULTITHREADED_BUILD
#include <annoylib.h>
#include <kissrandom.h>
#pragma warning(pop)

using AnnoyIndex = Annoy::AnnoyIndex<int, float, Annoy::Angular, Annoy::Kiss32Random, Annoy::AnnoyIndexMultiThreadedBuildPolicy>;

typedef knncpp::Matrixi Matrixi;
using DataMatrix = Eigen::MatrixXf;

using KdTree = knncpp::KDTreeMinkowskiX<float, knncpp::ManhattenDistance<float>>;

class KnnGraph
{
public:
    void build(const DataMatrix& data, KdTree* kdTree, int k);
    void build(const DataMatrix& data, AnnoyIndex* index, int k);

    std::vector<std::vector<int>> neighbours;
    int numNeighbours;
};
