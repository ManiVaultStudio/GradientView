#pragma once

#include <cstdlib>
#include <stdint.h>

using idx_t = int64_t;

#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/gpu/GpuIndexFlat.h>
#include <faiss/gpu/StandardGpuResources.h>
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
using Brute = knncpp::BruteForce<float, knncpp::ManhattenDistance<float>>;

class KnnGraph
{
public:
    void build(const DataMatrix& data, Brute* brute, int k);
    void build(const DataMatrix& data, KdTree* kdTree, int k);
    //void build(const DataMatrix& data, flann::Index<flann::L2<float>>* index, int k);
    void build(const DataMatrix& data, faiss::IndexFlat* index, int k);
    void build(const DataMatrix& data, faiss::IndexIVFFlat* index, int k);
    void build(const DataMatrix& data, faiss::gpu::GpuIndexFlatL2* index, int k);
    void build(const DataMatrix& data, AnnoyIndex* index, int k);

    std::vector<std::vector<int>> neighbours;
    int numNeighbours;
};
