#pragma once

#include <knncpp.h>

#pragma warning(push, 0) // annoylib.h has some warnings that are 'annoy'ing
#define ANNOYLIB_MULTITHREADED_BUILD
#include <annoylib.h>
#include <kissrandom.h>
#pragma warning(pop)

#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/gpu/GpuIndexFlat.h>
#include <faiss/gpu/StandardGpuResources.h>

using BruteIndexE = knncpp::BruteForce<float, knncpp::EuclideanDistance<float>>;
using BruteIndexM = knncpp::BruteForce<float, knncpp::ManhattenDistance<float>>;
using KdTree = knncpp::KDTreeMinkowskiX<float, knncpp::ManhattenDistance<float>>;
using AnnoyIndex = Annoy::AnnoyIndex<int, float, Annoy::Angular, Annoy::Kiss32Random, Annoy::AnnoyIndexMultiThreadedBuildPolicy>;

typedef knncpp::Matrixi Matrixi;
using DataMatrix = Eigen::MatrixXf;

using idx_t = int64_t;

namespace knn
{
    enum class Metric
    {
        EUCLIDEAN, MANHATTAN, COSINE, ANGULAR
    };

    class Index
    {
    public:
        Index();

        void create(int numDimensions, Metric metric);
        void addData(const DataMatrix& data);
        void search(const DataMatrix& data, int numNeighbours, std::vector<int>& indices, std::vector<float>& distances) const;

    private:
        //BruteIndexE*                        _bruteEuclidean = nullptr;
        //BruteIndexM*                        _bruteManhattan = nullptr;
        //faiss::IndexFlat*                   _faissIndex     = nullptr;
        //faiss::IndexFlat*                   _quantizer      = nullptr;
        //faiss::IndexIVFFlat*                _faissIvfIndex  = nullptr;
        AnnoyIndex*                         _annoyIndex     = nullptr;
        faiss::gpu::StandardGpuResources*   _res            = nullptr;
        faiss::gpu::GpuIndexFlatL2*         _faissGpuIndex  = nullptr;

        bool _hasCudaCapableGpu = false;
    };
}
