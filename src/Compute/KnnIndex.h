#pragma once

#include "DataMatrix.h"

#pragma warning(push, 0) // annoylib.h has some warnings that are 'annoy'ing
#define ANNOYLIB_MULTITHREADED_BUILD
#include <annoylib.h>
#include <kissrandom.h>
#pragma warning(pop)

#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/gpu/GpuIndexFlat.h>
#include <faiss/gpu/StandardGpuResources.h>

using AnnoyIndex = Annoy::AnnoyIndex<int, float, Annoy::Angular, Annoy::Kiss32Random, Annoy::AnnoyIndexMultiThreadedBuildPolicy>;

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
        void linearizeData(const DataMatrix& data, std::vector<float>& highDimArray) const;

    private:
        AnnoyIndex*                         _annoyIndex     = nullptr;
        faiss::gpu::StandardGpuResources*   _res            = nullptr;
        faiss::gpu::GpuIndexFlat*           _faissGpuIndex  = nullptr;

        Metric                              _metric;

        bool _usingGpuCompute = false;
    };
}
