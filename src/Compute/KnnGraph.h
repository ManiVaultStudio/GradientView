#pragma once

#include "KnnIndex.h"

//#include <cstdlib>
//#include <stdint.h>

#include <Eigen/Eigen>

class KnnGraph
{
public:
    const std::vector<std::vector<int>>& getNeighbours() const { return _neighbours; }
    int getNumNeighbours() const { return _numNeighbours; }

    //void build(const DataMatrix& data, BruteIndexM* brute, int k);
    //void build(const DataMatrix& data, KdTree* kdTree, int k);
    ////void build(const DataMatrix& data, flann::Index<flann::L2<float>>* index, int k);
    //void build(const DataMatrix& data, faiss::IndexFlat* index, int k);
    //void build(const DataMatrix& data, faiss::IndexIVFFlat* index, int k);
    //void build(const DataMatrix& data, faiss::gpu::GpuIndexFlatL2* index, int k);
    //void build(const DataMatrix& data, AnnoyIndex* index, int k);

    void build(const DataMatrix& data, const knn::Index& index, int numNeighbours);

    std::vector<std::vector<int>> _neighbours;
    int _numNeighbours;
};
