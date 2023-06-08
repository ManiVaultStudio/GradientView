#include "KnnIndex.h"

#include "CudaCheck.h"

#include <QDebug>
#include <iostream>
#include <iomanip>

//void createBruteGraph(const DataMatrix& highDim, BruteIndexM* brute)
//{
//    Eigen::MatrixXf tarray = highDim.transpose();
//    brute = new knncpp::BruteForce<float, knncpp::ManhattenDistance<float>>(tarray, true);
//    brute->setThreads(-1);
//    brute->build();
//
//    qDebug() << "Rows:" << tarray.rows() << tarray.cols();
//}
//
//void createFaissGraph(const DataMatrix& highDim, faiss::IndexFlat* index)
//{
//    Eigen::MatrixXf tarray = highDim.transpose();
//
//    qDebug() << "Rows:" << tarray.rows() << tarray.cols();
//
//    std::vector<float> highDimArray(highDim.rows() * highDim.cols());
//    //Eigen::Map<Eigen::MatrixXf>(highDimArray.data(), highDim.rows(), highDim.cols()) = highDim;
//
//    int idx = 0;
//    for (int i = 0; i < highDim.rows(); i++)
//    {
//        for (int d = 0; d < highDim.cols(); d++)
//        {
//            highDimArray[idx++] = highDim(i, d);
//        }
//    }
//
//    index = new faiss::IndexFlat(highDim.cols(), faiss::METRIC_L2);
//    printf("is_trained = %s\n", index->is_trained ? "true" : "false");
//    index->add(highDim.rows(), highDimArray.data());                     // add vectors to the index
//    printf("ntotal = %ld\n", index->ntotal);
//}
//
//void createFaissIVFGraph(const DataMatrix& highDim, faiss::IndexFlat* quantizer, faiss::IndexIVFFlat* index)
//{
//    Eigen::MatrixXf tarray = highDim.transpose();
//
//    qDebug() << "Rows:" << tarray.rows() << tarray.cols();
//
//    std::vector<float> highDimArray(highDim.rows() * highDim.cols());
//    //Eigen::Map<Eigen::MatrixXf>(highDimArray.data(), highDim.rows(), highDim.cols()) = highDim;
//
//    int idx = 0;
//    for (int i = 0; i < highDim.rows(); i++)
//    {
//        for (int d = 0; d < highDim.cols(); d++)
//        {
//            highDimArray[idx++] = highDim(i, d);
//        }
//    }
//
//    int numPoints = highDim.rows();
//    int numDimensions = highDim.cols();
//
//    quantizer = new faiss::IndexFlatL2(numDimensions);
//    index = new faiss::IndexIVFFlat(quantizer, numDimensions, 100);
//    assert(!_faissIvfIndex->is_trained);
//    index->train(numPoints, highDimArray.data());
//    assert(index->is_trained);
//    printf("is_trained = %s\n", index->is_trained ? "true" : "false");
//    index->add(numPoints, highDimArray.data());                     // add vectors to the index
//    printf("ntotal = %ld\n", index->ntotal);
//}

void createFaissGpuGraph(faiss::gpu::StandardGpuResources*& res, faiss::gpu::GpuIndexFlat*& index, int numDimensions, knn::Metric metric)
{
    res = new faiss::gpu::StandardGpuResources();

    switch (metric)
    {
    case knn::Metric::MANHATTAN:
        index = new faiss::gpu::GpuIndexFlat(res, numDimensions, faiss::METRIC_L1); break;
    case knn::Metric::EUCLIDEAN:
        index = new faiss::gpu::GpuIndexFlat(res, numDimensions, faiss::METRIC_L2); break;
    }
}

void createAnnoyIndex(AnnoyIndex*& index, int numDimensions)
{
    index = new AnnoyIndex(numDimensions);
}

namespace knn
{
    Index::Index()
    {
        // Establish whether CUDA is available
        _hasCudaCapableGpu = checkCuda();
        if (_hasCudaCapableGpu)
            qDebug() << "CUDA capable device detected, using fast KNN.";
        else
            qDebug() << "No CUDA device detected, using CPU KNN.";
    }

    void Index::create(int numDimensions, Metric metric)
    {
        if (_hasCudaCapableGpu)
        {
            createFaissGpuGraph(_res, _faissGpuIndex, numDimensions, metric);
        }
        else
        {
            createAnnoyIndex(_annoyIndex, numDimensions);
        }
    }

    void Index::addData(const DataMatrix& data)
    {
        size_t numPoints = data.rows();
        size_t numDimensions = data.cols();

        size_t vectorSize = numPoints * numDimensions;
        std::vector<float> highDimArray(vectorSize);
        bool angular = numDimensions > 200;
        std::cout << "ANGULAR: " << angular;

        int idx = 0;
        for (int i = 0; i < numPoints; i++)
        {
            float len = 0;
            for (int d = 0; d < numDimensions; d++)
            {
                len += data(i, d) * data(i, d);
            }
            len = sqrt(len);
            for (int d = 0; d < numDimensions; d++)
            {
                if (angular)
                    highDimArray[idx++] = data(i, d) / len;
                else
                    highDimArray[idx++] = data(i, d);
            }
        }

        if (_hasCudaCapableGpu)
        {
            _faissGpuIndex->add(numPoints, highDimArray.data());
        }
        else
        {
            for (size_t i = 0; i < numPoints; ++i) {
                _annoyIndex->add_item(i, highDimArray.data() + (i * numDimensions));
                if (i % 1000 == 0)
                    std::cout << "Loading objects ...\t object: " << i + 1 << "\tProgress:" << std::fixed << std::setprecision(2) << (double)i / (double)(numPoints + 1) * 100 << "%\r";
            }

            std::cout << "Building index.." << std::endl;
            _annoyIndex->build(20 * numDimensions);
        }
    }

    void Index::search(const DataMatrix& data, int k, std::vector<int>& indices, std::vector<float>& distances) const
    {
        int numPoints = data.rows();
        int numDimensions = data.cols();

        // Put eigen matrix into flat float vector
        std::vector<float> highDimArray(numPoints * numDimensions);

        int idx = 0;
        for (int i = 0; i < numPoints; i++)
            for (int d = 0; d < numDimensions; d++)
                highDimArray[idx++] = data(i, d);

        // Initialize result vectors
        size_t resultSize = numPoints * k;
        indices.resize(resultSize);
        distances.resize(resultSize);

        if (_hasCudaCapableGpu)
        {
            idx_t* I = new idx_t[k * data.rows()];

            _faissGpuIndex->search(data.rows(), highDimArray.data(), k, distances.data(), I);

            indices.assign(I, I + resultSize);

            delete[] I;
        }
        else
        {
#pragma omp parallel for
            for (int i = 0; i < data.rows(); i++)
            {
                //if (i % 1000 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
                std::vector<int> closest;
                std::vector<float> closestD;

                _annoyIndex->get_nns_by_item(i, k, -1, &closest, &closestD);

                int startIndex = i * k;
                for (int j = 0; j < closest.size(); j++)
                {
                    indices[startIndex + j] = closest[j];
                    distances[startIndex + j] = closestD[j];
                }
            }
        }
    }
}
