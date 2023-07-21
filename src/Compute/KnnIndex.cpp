#include "KnnIndex.h"

#include "CudaCheck.h"

#include <QDebug>
#include <iostream>
#include <iomanip>

void createFaissGpuIndex(faiss::gpu::StandardGpuResources*& res, faiss::gpu::GpuIndexFlat*& index, int numDimensions, knn::Metric metric)
{
    res = new faiss::gpu::StandardGpuResources();

    switch (metric)
    {
    case knn::Metric::MANHATTAN:
        index = new faiss::gpu::GpuIndexFlat(res, numDimensions, faiss::METRIC_L1); break;
    case knn::Metric::EUCLIDEAN:
        index = new faiss::gpu::GpuIndexFlat(res, numDimensions, faiss::METRIC_L2); break;
    case knn::Metric::COSINE:
        index = new faiss::gpu::GpuIndexFlat(res, numDimensions, faiss::METRIC_INNER_PRODUCT); break;
    }
}

void createAnnoyIndex(AnnoyIndex*& index, int numDimensions)
{
    index = new AnnoyIndex(numDimensions);
}

namespace knn
{
    Index::Index() :
        _metric(Metric::EUCLIDEAN),
        _usingGpuCompute(false)
    {
        // Establish whether CUDA is available
        bool hasCudaCapableGpu = checkCuda();
        if (hasCudaCapableGpu)
        {
            qDebug() << "CUDA capable device detected, using fast KNN.";
            _usingGpuCompute = true;
        }
        else
            qDebug() << "No CUDA device detected, using CPU KNN.";
    }

    void Index::create(int numDimensions, Metric metric)
    {
        _metric = metric;
        if (_usingGpuCompute)
            createFaissGpuIndex(_res, _faissGpuIndex, numDimensions, metric);
        else
            createAnnoyIndex(_annoyIndex, numDimensions);
    }

    void Index::addData(const DataMatrix& data)
    {
        size_t numPoints = data.rows();
        size_t numDimensions = data.cols();

        std::vector<float> indexData;
        linearizeData(data, indexData);

        if (_usingGpuCompute)
        {
            _faissGpuIndex->add(numPoints, indexData.data());
        }
        else
        {
            _annoyIndex->load("test.ann");
            //writeDataMatrixToDisk(data);
            //_annoyIndex->on_disk_build("test.ann");
            //for (size_t i = 0; i < numPoints; ++i) {
            //    if (i % 10000 == 0) qDebug() << "Add Progress: " << i;
            //    _annoyIndex->add_item((int) i, highDimArray.data() + (i * numDimensions));
            //    if (i % 100 == 0)
            //        std::cout << "Loading objects ...\t object: " << i + 1 << "\tProgress:" << std::fixed << std::setprecision(2) << (double)i / (double)(numPoints + 1) * 100 << "%\r";
            //}

            //std::cout << "Building index.." << std::endl;
            //_annoyIndex->build((int) (10 * numDimensions));
            //_annoyIndex->save("test.ann");
        }
    }

    void Index::search(const DataMatrix& data, int k, std::vector<int>& indices, std::vector<float>& distances) const
    {
        int numPoints = data.rows();
        int numDimensions = data.cols();

        // Put eigen matrix into flat float vector
        std::vector<float> query;
        linearizeData(data, query);

        // Initialize result vectors
        size_t resultSize = numPoints * k;
        indices.resize(resultSize);
        distances.resize(resultSize);

        if (_usingGpuCompute)
        {
            idx_t* I = new idx_t[k * data.rows()];

            _faissGpuIndex->search(data.rows(), query.data(), k, distances.data(), I);

            indices.assign(I, I + resultSize);

            delete[] I;
        }
        else
        {
            std::vector<std::vector<int>> tempIndices(numPoints, std::vector<int>(k));
            std::vector<std::vector<float>> tempDistances(numPoints, std::vector<float>(k));

            int totalCount = 0;
#pragma omp parallel for
            for (int i = 0; i < data.rows(); i++)
            {
                _annoyIndex->get_nns_by_item(i, k, -1, &tempIndices[i], &tempDistances[i]);
                
                if (i % 1000 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
//#pragma omp critical
//                {
//                    totalCount += 1;
//                    if (totalCount % 1000 == 0) std::cout << "Querying neighbours: " << totalCount << "/" << data.rows() << std::endl;
//                }
            }
        }
    }

    void Index::linearizeData(const DataMatrix& data, std::vector<float>& highDimArray) const
    {
        size_t numPoints = data.rows();
        size_t numDimensions = data.cols();

        // Put eigen matrix into flat float vector
        highDimArray.resize(numPoints * numDimensions);

        int idx = 0;
        for (int i = 0; i < numPoints; i++)
        {
            if (_metric == Metric::COSINE)
            {
                double len = 0;
                for (int d = 0; d < numDimensions; d++)
                {
                    double dd = data(i, d);
                    len += dd * dd;
                }
                len = sqrt(len);

                for (int d = 0; d < numDimensions; d++)
                {
                    highDimArray[idx++] = data(i, d) / len;
                }
            }
            else
            {
                for (int d = 0; d < numDimensions; d++)
                    highDimArray[idx++] = data(i, d);
            }
            if (i % 10000 == 0) qDebug() << "Convert Data Progress: " << i;
        }
    }
}
