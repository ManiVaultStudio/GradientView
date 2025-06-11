#include "KnnIndex.h"

#include <QDebug>
#include <iostream>
#include <iomanip>

#include <fstream>
#include <sstream>

namespace
{
    void linearizeData(const DataMatrix& data, std::vector<float>& highDimArray, knn::Metric metric)
    {
        size_t numPoints = data.rows();
        size_t numDimensions = data.cols();

        // Put eigen matrix into flat float vector
        highDimArray.resize(numPoints * numDimensions);

        int idx = 0;
        for (int i = 0; i < numPoints; i++)
        {
            if (metric == knn::Metric::COSINE)
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
            if (i % 1000000 == 0) qDebug() << "Convert Data Progress: " << i;
        }
    }
}

void writeDataMatrixToDisk(const DataMatrix& dataMatrix)
{
    uint32_t numPoints = (uint32_t)dataMatrix.rows();
    uint32_t numDimensions = (uint32_t)dataMatrix.cols();

    // Linearize data
    std::vector<float> linearData(numPoints * numDimensions);
    int c = 0;
    for (int i = 0; i < numPoints; i++)
    {
        for (int j = 0; j < numDimensions; j++)
            linearData[c++] = dataMatrix(i, j);
    }
    std::cout << "Linearized data for export" << std::endl;
    // Write to file
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream fileName;
    fileName << "data_matrix";
    fileName << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    fileName << ".data";
    std::cout << "Writing to file: " << fileName.str() << std::endl;
    std::ofstream myfile(fileName.str(), std::ios::out | std::ios::binary);
    if (!myfile) {
        std::cout << "Cannot open file for writing data matrix!" << std::endl;
        return;
    }
    myfile.write((char*)&numPoints, sizeof(uint32_t));
    myfile.write((char*)&numDimensions, sizeof(uint32_t));

    for (size_t i = 0; i < linearData.size(); i++)
    {
        if (i % 10000 == 0) std::cout << "Progress: " << i << std::endl;
        myfile.write((char*)&linearData[i], sizeof(uint32_t));
    }

    myfile.close();
    std::cout << "Data matrix written to file" << std::endl;
}

void createFaissIndex(faiss::IndexFlat*& index, int numDimensions, knn::Metric metric)
{
    qDebug() << "Creating faiss index";
    switch (metric)
    {
    case knn::Metric::MANHATTAN:
        index = new faiss::IndexFlat(numDimensions, faiss::METRIC_L1); break;
    case knn::Metric::EUCLIDEAN:
        index = new faiss::IndexFlat(numDimensions, faiss::METRIC_L2); break;
    case knn::Metric::COSINE:
        index = new faiss::IndexFlat(numDimensions, faiss::METRIC_INNER_PRODUCT); break;
    }
}

void createFaissIVFIndex(const DataMatrix& data, faiss::IndexFlatL2*& quantizer, faiss::IndexIVFFlat*& index, int numDimensions, knn::Metric metric)
{
    qDebug() << "Creating IVF index";
    int nlist = sqrt(data.rows());
    quantizer = new faiss::IndexFlatL2(numDimensions);
    index = new faiss::IndexIVFFlat(quantizer, numDimensions, nlist);
    index->nprobe = 10;

    std::vector<float> indexData;
    linearizeData(data, indexData, metric);

    qDebug() << "Training IVF index..";
    index->train(data.rows(), indexData.data());
    qDebug() << "Trained IVF index.";
}

namespace knn
{
    Index::Index() :
        _metric(Metric::EUCLIDEAN)
    {

    }

    Index::~Index()
    {
        if (_annoyIndex != nullptr)
            delete _annoyIndex;
        if (_faissIndex != nullptr)
            delete _faissIndex;
    }

    void Index::create(const DataMatrix& data, int numDimensions, Metric metric)
    {
        _metric = metric;
        if (_preciseKnn)
            createFaissIndex(_faissIndex, numDimensions, metric);
        else
            createFaissIVFIndex(data, _ivfQuantizer, _faissIVFIndex, numDimensions, metric);
    }

    void Index::addData(const DataMatrix& data)
    {
        size_t numPoints = data.rows();
        size_t numDimensions = data.cols();

        std::vector<float> indexData;
        linearizeData(data, indexData, _metric);

        if (_preciseKnn)
        {
            //writeDataMatrixToDisk(data);
            _faissIndex->add(numPoints, indexData.data());
        }
        else
        {
            _faissIVFIndex->add(numPoints, indexData.data());
        }
    }

    void Index::search(const DataMatrix& data, int k, std::vector<int>& indices, std::vector<float>& distances) const
    {
        int numPoints = data.rows();
        int numDimensions = data.cols();

        // Put eigen matrix into flat float vector
        std::vector<float> query;
        linearizeData(data, query, _metric);

        // Initialize result vectors
        size_t resultSize = numPoints * k;
        indices.resize(resultSize);
        distances.resize(resultSize);

        if (_preciseKnn)
        {
            idx_t* I = new idx_t[resultSize];

            _faissIndex->search(numPoints, query.data(), k, distances.data(), I);

            indices.assign(I, I + resultSize);

            delete[] I;
        }
        else
        {
            idx_t* I = new idx_t[resultSize];
            qDebug() << "Performing search...";
            _faissIVFIndex->search(numPoints, query.data(), k, distances.data(), I);
            qDebug() << "Done with searching!";
            indices.assign(I, I + resultSize);

            delete[] I;
        }
    }
}
