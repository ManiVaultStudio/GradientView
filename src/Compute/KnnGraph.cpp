#include "KnnGraph.h"

#include <iostream>

void printIndices(std::string name, idx_t* indices, int k)
{
    std::cout << "Indices (5 first results): " << name << std::endl;

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < k+1; j++)
            printf("%5zd ", indices[i * (k + 1) + j]);
        printf("\n");
    }

    //printf("FAISS I (5 last results)=\n");
    //for (int i = data.rows() - 5; i < data.rows(); i++) {
    //    for (int j = 0; j < k; j++)
    //        printf("%5zd ", I[i * (k + 1) + j]);
    //    printf("\n");
    //}
}

void printIndices(std::string name, std::vector<int>& indices, int k)
{
    std::cout << "Indices (5 first results): " << name << std::endl;

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < k + 1; j++)
            printf("%5zd ", indices[i * (k + 1) + j]);
        printf("\n");
    }
}

void printDistances(std::string name, float* distances, int k)
{
    std::cout << "Distances (5 first results): " << name << std::endl;

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < k + 1; j++)
            printf("%7g ", distances[i * (k + 1) + j]);
        printf("\n");
    }
}

//void KnnGraph::build(const DataMatrix& data, BruteIndexM* brute, int k)
//{
//    numNeighbours = k;
//
//    Matrixi indices;
//    Eigen::MatrixXf distances;
//
//    Eigen::MatrixXf arrayt = data.transpose();
//    brute->query(arrayt, k + 1, indices, distances);
//
//    printf("Indices (5 first results)=: BRUTE\n");
//    for (int i = 0; i < 5; i++) {
//        std::cout << "i: " << i << std::endl;
//        for (int j = 0; j < k; j++)
//            printf("%5zd ", indices(j, i));
//        printf("\n");
//    }
//
//    std::cout << "Distances (5 first results): " << "BRUTE" << std::endl;
//    for (int i = 0; i < 5; i++) {
//        for (int j = 0; j < k + 1; j++)
//            printf("%7g ", distances(j, i));
//        printf("\n");
//    }
//
//    neighbours.resize(data.rows(), std::vector<int>(numNeighbours));
//
//    for (int i = 0; i < data.rows(); i++)
//    {
//        //if (i % 100 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
//        for (int j = 0; j < numNeighbours; j++)
//        {
//            neighbours[i][j] = indices(j + 1, i);
//        }
//    }
//}
//
//void KnnGraph::build(const DataMatrix& data, KdTree* kdTree, int k)
//{
//    numNeighbours = k;
//
//    Matrixi indices;
//    Eigen::MatrixXf distances;
//
//    Eigen::MatrixXf arrayt = data.transpose();
//    kdTree->query(arrayt, k+1, indices, distances);
//
//    printf("I (5 first results)=\n");
//    for (int i = 0; i < 5; i++) {
//        std::cout << "i: " << i << std::endl;
//        for (int j = 0; j < k; j++)
//            printf("%5zd ", indices(j, i));
//        printf("\n");
//    }
//
//    printf("I (5 last results)=\n");
//    for (int i = data.rows() - 5; i < data.rows(); i++) {
//        for (int j = 0; j < k; j++)
//            printf("%5zd ", indices(j, i));
//        printf("\n");
//    }
//
//    neighbours.resize(data.rows(), std::vector<int>(numNeighbours));
//
//    for (int i = 0; i < data.rows(); i++)
//    {
//        if (i % 100 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
//        for (int j = 0; j < numNeighbours; j++)
//        {
//            neighbours[i][j] = indices(j + 1, i);
//        }
//    }
//}
//
////void KnnGraph::build(const DataMatrix& data, flann::Index<flann::L2<float>>* index, int k)
////{
////    numNeighbours = k - 1;
////
////    Eigen::MatrixXf arrayt = data.transpose();
////
////    std::vector<float> highDimArray(data.rows() * data.cols());
////    //Eigen::Map<Eigen::MatrixXf>(highDimArray.data(), data.rows(), data.cols()) = data;
////
////    int idx = 0;
////    for (int i = 0; i < data.rows(); i++)
////    {
////        for (int d = 0; d < data.cols(); d++)
////        {
////            highDimArray[idx++] = data(i, d);
////        }
////    }
////    
////    std::vector<int> indices(data.rows() * k);
////    std::vector<float> distances_squared(data.rows() * k);
////    flann::Matrix<int> indices_mat(indices.data(), data.rows(), k);
////    flann::Matrix<float> dists_mat(distances_squared.data(), data.rows(), k);
////
////    flann::Matrix<float> query(highDimArray.data(), data.rows(), data.cols());
////
////    flann::SearchParams flann_params(128);
////    flann_params.cores = 0; //all cores
////    index->knnSearch(query, indices_mat, dists_mat, k, flann_params);
////
////    neighbours.resize(data.rows(), std::vector<int>(numNeighbours));
////
////    for (int i = 0; i < data.rows(); i++)
////    {
////        if (i % 100 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
////        for (int j = 0; j < numNeighbours; j++)
////        {
////            neighbours[i][j] = indices_mat[i][j + 1];
////        }
////    }
////}
//
//void KnnGraph::build(const DataMatrix& data, faiss::IndexFlat* index, int k)
//{
//    std::vector<float> highDimArray(data.rows() * data.cols());
//    //Eigen::Map<Eigen::MatrixXf>(highDimArray.data(), data.rows(), data.cols()) = data;
//
//    int idx = 0;
//    for (int i = 0; i < data.rows(); i++)
//    {
//        for (int d = 0; d < data.cols(); d++)
//        {
//            highDimArray[idx++] = data(i, d);
//        }
//    }
//
//    idx_t* I = new idx_t[(k+1) * data.rows()];
//    float* D = new float[(k+1) * data.rows()];
//
//    index->search(data.rows(), highDimArray.data(), k+1, D, I);
//
//    // print results
//    printIndices("FAISS CPU", I, k);
//    printDistances("FAISS CPU", D, k);
//
//    numNeighbours = k;
//    neighbours.clear();
//    neighbours.resize(data.rows(), std::vector<int>(numNeighbours));
//
//    for (int i = 0; i < data.rows(); i++)
//    {
//        //if (i % 100 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
//        for (int j = 0; j < numNeighbours; j++)
//        {
//            neighbours[i][j] = I[i * (k+1) + j + 1];
//        }
//    }
//
//    delete[] I;
//    delete[] D;
//}
//
//void KnnGraph::build(const DataMatrix& data, faiss::IndexIVFFlat* index, int k)
//{
//    std::vector<float> highDimArray(data.rows() * data.cols());
//    //Eigen::Map<Eigen::MatrixXf>(highDimArray.data(), data.rows(), data.cols()) = data;
//
//    int idx = 0;
//    for (int i = 0; i < data.rows(); i++)
//    {
//        for (int d = 0; d < data.cols(); d++)
//        {
//            highDimArray[idx++] = data(i, d);
//        }
//    }
//
//    idx_t* I = new idx_t[k * data.rows()];
//    float* D = new float[k * data.rows()];
//
//    index->nprobe = 10;
//    index->search(data.rows(), highDimArray.data(), k, D, I);
//
//    // print results
//    printf("I (5 first results)=\n");
//    for (int i = 0; i < 5; i++) {
//        for (int j = 0; j < k; j++)
//            printf("%5zd ", I[i * k + j]);
//        printf("\n");
//    }
//
//    printf("I (5 last results)=\n");
//    for (int i = data.rows() - 5; i < data.rows(); i++) {
//        for (int j = 0; j < k; j++)
//            printf("%5zd ", I[i * k + j]);
//        printf("\n");
//    }
//
//    numNeighbours = k - 1;
//    neighbours.clear();
//    neighbours.resize(data.rows(), std::vector<int>(numNeighbours));
//
//    for (int i = 0; i < data.rows(); i++)
//    {
//        if (i % 100 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
//        for (int j = 0; j < numNeighbours; j++)
//        {
//            neighbours[i][j] = I[i * k + j + 1];
//        }
//    }
//
//    delete[] I;
//    delete[] D;
//}

// Build KNN sub-graph from bigger graph
void KnnGraph::build(const KnnGraph& graph, int numNeighbours)
{
    assert(graph.getNumNeighbours() > numNeighbours);

    _numNeighbours = numNeighbours;

    const auto& graphNeighbours = graph.getNeighbours();
    _neighbours.resize(graphNeighbours.size(), std::vector<int>(numNeighbours));

    for (int i = 0; i < _neighbours.size(); i++)
    {
        for (int k = 0; k < numNeighbours; k++)
        {
            _neighbours[i][k] = graphNeighbours[i][k];
        }
    }
}

void KnnGraph::build(const DataMatrix& data, const knn::Index& index, int numNeighbours)
{
    std::vector<int> indices;
    std::vector<float> distances;

    int k = numNeighbours + 1; // Plus one to account for the query point itself being in the results

    index.search(data, k, indices, distances);

    // print results
    printIndices("INDEX", indices, k);
    printDistances("INDEX", distances.data(), k);

    _numNeighbours = numNeighbours;
    _neighbours.clear();
    _neighbours.resize(data.rows(), std::vector<int>(_numNeighbours));

#pragma omp parallel for
    for (int i = 0; i < data.rows(); i++)
    {
        //if (i % 100 == 0) std::cout << "Building graph: " << i << "/" << data.rows() << std::endl;
        for (int j = 0; j < _numNeighbours; j++)
        {
            _neighbours[i][j] = indices[i * k + j + 1];
        }
    }
}
