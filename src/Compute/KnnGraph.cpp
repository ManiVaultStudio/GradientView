#include "KnnGraph.h"

#include <iostream>

void KnnGraph::build(const DataMatrix& data, KdTree* kdTree, int k)
{
    numNeighbours = k - 1;

    Matrixi indices;
    Eigen::MatrixXf distances;

    Eigen::MatrixXf arrayt = data.transpose();
    kdTree->query(arrayt, k, indices, distances);

    neighbours.resize(data.rows(), std::vector<int>(numNeighbours));

    for (int i = 0; i < data.rows(); i++)
    {
        if (i % 100 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
        for (int j = 0; j < numNeighbours; j++)
        {
            neighbours[i][j] = indices(j + 1, i);
        }
    }
}
void KnnGraph::build(const DataMatrix& data, AnnoyIndex* index, int k)
{
    std::vector<float> highDimArray(data.rows() * data.cols());
    //Eigen::Map<Eigen::MatrixXf>(highDimArray.data(), data.rows(), data.cols()) = data;

    int idx = 0;
    for (int i = 0; i < data.rows(); i++)
    {
        for (int d = 0; d < data.cols(); d++)
        {
            highDimArray[idx++] = data(i, d);
        }
    }

    std::vector<int> indices(data.rows() * (k + 1));
    std::vector<float> distances(data.rows() * (k + 1));

#pragma omp parallel for
    for (int i = 0; i < data.rows(); i++)
    {
        //if (i % 1000 == 0) std::cout << "Querying neighbours: " << i << "/" << data.rows() << std::endl;
        std::vector<int> closest;
        std::vector<float> closestD;
        
        index->get_nns_by_item(i, k+1, -1, &closest, &closestD);

        int startIndex = i * (k + 1);
        for (int j = 0; j < closest.size(); j++)
        {
            indices[startIndex + j] = closest[j];
            distances[startIndex + j] = closestD[j];
        }
    }

    printIndices("ANNOY", indices, k);
    printDistances("ANNOY", distances.data(), k);

    //printf("I (5 last results)=\n");
    //for (int i = data.rows() - 5; i < data.rows(); i++) {
    //    for (int j = 0; j < k; j++)
    //        printf("%5zd ", indices[i * (k + 1) + j]);
    //    printf("\n");
    //}

    numNeighbours = k;
    neighbours.clear();
    neighbours.resize(data.rows(), std::vector<int>(numNeighbours));

#pragma omp parallel for
    for (int i = 0; i < data.rows(); i++)
    {
        //if (i % 100 == 0) std::cout << "Building neighbours: " << i << "/" << data.rows() << std::endl;
        for (int j = 0; j < numNeighbours; j++)
        {
            neighbours[i][j] = indices[i * (k+1) + j + 1];
        }
    }

}
