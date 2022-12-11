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
