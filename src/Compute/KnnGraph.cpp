#include "KnnGraph.h"

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
        for (int j = 0; j < k - 1; j++)
        {
            neighbours[i][j] = indices(j + 1, i);
        }
    }
}
