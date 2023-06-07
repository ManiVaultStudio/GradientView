#pragma once

#include "KnnGraph.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

void exportKnnGraph(const KnnGraph& graph)
{
    const std::vector<std::vector<int>>& neighbours = graph.getNeighbours();
    int numPoints = neighbours.size();
    int numNeighbours = graph.getNumNeighbours();

    // Linearize data
    std::vector<int> linearNeighbours(numPoints * numNeighbours);
    int c = 0;
    for (int i = 0; i < numPoints; i++)
    {
        const std::vector<int>& n = neighbours[i];
        for (int j = 0; j < numNeighbours; j++)
            linearNeighbours[c++] = neighbours[i][j];
    }
    std::cout << "Linearized data for export" << std::endl;
    // Write to file
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream fileName;
    fileName << "knngraph";
    fileName << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    fileName << ".knn";
    std::cout << "Writing to file: " << fileName.str() << std::endl;
    std::ofstream myfile(fileName.str(), std::ios::out | std::ios::binary);
    if (!myfile) {
        std::cout << "Cannot open file for writing KNN graph!" << std::endl;
        return;
    }
    myfile.write((char*) &numPoints, sizeof(int));
    myfile.write((char*) &numNeighbours, sizeof(int));

    for (int i = 0; i < linearNeighbours.size(); i++)
    {
        if (i % 10000 == 0) std::cout << "Progress: " << i << std::endl;
        myfile.write((char*) &linearNeighbours[i], sizeof(int));
    }

    myfile.close();
    std::cout << "Knn graph written to file" << std::endl;
}
