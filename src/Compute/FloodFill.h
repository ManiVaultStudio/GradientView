#pragma once

#include "KnnGraph.h"

#include <vector>

class FloodFill
{
public:
    FloodFill(int numWaves);

    void compute(const KnnGraph& knnGraph, int selectedPoint);
    void recompute();

    void setNumWaves(int numWaves);

    size_t getNumWaves() const { return _waves.size(); }
    size_t getTotalNumNodes() const { return _allNodes.size(); }

    std::vector<std::vector<int>>& getWaves() { return _waves; }
    const std::vector<std::vector<int>>& getWaves() const { return _waves; }
    const std::vector<int>& getAllNodes() const { return _allNodes; }

private:
    int _numWaves;

    std::vector<std::vector<int>> _waves;

    std::vector<int> _allNodes;

    // Store knn graph for recomputation
    const KnnGraph* _lastKnnGraph;
    int _lastSelectedPoint;
    int _lastNumWaves;
};
