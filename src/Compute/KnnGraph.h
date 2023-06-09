#pragma once

#include "KnnIndex.h"

#include <QString>

class KnnGraphImporter;
class KnnGraphExporter;

class KnnGraph
{
public:
    KnnGraph();

    const std::vector<std::vector<int>>& getNeighbours() const { return _neighbours; }
    size_t getNumNeighbours() const { return _numNeighbours; }

    void build(const KnnGraph& graph, size_t numNeighbours);
    void build(const DataMatrix& data, const knn::Index& index, size_t numNeighbours);
    void build(const KnnGraph& graph, size_t numNeighbours, bool shared);

    void readFromFile(QString filePath);
    void writeToFile();

private:
    std::vector<std::vector<int>> _neighbours;
    size_t _numNeighbours;

    friend class ScatterplotPlugin;
    friend class KnnGraphImporter;
    friend class KnnGraphExporter;
};
