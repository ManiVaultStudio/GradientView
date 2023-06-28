#pragma once

#include "DataMatrix.h"

// Eigen::IndexedView<Eigen::MatrixXf, std::vector<int>, Eigen::internal::AllRange<-1>>

class DataStorage
{
public:
    DataMatrix& getData() { return _currentDataView; }
    DataMatrix& getFullData() { return _dataMatrix; }
    DataMatrix& getProjection() { return _currentProjView; }
    DataMatrix& getFullProjection() { return _currentFullProjView; }
    std::vector<float>& getVariances() { return _variancesView; }

    const std::vector<int>& getViewIndices() const { return _viewIndices; }

    int getNumPoints() { return _currentDataView.rows(); }
    int getNumDimensions() { return _currentDataView.cols(); }

    void setProjectionSize(float projectionSize) { _projectionSize = projectionSize; }
    float getProjectionSize() { return _projectionSize; }

    void storeState()
    {
        _dataMatrix = _currentDataView;
        _fullProjMatrix = _currentFullProjView;
        _projMatrix = _currentProjView;

        //_variances = _variancesView;
    }

    void createDataView(const std::vector<int>& indices)
    {
        _currentDataView = _dataMatrix(indices, Eigen::all);
        _currentFullProjView = _fullProjMatrix(indices, Eigen::all);
        _currentProjView = _projMatrix(indices, Eigen::all);

        //_variancesView.resize(indices.size());
        //for (int i = 0; i < indices.size(); i++)
        //{
        //    _variancesView[i] = _variances[indices[i]];
        //}

        _viewIndices = indices;
    }

private:
    // Stored state
    // Data
    DataMatrix                      _dataMatrix;
    std::vector<float>              _variances;

    // Projection
    DataMatrix                      _fullProjMatrix;
    DataMatrix                      _projMatrix;
    float                           _projectionSize = 0;

    // View
    DataMatrix                      _currentDataView;
    std::vector<float>              _variancesView;

    DataMatrix                      _currentFullProjView;
    DataMatrix                      _currentProjView;

    std::vector<int>                _viewIndices;
};
