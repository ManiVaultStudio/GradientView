#pragma once

#include "Set.h"
#include "PointData/PointData.h"

#include <Eigen/Eigen>
#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>

using DataMatrix = Eigen::MatrixXf;

void convertToEigenMatrix(hdps::Dataset<Points> dataset, hdps::Dataset<Points> sourceDataset, DataMatrix& dataMatrix);

void convertToEigenMatrixProjection(hdps::Dataset<Points> dataset, DataMatrix& dataMatrix);
