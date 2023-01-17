#pragma once

#include <Eigen/Eigen>

#include <iostream>

using DataMatrix = Eigen::MatrixXf;

void standardizeData(DataMatrix& dataMatrix, std::vector<float>& variances);
void normalizeData(const DataMatrix& dataMatrix, std::vector<std::vector<float>>& normalizedData);
