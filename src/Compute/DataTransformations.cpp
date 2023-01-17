#include "DataTransformations.h"

void standardizeData(DataMatrix& dataMatrix, std::vector<float>& variances)
{
    // Compute means
    std::vector<float> means(dataMatrix.cols());
    for (int d = 0; d < dataMatrix.cols(); d++)
    {
        means[d] = 0;
        for (int i = 0; i < dataMatrix.rows(); i++)
        {
            means[d] += dataMatrix(i, d);
        }
        means[d] /= dataMatrix.rows();
        std::cout << "Mean " << d << " " << means[d] << std::endl;
    }

    // Compute variances
    variances.resize(dataMatrix.cols());
    for (int d = 0; d < dataMatrix.cols(); d++)
    {
        variances[d] = 0;
        for (int i = 0; i < dataMatrix.rows(); i++)
        {
            variances[d] += (dataMatrix(i, d) - means[d]) * (dataMatrix(i, d) - means[d]);
        }
        variances[d] /= dataMatrix.rows();
        std::cout << "Variance " << d << " " << variances[d] << std::endl;
    }
    for (int d = 0; d < dataMatrix.cols(); d++)
    {
        if (variances[d] <= 0) continue;
        for (int i = 0; i < dataMatrix.rows(); i++)
        {
            dataMatrix(i, d) -= means[d];
            dataMatrix(i, d) /= sqrt(variances[d]);
        }
    }
}

void normalizeData(const DataMatrix& dataMatrix, std::vector<std::vector<float>>& normalizedData)
{
    normalizedData.resize(dataMatrix.cols(), std::vector<float>(dataMatrix.rows()));
    for (int d = 0; d < dataMatrix.cols(); d++)
    {
        auto col = dataMatrix.col(d);

        float minVal = *std::min_element(col.begin(), col.end());
        float maxVal = *std::max_element(col.begin(), col.end());
        float range = maxVal - minVal;
        if (range == 0) range = 1;

        for (int i = 0; i < dataMatrix.rows(); i++)
            normalizedData[d][i] = std::min(0.99999f, (col(i) - minVal) / range);
    }
}
