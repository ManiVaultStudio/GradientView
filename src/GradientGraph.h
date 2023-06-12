#pragma once

#include "Types.h"

#include <QWidget>

#include <Eigen/Eigen>

#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QChartView>
#include <QChart>

#include <vector>

using DataMatrix = Eigen::MatrixXf;

class LineSeries : public QLineSeries
{
    Q_OBJECT
public:
    LineSeries(QObject* parent, int dim);

signals:
    void lineClicked(dint dim);

private slots:
    void onHover(const QPointF& point, bool state);
    void onClicked(const QPointF& point);

private:
    int _penWidth;
    dint _dim;
};

class GradientGraph : public QWidget
{
    Q_OBJECT
public:
    GradientGraph();

    void setNumDimensions(dint numDimensions);
    void setValues(const std::vector<std::vector<float>>& values);
    void setBins(const std::vector<std::vector<int>>& bins);
    void setTopDimensions(dint dimension1, dint dimension2);

    void updateChartColors();

signals:
    void lineClicked(dint dim);

private slots:
    void onLineClicked(dint dim);

private:
    dint _numDimensions;
    std::vector<LineSeries*> _seriesArray;
    QChart* _chart;
    QChartView* _chartView;
    QValueAxis* _xAxis;
    QValueAxis* _yAxis;
    dint _topDimension1;
    dint _topDimension2;
    dint _selectedDimension;
};
