#pragma once

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
    void lineClicked(int dim);

private slots:
    void onHover(const QPointF& point, bool state);
    void onClicked(const QPointF& point);

private:
    int _penWidth;
    int _dim;
};

class GradientGraph : public QWidget
{
    Q_OBJECT
public:
    GradientGraph();

    void setNumDimensions(int numDimensions);
    void setValues(const std::vector<std::vector<float>>& values);
    void setBins(const std::vector<std::vector<int>>& bins);
    void setTopDimensions(int dimension1, int dimension2);

    void updateChartColors();

signals:
    void lineClicked(int dim);

private slots:
    void onLineClicked(int dim);

private:
    int _numDimensions;
    std::vector<LineSeries*> _seriesArray;
    QChart* _chart;
    QChartView* _chartView;
    QValueAxis* _xAxis;
    QValueAxis* _yAxis;
    int _topDimension1;
    int _topDimension2;
    int _selectedDimension;
};
