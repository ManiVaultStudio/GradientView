#pragma once

#include <QWidget>

#include <Eigen/Eigen>

#include <QtCharts/QLineSeries>
#include <QtCharts/QChartView>
#include <QChart>

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

    void setDimension(const DataMatrix& data, int dim);
    void setValues(const std::vector<std::vector<float>>& values);

signals:
    void lineClicked(int dim);

private slots:
    void onLineClicked(int dim);

private:
    QLineSeries* _series;
    QChart* _chart;
    QChartView* _chartView;
};
