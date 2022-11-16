#include "GradientGraph.h"

#include <QVBoxLayout>

#include <iostream>
#include <algorithm>

LineSeries::LineSeries(QObject* parent, int dim) :
    QLineSeries(parent),
    _penWidth(2),
    _dim(dim)
{
    connect(this, &QXYSeries::hovered, this, &LineSeries::onHover);
    connect(this, &QXYSeries::clicked, this, &LineSeries::onClicked);
}

void LineSeries::onHover(const QPointF& point, bool state)
{
    QPen p = pen();
    p.setWidth(state ? _penWidth * 2 : _penWidth);
    setPen(p);
}

void LineSeries::onClicked(const QPointF& point)
{
    emit lineClicked(_dim);
}

GradientGraph::GradientGraph() :
    _numDimensions(0),
    _chart(new QChart()),
    _chartView(nullptr)
{
    _chart->legend()->hide();

    _xAxis = new QValueAxis();
    _xAxis->setRange(0, 1);
    _xAxis->setTickCount(10);
    _xAxis->setLabelFormat("%d");
    _chart->addAxis(_xAxis, Qt::AlignBottom);

    _yAxis = new QValueAxis();
    _yAxis->setRange(0, 1);
    _yAxis->setTickCount(10);
    _yAxis->setLabelFormat("%d");
    _chart->addAxis(_yAxis, Qt::AlignLeft);

    _chart->setTitle("Gradient Graph");

    _chartView = new QChartView(_chart);
    _chartView->setRenderHint(QPainter::Antialiasing);

    QVBoxLayout* layout = new QVBoxLayout();
    layout->addWidget(_chartView);
    setLayout(layout);
}

void GradientGraph::setNumDimensions(int numDimensions)
{
    _numDimensions = numDimensions;

    _seriesArray.resize(numDimensions);

    // Add numDimensions lineseries to the chart
    _chart->removeAllSeries();
    for (int d = 0; d < _numDimensions; d++)
    {
        // Create series
        _seriesArray[d] = new LineSeries(nullptr, d);
        _seriesArray[d]->setUseOpenGL(); // Speed up rendering

        connect(_seriesArray[d], &LineSeries::lineClicked, this, &GradientGraph::onLineClicked);

        // Add series to chart
        _chart->addSeries(_seriesArray[d]);

        // Attach chart axes to each series
        _seriesArray[d]->attachAxis(_xAxis);
        _seriesArray[d]->attachAxis(_yAxis);
    }
}

void GradientGraph::setValues(const std::vector<std::vector<float>>& values)
{
    // Compute the maximum value in the data and store values in a QList<QPointF>
    float maxValue = 0;
    std::vector<QList<QPointF>> pointLists(values.size());
    for (int d = 0; d < values.size(); d++)
    {
        float maxVal = *std::max_element(values[d].begin(), values[d].end());
        maxValue = maxVal > maxValue ? maxVal : maxValue;

        for (int i = 0; i < values[d].size(); i++)
        {
            pointLists[d].append(QPointF(i, values[d][i]));
        }
    }

    // Replace series values with new values
    for (int d = 0; d < values.size(); d++)
        _seriesArray[d]->replace(pointLists[d]);

    // Update axes to proper values
    _xAxis->setRange(0, values[0].size());
    _yAxis->setRange(0, maxValue);
    
    _chartView->update();
}

void GradientGraph::onLineClicked(int dim)
{
    emit lineClicked(dim);
}
