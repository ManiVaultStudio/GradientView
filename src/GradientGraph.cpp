#include "GradientGraph.h"

#include <QVBoxLayout>

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
    _series(new QLineSeries()),
    _chart(new QChart()),
    _chartView(nullptr)
{
    _series->append(0, 6);
    _series->append(2, 4);
    _series->append(3, 8);
    _series->append(7, 4);
    _series->append(10, 5);
    *_series << QPointF(11, 1) << QPointF(13, 3) << QPointF(17, 6) << QPointF(18, 3) << QPointF(20, 2);

    _chart->legend()->hide();
    _chart->addSeries(_series);
    _chart->createDefaultAxes();
    _chart->setTitle("Gradient Graph");

    _chartView = new QChartView(_chart);
    _chartView->setRenderHint(QPainter::Antialiasing);

    QVBoxLayout* layout = new QVBoxLayout();
    layout->addWidget(_chartView);
    setLayout(layout);
}

void GradientGraph::setDimension(const DataMatrix& data, int dim)
{
    auto dimValues = data.col(dim);

    _chart->removeAllSeries();
    _series = new QLineSeries();
    for (int i = 0; i < dimValues.size(); i++)
    {
        _series->append(i, dimValues(i));
    }
    _chart->addSeries(_series);
    _chart->axisX()->setRange(0, dimValues.size());

    _chartView->update();
}

void GradientGraph::setValues(const std::vector<std::vector<float>>& values)
{
    _chart->removeAllSeries();

    for (int d = 0; d < values.size(); d++)
    {
        LineSeries* series = new LineSeries(nullptr, d);
        for (int i = 0; i < values[d].size(); i++)
        {
            series->append(i, values[d][i]);
        }
        connect(series, &LineSeries::lineClicked, this, &GradientGraph::onLineClicked);
        _chart->addSeries(series);
    }

    _chart->axisX()->setRange(0, values[0].size());
    _chart->axisY()->setRange(0, 100);

    _chartView->update();
}

void GradientGraph::onLineClicked(int dim)
{
    emit lineClicked(dim);
}
