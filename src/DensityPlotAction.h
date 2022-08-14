#pragma once

#include "PluginAction.h"

#include <QLabel>

using namespace hdps::gui;

class DensityPlotAction : public PluginAction
{
protected:
    class Widget : public WidgetActionWidget {
    public:
        Widget(QWidget* parent, DensityPlotAction* densityPlotAction);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this);
    };

public:
    DensityPlotAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

protected:
    DecimalAction       _sigmaAction;
    ToggleAction        _continuousUpdatesAction;

    static constexpr double DEFAULT_SIGMA = 0.15;
    static constexpr bool DEFAULT_CONTINUOUS_UPDATES = true;

    friend class Widget;
    friend class PlotAction;
};