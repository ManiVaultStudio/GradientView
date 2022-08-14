#pragma once

#include "PluginAction.h"
#include "PointPlotAction.h"
#include "DensityPlotAction.h"

using namespace hdps::gui;

class PlotAction : public PluginAction
{
protected: // Widget

    class Widget : public WidgetActionWidget {
    public:
        Widget(QWidget* parent, PlotAction* plotAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this, widgetFlags);
    };

public:
    PlotAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

public: // Action getters

    PointPlotAction& getPointPlotAction() { return _pointPlotAction; }
    DensityPlotAction& getDensityPlotAction() { return _densityPlotAction; }

protected:
    PointPlotAction     _pointPlotAction;
    DensityPlotAction   _densityPlotAction;

    friend class Widget;
};