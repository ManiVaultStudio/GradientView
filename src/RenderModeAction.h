#pragma once

#include "PluginAction.h"

#include <QActionGroup>
#include <QHBoxLayout>

using namespace hdps::gui;

class QMenu;

class RenderModeAction : public PluginAction
{
protected: // Widget

    class Widget : public WidgetActionWidget {
    public:
        Widget(QWidget* parent, RenderModeAction* renderModeAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this, widgetFlags);
    };

public:
    RenderModeAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

protected:
    ToggleAction    _scatterPlotAction;
    ToggleAction    _densityPlotAction;
    ToggleAction    _contourPlotAction;
    QActionGroup    _actionGroup;

    friend class Widget;
};