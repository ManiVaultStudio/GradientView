#pragma once

#include "PluginAction.h"

using namespace hdps::gui;

class SubsetAction : public PluginAction
{
protected: // Widget

    class Widget : public WidgetActionWidget {
    public:
        Widget(QWidget* parent, SubsetAction* subsetAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this, widgetFlags);
    };

public:
    SubsetAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

protected:
    StringAction     _subsetNameAction;
    TriggerAction    _createSubsetAction;
    OptionAction     _sourceDataAction;

    friend class Widget;
};