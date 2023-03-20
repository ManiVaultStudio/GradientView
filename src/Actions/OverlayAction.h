#pragma once

#include "PluginAction.h"

#include "actions/GroupAction.h"

using namespace hdps::gui;

class OverlayAction : public PluginAction
{
protected: // Widget
    class Widget : public WidgetActionWidget
    {
    public:
        Widget(QWidget* parent, OverlayAction* overlayAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override
    {
        return new Widget(parent, this, widgetFlags);
    }

public:
    OverlayAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

    /**
     *
     *
     */

public: // Action getters
    IntegralAction& getFloodDecimalAction() { return _floodDecimal; }
    IntegralAction& getFloodStepsAction() { return _floodStepsAction; }
    ToggleAction& getSharedDistAction() { return _sharedDistAction; }

    TriggerAction& getFloodOverlayAction() { return _floodOverlayAction; }
    TriggerAction& getDimensionOverlayAction() { return _dimensionOverlayAction; }
    TriggerAction& getDimensionalityOverlayAction() { return _dimensionalityOverlayAction; }

    //GroupAction& getOverlayGroupAction() { return _overlayGroupAction; }

protected:
    IntegralAction                      _floodDecimal;
    IntegralAction                      _floodStepsAction;
    ToggleAction                        _sharedDistAction;

    TriggerAction                       _floodOverlayAction;
    TriggerAction                       _dimensionOverlayAction;
    TriggerAction                       _dimensionalityOverlayAction;

    //QVector<TriggersAction::Trigger>    _triggers;
    //GroupAction _overlayGroupAction;
};
