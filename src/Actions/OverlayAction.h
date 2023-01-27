#pragma once

#include "PluginAction.h"

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
    //DecimalAction& getInnerFilterSizeAction() { return _innerFilterSizeAction; }
    //DecimalAction& getOuterFilterSizeAction() { return _outerFilterSizeAction; }

    //IntegralAction& getHDInnerFilterSizeAction() { return _hdInnerFilterSizeAction; }
    //IntegralAction& getHDOuterFilterSizeAction() { return _hdOuterFilterSizeAction; }
    IntegralAction& getFloodDecimalAction() { return _floodDecimal; }
    IntegralAction& getFloodStepsAction() { return _floodStepsAction; }
    ToggleAction& getSharedDistAction() { return _sharedDistAction; }


protected:
    IntegralAction                      _floodDecimal;
    IntegralAction                      _floodStepsAction;
    ToggleAction                        _sharedDistAction;

    QVector<TriggersAction::Trigger>    _triggers;

    //DecimalAction       _innerFilterSizeAction;
    //DecimalAction       _outerFilterSizeAction;

    //IntegralAction      _hdInnerFilterSizeAction;
    //IntegralAction      _hdOuterFilterSizeAction;
};
