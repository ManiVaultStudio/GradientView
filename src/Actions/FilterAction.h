#pragma once

#include "PluginAction.h"

class FilterAction : public PluginAction
{
protected: // Widget
    class Widget : public WidgetActionWidget
    {
    public:
        Widget(QWidget* parent, FilterAction* filterAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override
    {
        return new Widget(parent, this, widgetFlags);
    }

public:
    FilterAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

    /**
     *
     *
     */

public: // Action getters
    TriggerAction& getSpatialPeakFilterAction() { return _spatialPeakFilterAction; }
    TriggerAction& getHDPeakFilterAction() { return _hdPeakFilterAction; }

    DecimalAction& getInnerFilterSizeAction() { return _innerFilterSizeAction; }
    DecimalAction& getOuterFilterSizeAction() { return _outerFilterSizeAction; }

    IntegralAction& getHDInnerFilterSizeAction() { return _hdInnerFilterSizeAction; }
    IntegralAction& getHDOuterFilterSizeAction() { return _hdOuterFilterSizeAction; }

protected:
    TriggerAction       _spatialPeakFilterAction;
    TriggerAction       _hdPeakFilterAction;

    DecimalAction       _innerFilterSizeAction;
    DecimalAction       _outerFilterSizeAction;

    IntegralAction      _hdInnerFilterSizeAction;
    IntegralAction      _hdOuterFilterSizeAction;
};
