#pragma once

#include "PluginAction.h"

#include "ScalarSourceModel.h"

#include <DimensionPickerAction.h>

using namespace hdps::gui;

/**
 * Scalar source action class
 *
 * Action class picking a scalar source (scalar by constant or dataset dimension)
 *
 * @author Thomas Kroes
 */
class ScalarSourceAction : public PluginAction
{

    Q_OBJECT

protected: // Widget

    /** Widget class for scalar source action */
    class Widget : public WidgetActionWidget {
    public:

        /**
         * Constructor
         * @param parent Pointer to parent widget
         * @param scalarSourceAction Pointer to scalar source action
         */
        Widget(QWidget* parent, ScalarSourceAction* scalarSourceAction);
    };

protected:

    /**
     * Get widget representation of the scalar source action
     * @param parent Pointer to parent widget
     * @param widgetFlags Widget flags for the configuration of the widget (type)
     */
    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this);
    };

public:

    /**
     * Constructor
     * @param scatterplotPlugin Pointer to scatter plot plugin
     * @param title Title
     */
    ScalarSourceAction(ScatterplotPlugin* scatterplotPlugin, const QString& title);

    /** Get the scalar source model */
    ScalarSourceModel& getModel();

    /** Update scalar range */
    void updateScalarRange();

public: // Action getters

    OptionAction& getPickerAction() { return _pickerAction; }
    DimensionPickerAction& getDimensionPickerAction() { return _dimensionPickerAction; }
    DecimalAction& getOffsetAction() { return _offsetAction; }
    DecimalRangeAction& getRangeAction() { return _rangeAction; }

    const OptionAction& getPickerAction() const { return _pickerAction; }
    const DimensionPickerAction& getDimensionPickerAction() const { return _dimensionPickerAction; }
    const DecimalAction& getOffsetAction() const { return _offsetAction; }
    const DecimalRangeAction& getRangeAction() const { return _rangeAction; }

signals:

    /**
     * Signals that the scalar range changed
     * @param minimum Scalar range minimum
     * @param maximum Scalar range maximum
     */
    void scalarRangeChanged(const float& minimum, const float& maximum);

protected:
    ScalarSourceModel       _model;                     /** Scalar model */
    OptionAction            _pickerAction;              /** Source picker action */
    DimensionPickerAction   _dimensionPickerAction;     /** Dimension picker action */
    DecimalAction           _offsetAction;              /** Scalar source offset action */
    DecimalRangeAction      _rangeAction;               /** Range action */
};