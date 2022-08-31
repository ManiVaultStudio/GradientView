#pragma once

#include "PluginAction.h"

#include <DimensionPickerAction.h>

#include <QHBoxLayout>
#include <QLabel>

using namespace hdps::gui;

/**
 * Position action class
 *
 * Action class for picking data dimensions for the point positions
 *
 * @author Thomas Kroes
 */
class PositionAction : public PluginAction
{
protected: // Widget

    /** Widget class for position action */
    class Widget : public WidgetActionWidget {
    public:

        /**
         * Constructor
         * @param parent Pointer to parent widget
         * @param positionAction Pointer to position action
         */
        Widget(QWidget* parent, PositionAction* positionAction, const std::int32_t& widgetFlags);
    };

protected:

    /**
     * Get widget representation of the position action
     * @param parent Pointer to parent widget
     * @param widgetFlags Widget flags for the configuration of the widget (type)
     */
    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this, widgetFlags);
    };

public:

    /**
     * Constructor
     * @param scatterplotPlugin Pointer to scatter plot plugin
     */
    PositionAction(ScatterplotPlugin* scatterplotPlugin);

    /**
     * Get the context menu for the action
     * @param parent Parent widget
     * @return Context menu
     */
    QMenu* getContextMenu(QWidget* parent = nullptr) override;

    /** Get current x-dimension */
    std::int32_t getDimensionX() const;

    /** Get current y-dimension */
    std::int32_t getDimensionY() const;

protected:
    DimensionPickerAction    _xDimensionPickerAction;   /** X-dimension picker action */
    DimensionPickerAction    _yDimensionPickerAction;   /** Y-dimension picker action */

    friend class Widget;
};