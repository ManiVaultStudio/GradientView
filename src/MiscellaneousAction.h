#pragma once

#include "PluginAction.h"

#include <QActionGroup>

using namespace hdps::gui;

class QMenu;

class MiscellaneousAction : public PluginAction
{
protected: // Widget

    class Widget : public WidgetActionWidget {
    public:
        Widget(QWidget* parent, MiscellaneousAction* miscellaneousAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this, widgetFlags);
    };

public:
    MiscellaneousAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

protected:
    ColorAction  _backgroundColorAction;

    static const QColor DEFAULT_BACKGROUND_COLOR;

    friend class Widget;
};