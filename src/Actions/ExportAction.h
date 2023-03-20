#pragma once

#include "PluginAction.h"

#include "actions/GroupAction.h"

using namespace hdps::gui;

class ExportAction : public PluginAction
{
protected: // Widget
    class Widget : public WidgetActionWidget
    {
    public:
        Widget(QWidget* parent, ExportAction* exportAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override
    {
        return new Widget(parent, this, widgetFlags);
    }

public:
    ExportAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

    /**
     *
     *
     */

public: // Action getters
    TriggerAction& getExportRankingsAction() { return _exportRankingsAction; }
    TriggerAction& getExportFloodnodesAction() { return _exportFloodnodesAction; }

protected:
    TriggerAction       _exportRankingsAction;
    TriggerAction       _exportFloodnodesAction;
};
