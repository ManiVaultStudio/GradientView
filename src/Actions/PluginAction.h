#pragma once

#include "actions/Actions.h"

class ScatterplotPlugin;
class ScatterplotWidget;
class ProjectionView;

class PluginAction : public hdps::gui::WidgetAction
{
public:
    PluginAction(ScatterplotPlugin* scatterplotPlugin, const QString& title);

    ScatterplotWidget& getScatterplotWidget();
    std::vector<ProjectionView*>& getProjectionViews();

protected:
    ScatterplotPlugin*  _scatterplotPlugin;
};
