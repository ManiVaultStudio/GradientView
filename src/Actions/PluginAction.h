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

    /**
     * Create collapsed widget
     * @param parent Parent widget
     * @return Pointer to collapsed widget
     */
    QWidget* createCollapsedWidget(QWidget* parent);

protected:
    ScatterplotPlugin*  _scatterplotPlugin;
};
