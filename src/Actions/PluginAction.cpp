#include "PluginAction.h"

#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"
#include "ProjectionView.h"

using namespace hdps::gui;

PluginAction::PluginAction(ScatterplotPlugin* scatterplotPlugin, const QString& title) :
    WidgetAction(reinterpret_cast<QObject*>(scatterplotPlugin)),
    _scatterplotPlugin(scatterplotPlugin)
{
    _scatterplotPlugin->getWidget().addAction(this);

    setText(title);
    setToolTip(title);
}

ScatterplotWidget& PluginAction::getScatterplotWidget()
{
    Q_ASSERT(_scatterplotPlugin != nullptr);

    return _scatterplotPlugin->getScatterplotWidget();
}

std::vector<ProjectionView*>& PluginAction::getProjectionViews()
{
    Q_ASSERT(_scatterplotPlugin != nullptr);

    return _scatterplotPlugin->getProjectionViews();
}
