#include "DensityPlotAction.h"
#include "Application.h"

#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"

using namespace hdps::gui;

DensityPlotAction::DensityPlotAction(ScatterplotPlugin* scatterplotPlugin) :
    PluginAction(scatterplotPlugin, "Density"),
    _sigmaAction(this, "Sigma", 0.01, 0.5, DEFAULT_SIGMA, DEFAULT_SIGMA, 3),
    _continuousUpdatesAction(this, "Live Updates", DEFAULT_CONTINUOUS_UPDATES, DEFAULT_CONTINUOUS_UPDATES)
{
    setToolTip("Density plot settings");

    _scatterplotPlugin->getWidget().addAction(&_sigmaAction);
    _scatterplotPlugin->getWidget().addAction(&_continuousUpdatesAction);

    const auto computeDensity = [this]() -> void {
        getScatterplotWidget().setSigma(_sigmaAction.getValue());

        const auto maxDensity = getScatterplotWidget().getDensityRenderer().getMaxDensity();

        if (maxDensity > 0)
            _scatterplotPlugin->getSettingsAction().getColoringAction().getColorMapAction().getSettingsAction().getHorizontalAxisAction().getRangeAction().setRange(0.0f, maxDensity);
    };

    connect(&_sigmaAction, &DecimalAction::valueChanged, this, [this, computeDensity](const double& value) {
        computeDensity();
    });

    const auto updateSigmaAction = [this]() {
        _sigmaAction.setUpdateDuringDrag(_continuousUpdatesAction.isChecked());
    };

    connect(&_continuousUpdatesAction, &ToggleAction::toggled, updateSigmaAction);

    connect(&_scatterplotPlugin->getPositionDataset(), &Dataset<Points>::changed, this, [this, updateSigmaAction, computeDensity](DatasetImpl* dataset) {
        updateSigmaAction();
        computeDensity();
    });

    connect(&getScatterplotWidget(), &ScatterplotWidget::renderModeChanged, this, [this, computeDensity](const ScatterplotWidget::RenderMode& renderMode) {
        computeDensity();
    });

    updateSigmaAction();
    computeDensity();
}

QMenu* DensityPlotAction::getContextMenu()
{
    auto menu = new QMenu("Plot settings");

    const auto renderMode = getScatterplotWidget().getRenderMode();

    const auto addActionToMenu = [menu](QAction* action) {
        auto actionMenu = new QMenu(action->text());

        actionMenu->addAction(action);

        menu->addMenu(actionMenu);
    };

    addActionToMenu(&_sigmaAction);
    addActionToMenu(&_continuousUpdatesAction);

    return menu;
}

DensityPlotAction::Widget::Widget(QWidget* parent, DensityPlotAction* densityPlotAction) :
    WidgetActionWidget(parent, densityPlotAction)
{
    auto layout = new QHBoxLayout();

    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(densityPlotAction->_sigmaAction.createLabelWidget(this));
    layout->addWidget(densityPlotAction->_sigmaAction.createWidget(this));
    layout->addWidget(densityPlotAction->_continuousUpdatesAction.createWidget(this));

    setLayout(layout);
}