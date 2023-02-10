#include "FilterAction.h"

#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"

#include <QMenu>
#include <QGroupBox>

FilterAction::FilterAction(ScatterplotPlugin* scatterplotPlugin) :
    PluginAction(scatterplotPlugin, "Filter"),
    _spatialPeakFilterAction(scatterplotPlugin, "Spatial Peak Filter"),
    _hdPeakFilterAction(scatterplotPlugin, "HD Peak Filter"),
    _restrictToFloodAction(scatterplotPlugin, "Restrict to flood nodes", true, true),
    _innerFilterSizeAction(scatterplotPlugin, "Inner Filter Radius", 1, 10, 2.5f, 2.5f),
    _outerFilterSizeAction(scatterplotPlugin, "Outer Filter Radius", 2, 20, 5, 5),
    _hdInnerFilterSizeAction(scatterplotPlugin, "HD Inner Filter Size", 1, 10, 5, 5)
    //_hdOuterFilterSizeAction(scatterplotPlugin, "HD Outer Filter Size", 2, 10, 10, 10)
{
    setIcon(hdps::Application::getIconFont("FontAwesome").getIcon("bullseye"));

    auto& spatialPeakFilter = scatterplotPlugin->getSpatialPeakFilter();
    auto& hdPeakFilter = scatterplotPlugin->getHDPeakFilter();

    connect(&_spatialPeakFilterAction, &TriggerAction::triggered, this, [scatterplotPlugin]() {
        scatterplotPlugin->setFilterType(filters::FilterType::SPATIAL_PEAK);
        scatterplotPlugin->setFilterLabelText("Spatial Peak Ranking");
        scatterplotPlugin->getScatterplotWidget().showFiltersCircles(true);
        scatterplotPlugin->getScatterplotWidget().update();
    });
    connect(&_hdPeakFilterAction, &TriggerAction::triggered, this, [scatterplotPlugin]() {
        scatterplotPlugin->setFilterType(filters::FilterType::HD_PEAK);
        scatterplotPlugin->setFilterLabelText("HD Peak Ranking");
        scatterplotPlugin->getScatterplotWidget().showFiltersCircles(false);
        scatterplotPlugin->getScatterplotWidget().update();
    });

    connect(&_innerFilterSizeAction, &DecimalAction::valueChanged, [scatterplotPlugin, &spatialPeakFilter](const float& value) {
        float projSize = scatterplotPlugin->getProjectionSize();
        spatialPeakFilter.setInnerFilterRadius(value * 0.01f);
        scatterplotPlugin->getScatterplotWidget().setFilterRadii(Vector2f(spatialPeakFilter.getInnerFilterRadius() * projSize, spatialPeakFilter.getOuterFilterRadius() * projSize));
        scatterplotPlugin->onPointSelection();
    });
    connect(&_outerFilterSizeAction, &DecimalAction::valueChanged, [scatterplotPlugin, &spatialPeakFilter](const float& value) {
        float projSize = scatterplotPlugin->getProjectionSize();
        spatialPeakFilter.setOuterFilterRadius(value * 0.01f);
        scatterplotPlugin->getScatterplotWidget().setFilterRadii(Vector2f(spatialPeakFilter.getInnerFilterRadius() * projSize, spatialPeakFilter.getOuterFilterRadius() * projSize));
        scatterplotPlugin->onPointSelection();
    });

    connect(&_hdInnerFilterSizeAction, &IntegralAction::valueChanged, [&hdPeakFilter](int value) { hdPeakFilter.setInnerFilterSize(value); });
    //connect(&_hdOuterFilterSizeAction, &IntegralAction::valueChanged, [&hdPeakFilter](int value) { hdPeakFilter.setOuterFilterSize(value); });
}

QMenu* FilterAction::getContextMenu()
{
    QMenu* menu = new QMenu("Filter settings");

    const auto addActionToMenu = [menu](QAction* action) {
        auto actionMenu = new QMenu(action->text());

        actionMenu->addAction(action);

        menu->addMenu(actionMenu);
    };

    addActionToMenu(&_innerFilterSizeAction);
    addActionToMenu(&_outerFilterSizeAction);

    return menu;
}

void FilterAction::setFloodSteps(int numSteps)
{
    if (_hdInnerFilterSizeAction.getValue() > numSteps - 1)
        _hdInnerFilterSizeAction.setValue(numSteps - 1);

    _hdInnerFilterSizeAction.setMaximum(numSteps-1);
}

FilterAction::Widget::Widget(QWidget* parent, FilterAction* filterAction, const std::int32_t& widgetFlags) :
    WidgetActionWidget(parent, filterAction, widgetFlags)
{
    setToolTip("Filter settings");
    //setStyleSheet("QToolButton { width: 36px; height: 36px; qproperty-iconSize: 18px; }");
    
    // Add widgets
    if (widgetFlags & PopupLayout)
    {
        auto layout = new QGridLayout();

        layout->setContentsMargins(4, 4, 4, 4);

        layout->addWidget(filterAction->getSpatialPeakFilterAction().createWidget(this), 0, 0);
        layout->addWidget(filterAction->getHDPeakFilterAction().createWidget(this), 0, 1);

        layout->addWidget(filterAction->getRestrictToFloodAction().createWidget(this), 1, 0);

        layout->addWidget(filterAction->getInnerFilterSizeAction().createLabelWidget(this), 2, 0);
        QWidget* w = filterAction->getInnerFilterSizeAction().createWidget(this);
        w->setMinimumWidth(200);
        layout->addWidget(w, 2, 1);

        layout->addWidget(filterAction->getOuterFilterSizeAction().createLabelWidget(this), 3, 0);
        layout->addWidget(filterAction->getOuterFilterSizeAction().createWidget(this), 3, 1);

        layout->addWidget(filterAction->getHDInnerFilterSizeAction().createLabelWidget(this), 4, 0);
        layout->addWidget(filterAction->getHDInnerFilterSizeAction().createWidget(this), 4, 1);

        //layout->addWidget(filterAction->getHDOuterFilterSizeAction().createLabelWidget(this), 5, 0);
        //layout->addWidget(filterAction->getHDOuterFilterSizeAction().createWidget(this), 5, 1);

        auto mainLayout = new QVBoxLayout();

        mainLayout->setContentsMargins(4, 4, 4, 4);

        setLayout(mainLayout);

        auto groupBox = new QGroupBox("Filter Settings");

        groupBox->setLayout(layout);
        groupBox->setCheckable(false);

        mainLayout->addWidget(groupBox);
    }
    else
    {
        auto layout = new QHBoxLayout();

        layout->setContentsMargins(0, 0, 0, 0);

        // Big widget to force all other actions to become pop-ups
        QWidget* w = new QWidget();
        w->setMinimumWidth(2000);
        layout->addWidget(w);
        //layout->addWidget(filterAction->getInnerFilterSizeAction().createLabelWidget(this));
        //layout->addWidget(filterAction->getInnerFilterSizeAction().createWidget(this));

        //layout->addWidget(filterAction->getOuterFilterSizeAction().createLabelWidget(this));
        //layout->addWidget(filterAction->getOuterFilterSizeAction().createWidget(this));

        //layout->addWidget(filterAction->getHDInnerFilterSizeAction().createLabelWidget(this));
        //layout->addWidget(filterAction->getHDInnerFilterSizeAction().createWidget(this));

        //layout->addWidget(filterAction->getHDOuterFilterSizeAction().createLabelWidget(this));
        //layout->addWidget(filterAction->getHDOuterFilterSizeAction().createWidget(this));

        setLayout(layout);
    }
}
