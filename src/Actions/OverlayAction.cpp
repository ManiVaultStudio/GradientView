#include "OverlayAction.h"

#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"

#include <QMenu>
#include <QGroupBox>

OverlayAction::OverlayAction(ScatterplotPlugin* scatterplotPlugin) :
    PluginAction(scatterplotPlugin, "Filter"),
    _floodDecimal(this, "Flood nodes", 10, 500, 10, 10),
    _floodStepsAction(this, "Flood steps", 3, 50, 10, 10),
    _sharedDistAction(this, "Shared distances", false, false)
{
    setIcon(hdps::Application::getIconFont("FontAwesome").getIcon("image"));

    connect(&_floodDecimal, &IntegralAction::valueChanged, this, [scatterplotPlugin](int32_t value)
    {
        scatterplotPlugin->rebuildKnnGraph(value);
    });

    connect(&_floodStepsAction, &IntegralAction::valueChanged, this, [scatterplotPlugin](int32_t value)
    {
        scatterplotPlugin->setFloodSteps(value);
        scatterplotPlugin->onPointSelection();
    });

    _triggers << TriggersAction::Trigger("Flood Steps", "Color flood points by closeness to seed point in HD space");
    _triggers << TriggersAction::Trigger("Top Dimension Values", "Color flood points by values of top ranked dimension");
    _triggers << TriggersAction::Trigger("Local Dimensionality", "Color flood points by local intrinsic dimensionality");
    _triggers << TriggersAction::Trigger("Directions", "Show major eigenvector directions over flood points");

    auto overlayGroupAction = new GroupAction(&scatterplotPlugin->getWidget(), true);
    overlayGroupAction->setText("Flood Nodes Overlay");
    overlayGroupAction->setShowLabels(false);

    TriggersAction* overlayTriggers = new TriggersAction(overlayGroupAction, "Overlay Triggers", _triggers);

    connect(overlayTriggers, &TriggersAction::triggered, this, [scatterplotPlugin](int32_t triggerIndex)
    {
        scatterplotPlugin->getScatterplotWidget().showDirections(false);
        switch (triggerIndex)
        {
        case 0: scatterplotPlugin->setOverlayType(OverlayType::NONE); break;
        case 1: scatterplotPlugin->setOverlayType(OverlayType::DIM_VALUES); break;
        case 2: scatterplotPlugin->setOverlayType(OverlayType::LOCAL_DIMENSIONALITY); break;
        case 3: {scatterplotPlugin->setOverlayType(OverlayType::DIRECTIONS); scatterplotPlugin->getScatterplotWidget().showDirections(true); break; }
        }
    });
}

QMenu* OverlayAction::getContextMenu()
{
    QMenu* menu = new QMenu("Overlay settings");

    const auto addActionToMenu = [menu](QAction* action) {
        auto actionMenu = new QMenu(action->text());

        actionMenu->addAction(action);

        menu->addMenu(actionMenu);
    };

    addActionToMenu(&_floodDecimal);
    addActionToMenu(&_floodStepsAction);
    addActionToMenu(&_sharedDistAction);

    return menu;
}

OverlayAction::Widget::Widget(QWidget* parent, OverlayAction* overlayAction, const std::int32_t& widgetFlags) :
    WidgetActionWidget(parent, overlayAction, widgetFlags)
{
    setToolTip("Overlay settings");
    //setStyleSheet("QToolButton { width: 36px; height: 36px; qproperty-iconSize: 18px; }");

    // Add widgets
    if (widgetFlags & PopupLayout)
    {
        auto layout = new QGridLayout();

        layout->setContentsMargins(4, 4, 4, 4);

        layout->addWidget(overlayAction->getFloodDecimalAction().createLabelWidget(this), 0, 0);
        layout->addWidget(overlayAction->getFloodDecimalAction().createWidget(this), 0, 1);

        layout->addWidget(overlayAction->getFloodStepsAction().createLabelWidget(this));
        layout->addWidget(overlayAction->getFloodStepsAction().createWidget(this));

        layout->addWidget(overlayAction->getSharedDistAction().createLabelWidget(this));
        layout->addWidget(overlayAction->getSharedDistAction().createWidget(this));


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
