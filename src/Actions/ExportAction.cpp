#include "ExportAction.h"

#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"

#include <QMenu>
#include <QGroupBox>

using namespace hdps::gui;

ExportAction::ExportAction(ScatterplotPlugin* scatterplotPlugin) :
    PluginAction(scatterplotPlugin, "Export"),
    _exportRankingsAction(scatterplotPlugin, "Export rankings"),
    _exportFloodnodesAction(scatterplotPlugin, "Export floodnodes")
{
    setIcon(hdps::Application::getIconFont("FontAwesome").getIcon("file-export"));

    connect(&_exportRankingsAction, &TriggerAction::triggered, this, [scatterplotPlugin]() {
        scatterplotPlugin->exportRankings();
    });
    connect(&_exportFloodnodesAction, &TriggerAction::triggered, this, [scatterplotPlugin]() {
        scatterplotPlugin->exportFloodnodes();
    });
}

QMenu* ExportAction::getContextMenu()
{
    QMenu* menu = new QMenu("Export settings");

    const auto addActionToMenu = [menu](QAction* action) {
        auto actionMenu = new QMenu(action->text());

        actionMenu->addAction(action);

        menu->addMenu(actionMenu);
    };

    addActionToMenu(&_exportRankingsAction);
    addActionToMenu(&_exportFloodnodesAction);

    return menu;
}

ExportAction::Widget::Widget(QWidget* parent, ExportAction* exportAction, const std::int32_t& widgetFlags) :
    WidgetActionWidget(parent, exportAction, widgetFlags)
{
    setToolTip("Export settings");
    //setStyleSheet("QToolButton { width: 36px; height: 36px; qproperty-iconSize: 18px; }");

    // Add widgets
    if (widgetFlags & PopupLayout)
    {
        auto layout = new QGridLayout();

        layout->setContentsMargins(4, 4, 4, 4);

        layout->addWidget(exportAction->getExportRankingsAction().createLabelWidget(this), 0, 0);
        layout->addWidget(exportAction->getExportRankingsAction().createWidget(this), 0, 1);
        layout->addWidget(exportAction->getExportFloodnodesAction().createLabelWidget(this), 1, 0);
        layout->addWidget(exportAction->getExportFloodnodesAction().createWidget(this), 1, 1);

        auto mainLayout = new QVBoxLayout();

        mainLayout->setContentsMargins(4, 4, 4, 4);

        setLayout(mainLayout);

        auto groupBox = new QGroupBox("Export Settings");

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
