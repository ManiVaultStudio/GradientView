#include "ExportAction.h"

#include "GradientExplorerPlugin.h"
#include "Widgets/MainView.h"

#include <QMenu>
#include <QGroupBox>

using namespace mv::gui;

ExportAction::ExportAction(QObject* parent, const QString& title) :
    WidgetAction(parent, "Export Settings"),
    _exportRankingsAction(this, "Export rankings"),
    _exportFloodnodesAction(this, "Export floodnodes"),
    _importKnnGraphAction(this, "Import KNN Graph")
{
    setIconByName("file-export");

    setConfigurationFlag(WidgetAction::ConfigurationFlag::ForceCollapsedInGroup);
}

void ExportAction::initialize(GradientExplorerPlugin* scatterplotPlugin)
{
    connect(&_exportRankingsAction, &TriggerAction::triggered, this, [scatterplotPlugin]() {
        scatterplotPlugin->exportDimRankings();
    });
    connect(&_exportFloodnodesAction, &TriggerAction::triggered, this, [scatterplotPlugin]() {
        scatterplotPlugin->exportFloodnodes();
    });
    connect(&_importKnnGraphAction, &TriggerAction::triggered, this, [scatterplotPlugin]() {
        scatterplotPlugin->importKnnGraph();
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
    addActionToMenu(&_importKnnGraphAction);

    return menu;
}

ExportAction::Widget::Widget(QWidget* parent, ExportAction* exportAction, const std::int32_t& widgetFlags) :
    WidgetActionWidget(parent, exportAction, widgetFlags)
{
    setToolTip("Export settings");
    //setStyleSheet("QToolButton { width: 36px; height: 36px; qproperty-iconSize: 18px; }");

    auto layout = new QGridLayout();

    layout->setContentsMargins(4, 4, 4, 4);

    layout->addWidget(exportAction->getExportRankingsAction().createLabelWidget(this), 0, 0);
    layout->addWidget(exportAction->getExportRankingsAction().createWidget(this), 0, 1);
    layout->addWidget(exportAction->getExportFloodnodesAction().createLabelWidget(this), 1, 0);
    layout->addWidget(exportAction->getExportFloodnodesAction().createWidget(this), 1, 1);
    layout->addWidget(exportAction->getImportKnnGraphAction().createLabelWidget(this), 2, 0);
    layout->addWidget(exportAction->getImportKnnGraphAction().createWidget(this), 2, 1);

    setLayout(layout);
}
