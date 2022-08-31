#include "PositionAction.h"
#include "Application.h"

#include "ScatterplotPlugin.h"

#include <QMenu>
#include <QComboBox>

using namespace hdps::gui;

PositionAction::PositionAction(ScatterplotPlugin* scatterplotPlugin) :
    PluginAction(scatterplotPlugin, "Position"),
    _xDimensionPickerAction(this, "X"),
    _yDimensionPickerAction(this, "Y")
{
    setIcon(hdps::Application::getIconFont("FontAwesome").getIcon("ruler-combined"));

    // Add actions to scatter plot plugin (for shortcuts)
    _scatterplotPlugin->getWidget().addAction(&_xDimensionPickerAction);
    _scatterplotPlugin->getWidget().addAction(&_yDimensionPickerAction);

    // Set tooltips
    _xDimensionPickerAction.setToolTip("X dimension");
    _yDimensionPickerAction.setToolTip("Y dimension");

    // Update scatter plot when the x-dimension changes
    connect(&_xDimensionPickerAction, &DimensionPickerAction::currentDimensionIndexChanged, [this, scatterplotPlugin](const std::uint32_t& currentDimensionIndex) {
        scatterplotPlugin->setXDimension(currentDimensionIndex);
    });

    // Update scatter plot when the y-dimension changes
    connect(&_yDimensionPickerAction, &DimensionPickerAction::currentDimensionIndexChanged, [this, scatterplotPlugin](const std::uint32_t& currentDimensionIndex) {
        scatterplotPlugin->setYDimension(currentDimensionIndex);
    });

    // Set dimension defaults when the position dataset changes
    connect(&scatterplotPlugin->getPositionDataset(), &Dataset<Points>::changed, this, [this]() {

        // Assign position dataset to x- and y dimension action
        _xDimensionPickerAction.setPointsDataset(_scatterplotPlugin->getPositionDataset());
        _yDimensionPickerAction.setPointsDataset(_scatterplotPlugin->getPositionDataset());

        // Assign current and default index to x-dimension action
        _xDimensionPickerAction.setCurrentDimensionIndex(0);
        _xDimensionPickerAction.setDefaultDimensionIndex(0);

        // Establish y-dimension
        const auto yIndex = _xDimensionPickerAction.getNumberOfDimensions() >= 2 ? 1 : 0;

        // Assign current and default index to y-dimension action
        _yDimensionPickerAction.setCurrentDimensionIndex(yIndex);
        _yDimensionPickerAction.setDefaultDimensionIndex(yIndex);
    });
}

QMenu* PositionAction::getContextMenu(QWidget* parent /*= nullptr*/)
{
    auto menu = new QMenu("Position", parent);

    auto xDimensionMenu = new QMenu("X dimension");
    auto yDimensionMenu = new QMenu("Y dimension");

    xDimensionMenu->addAction(&_xDimensionPickerAction);
    yDimensionMenu->addAction(&_yDimensionPickerAction);

    menu->addMenu(xDimensionMenu);
    menu->addMenu(yDimensionMenu);

    return menu;
}

std::int32_t PositionAction::getDimensionX() const
{
    return _xDimensionPickerAction.getCurrentDimensionIndex();
}

std::int32_t PositionAction::getDimensionY() const
{
    return _yDimensionPickerAction.getCurrentDimensionIndex();
}

PositionAction::Widget::Widget(QWidget* parent, PositionAction* positionAction, const std::int32_t& widgetFlags) :
    WidgetActionWidget(parent, positionAction, widgetFlags)
{
    auto xDimensionLabel    = positionAction->_xDimensionPickerAction.createLabelWidget(this);
    auto yDimensionLabel    = positionAction->_yDimensionPickerAction.createLabelWidget(this);
    auto xDimensionWidget   = positionAction->_xDimensionPickerAction.createWidget(this);
    auto yDimensionWidget   = positionAction->_yDimensionPickerAction.createWidget(this);

    xDimensionWidget->findChild<QComboBox*>("ComboBox")->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    yDimensionWidget->findChild<QComboBox*>("ComboBox")->setSizeAdjustPolicy(QComboBox::AdjustToContents);

    if (widgetFlags & PopupLayout) {
        auto layout = new QGridLayout();

        layout->addWidget(xDimensionLabel, 0, 0);
        layout->addWidget(xDimensionWidget, 0, 1);
        layout->addWidget(yDimensionLabel, 1, 0);
        layout->addWidget(yDimensionWidget, 1, 1);

        setPopupLayout(layout);
    }
    else {
        auto layout = new QHBoxLayout();

        layout->setContentsMargins(0, 0, 0, 0);
        layout->addWidget(xDimensionLabel);
        layout->addWidget(xDimensionWidget);
        layout->addWidget(yDimensionLabel);
        layout->addWidget(yDimensionWidget);

        setLayout(layout);
    }
}
