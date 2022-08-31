#include "SelectionAction.h"
#include "Application.h"
#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"

#include "util/PixelSelectionTool.h"

#include <QHBoxLayout>
#include <QPushButton>

using namespace hdps::gui;

const auto allowedPixelSelectionTypes = PixelSelectionTypes({
    PixelSelectionType::Rectangle,
    PixelSelectionType::Brush,
    PixelSelectionType::Lasso,
    PixelSelectionType::Polygon
});

SelectionAction::SelectionAction(ScatterplotPlugin& scatterplotPlugin) :
    PixelSelectionAction(&scatterplotPlugin, &scatterplotPlugin.getScatterplotWidget(), scatterplotPlugin.getScatterplotWidget().getPixelSelectionTool(), allowedPixelSelectionTypes),
    _scatterplotPlugin(scatterplotPlugin)
{
    setIcon(hdps::Application::getIconFont("FontAwesome").getIcon("mouse-pointer"));
    
    connect(&getSelectAllAction(), &QAction::triggered, [this]() {
        if (_scatterplotPlugin.getPositionDataset().isValid())
            _scatterplotPlugin.getPositionDataset()->selectAll();
    });

    connect(&getClearSelectionAction(), &QAction::triggered, [this]() {
        if (_scatterplotPlugin.getPositionDataset().isValid())
            _scatterplotPlugin.getPositionDataset()->selectNone();
    });

    connect(&getInvertSelectionAction(), &QAction::triggered, [this]() {
        if (_scatterplotPlugin.getPositionDataset().isValid())
            _scatterplotPlugin.getPositionDataset()->selectInvert();
    });
}

SelectionAction::Widget::Widget(QWidget* parent, SelectionAction* selectionAction, const std::int32_t& widgetFlags) :
    WidgetActionWidget(parent, selectionAction, widgetFlags)
{
    auto typeWidget                     = selectionAction->getTypeAction().createWidget(this);
    auto brushRadiusWidget              = selectionAction->getBrushRadiusAction().createWidget(this);
    auto modifierAddWidget              = selectionAction->getModifierAddAction().createWidget(this, ToggleAction::PushButtonIcon);
    auto modifierSubtractWidget         = selectionAction->getModifierSubtractAction().createWidget(this, ToggleAction::PushButtonIcon);
    auto clearSelectionWidget           = selectionAction->getClearSelectionAction().createWidget(this);
    auto selectAllWidget                = selectionAction->getSelectAllAction().createWidget(this);
    auto invertSelectionWidget          = selectionAction->getInvertSelectionAction().createWidget(this);
    auto notifyDuringSelectionWidget    = selectionAction->getNotifyDuringSelectionAction().createWidget(this);

    if (widgetFlags & PopupLayout) {
        const auto getTypeWidget = [&, this]() -> QWidget* {
            auto layout = new QHBoxLayout();

            layout->setContentsMargins(0, 0, 0, 0);
            layout->addWidget(typeWidget);
            layout->addWidget(modifierAddWidget);
            layout->addWidget(modifierSubtractWidget);
            layout->itemAt(0)->widget()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

            auto widget = new QWidget();

            widget->setLayout(layout);

            return widget;
        };

        const auto getSelectWidget = [&, this]() -> QWidget* {
            auto layout = new QHBoxLayout();

            layout->setContentsMargins(0, 0, 0, 0);
            layout->addWidget(clearSelectionWidget);
            layout->addWidget(selectAllWidget);
            layout->addWidget(invertSelectionWidget);
            layout->addStretch(1);

            auto widget = new QWidget();

            widget->setLayout(layout);

            return widget;
        };

        auto layout = new QGridLayout();

        layout->addWidget(selectionAction->getTypeAction().createLabelWidget(this), 0, 0);
        layout->addWidget(getTypeWidget(), 0, 1);
        layout->addWidget(selectionAction->_brushRadiusAction.createLabelWidget(this), 1, 0);
        layout->addWidget(brushRadiusWidget, 1, 1);
        layout->addWidget(getSelectWidget(), 2, 1);
        layout->addWidget(notifyDuringSelectionWidget, 3, 1);
        layout->itemAtPosition(1, 1)->widget()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

        setPopupLayout(layout);
    }
    else {
        auto layout = new QHBoxLayout();

        layout->setContentsMargins(0, 0, 0, 0);
        layout->addWidget(typeWidget);
        layout->addWidget(brushRadiusWidget);
        layout->addWidget(modifierAddWidget);
        layout->addWidget(modifierSubtractWidget);
        layout->addWidget(clearSelectionWidget);
        layout->addWidget(selectAllWidget);
        layout->addWidget(invertSelectionWidget);
        layout->addWidget(notifyDuringSelectionWidget);

        setLayout(layout);
    }
}
