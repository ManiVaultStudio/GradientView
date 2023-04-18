#include "ScalarSourceAction.h"
#include "Application.h"

#include "ScatterplotPlugin.h"
#include "PointData/PointData.h"

#include <QGridLayout>
#include <QHBoxLayout>

using namespace hdps::gui;

ScalarSourceAction::ScalarSourceAction(ScatterplotPlugin* scatterplotPlugin, const QString& title) :
    PluginAction(scatterplotPlugin, title),
    _model(this),
    _pickerAction(this, "Source"),
    _dimensionPickerAction(this, "Data dimension"),
    _offsetAction(this, "Offset", 0.0f, 100.0f, 0.0f, 0.0f, 2),
    _rangeAction(this, "Scalar range")
{
    _scatterplotPlugin->getWidget().addAction(&_pickerAction);
    _scatterplotPlugin->getWidget().addAction(&_dimensionPickerAction);

    // Configure scalar option picker
    _pickerAction.setCustomModel(&_model);
    _pickerAction.setToolTip("Pick scalar option");

    _offsetAction.setSuffix("px");

    // Invoked when the scalar source changed
    const auto scalarSourceChanged = [this]() -> void {

        // Current source index
        const auto sourceIndex = _pickerAction.getCurrentIndex();

        // Enable the dimension picker for datasets
        _dimensionPickerAction.setEnabled(sourceIndex >= ScalarSourceModel::DefaultRow::DatasetStart);

        // Assign the icon
        setIcon(_model.data(_model.index(sourceIndex, 0), Qt::DecorationRole).value<QIcon>());

        QSignalBlocker dimensionPickerActionBlocker(&_dimensionPickerAction);

        // Update dimension picker
        _dimensionPickerAction.setPointsDataset(Dataset<Points>(_model.getDataset(_pickerAction.getCurrentIndex())));

        // Disable offset action in constant mode
        _offsetAction.setEnabled(sourceIndex > ScalarSourceModel::DefaultRow::Constant);

        // Update the scalar range if possible
        updateScalarRange();
    };

    // Handle when the source index changes
    connect(&_pickerAction, &OptionAction::currentIndexChanged, this, scalarSourceChanged);

    // Update scalar range when dimension is picked
    connect(&_dimensionPickerAction, &DimensionPickerAction::currentDimensionIndexChanged, this, &ScalarSourceAction::updateScalarRange);

    // Notify others that the range changed when the user changes the range minimum
    connect(&_rangeAction.getRangeMinAction(), &DecimalAction::valueChanged, this, [this]() {
        emit scalarRangeChanged(_rangeAction.getRangeMinAction().getValue(), _rangeAction.getRangeMaxAction().getValue());
    });

    // Notify others that the range changed when the user changes the range maximum
    connect(&_rangeAction.getRangeMaxAction(), &DecimalAction::valueChanged, this, [this]() {
        emit scalarRangeChanged(_rangeAction.getRangeMinAction().getValue(), _rangeAction.getRangeMaxAction().getValue());
    });

    // Force initial update
    scalarSourceChanged();
}

ScalarSourceModel& ScalarSourceAction::getModel()
{
    return _model;
}

void ScalarSourceAction::updateScalarRange()
{
    // Get smart pointer to points dataset
    auto points = Dataset<Points>(_model.getDataset(_pickerAction.getCurrentIndex()));

    float minimum = 0.0f;
    float maximum = 1.0f;

    // Establish whether there is valid scalar range
    const auto hasScalarRange = points.isValid() && _pickerAction.getCurrentIndex() >= 1;

    if (hasScalarRange) {

        // Initialize minimum and maximum with extremes
        minimum = std::numeric_limits<float>::max();
        maximum = std::numeric_limits<float>::lowest();

        points->visitData([this, &minimum, &maximum, points](auto pointData) {

            // Get current dimension index
            const auto currentDimensionIndex = _dimensionPickerAction.getCurrentDimensionIndex();

            // Loop over all points to find minimum and maximum
            for (std::uint32_t pointIndex = 0; pointIndex < points->getNumPoints(); pointIndex++) {

                // Get point value for dimension
                auto pointValue = static_cast<float>(pointData[pointIndex][currentDimensionIndex]);
                
                // Adjust minimum and maximum
                minimum = std::min(minimum, pointValue);
                maximum = std::max(maximum, pointValue);
            }
        });
    }

    // Update range action
    _rangeAction.getRangeMinAction().initialize(minimum, maximum, minimum, 1);
    _rangeAction.getRangeMaxAction().initialize(minimum, maximum, maximum, 1);

    // Enable/disable range action
    _rangeAction.getRangeMinAction().setEnabled(hasScalarRange);
    _rangeAction.getRangeMaxAction().setEnabled(hasScalarRange);

    // Notify others the range ranged
    emit scalarRangeChanged(minimum, maximum);
}

ScalarSourceAction::Widget::Widget(QWidget* parent, ScalarSourceAction* scalarSourceAction) :
    WidgetActionWidget(parent, scalarSourceAction)
{
    //setFixedWidth(450);

    auto layout = new QGridLayout();

    // Create action widgets
    auto pickerLabelWidget          = scalarSourceAction->getPickerAction().createLabelWidget(this);
    auto pickerWidget               = scalarSourceAction->getPickerAction().createWidget(this);
    auto dimensionPickerLabelWidget = scalarSourceAction->getDimensionPickerAction().createLabelWidget(this);
    auto dimensionPickerWidget      = scalarSourceAction->getDimensionPickerAction().createWidget(this);
    auto offsetLabelWidget          = scalarSourceAction->getOffsetAction().createLabelWidget(this);
    auto offsetSpinBoxWidget        = scalarSourceAction->getOffsetAction().createWidget(this, DecimalAction::SpinBox);
    auto offsetSliderWidget         = scalarSourceAction->getOffsetAction().createWidget(this, DecimalAction::Slider);
    auto rangeMinLabelWidget        = scalarSourceAction->getRangeAction().getRangeMinAction().createLabelWidget(this);
    auto rangeMinSpinBoxWidget      = scalarSourceAction->getRangeAction().getRangeMinAction().createWidget(this, DecimalAction::SpinBox);
    auto rangeMinSliderWidget       = scalarSourceAction->getRangeAction().getRangeMinAction().createWidget(this, DecimalAction::Slider);
    auto rangeMaxLabelWidget        = scalarSourceAction->getRangeAction().getRangeMaxAction().createLabelWidget(this);
    auto rangeMaxSpinBoxWidget      = scalarSourceAction->getRangeAction().getRangeMaxAction().createWidget(this, DecimalAction::SpinBox);
    auto rangeMaxSliderWidget       = scalarSourceAction->getRangeAction().getRangeMaxAction().createWidget(this, DecimalAction::Slider);

    // Adjust size of the combo boxes to the contents
    pickerWidget->findChild<QComboBox*>("ComboBox")->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    dimensionPickerWidget->findChild<QComboBox*>("ComboBox")->setSizeAdjustPolicy(QComboBox::AdjustToContents);

    // Add source picker widget
    layout->addWidget(pickerLabelWidget, 0, 0);
    layout->addWidget(pickerWidget, 0, 1, 1, 2);

    // Add dimension picker widget
    layout->addWidget(dimensionPickerLabelWidget, 1, 0);
    layout->addWidget(dimensionPickerWidget, 1, 1, 1, 2);

    // Add scalar offset widget
    layout->addWidget(offsetLabelWidget, 2, 0);
    layout->addWidget(offsetSpinBoxWidget, 2, 1);
    layout->addWidget(offsetSliderWidget, 2, 2);

    // Add data clamp range minimum widget
    layout->addWidget(rangeMinLabelWidget, 3, 0);
    layout->addWidget(rangeMinSpinBoxWidget, 3, 1);
    layout->addWidget(rangeMinSliderWidget, 3, 2);

    // Add data clamp range maximum widget
    layout->addWidget(rangeMaxLabelWidget, 4, 0);
    layout->addWidget(rangeMaxSpinBoxWidget, 4, 1);
    layout->addWidget(rangeMaxSliderWidget, 4, 2);

    setPopupLayout(layout);
}
