#include "ScalarAction.h"
#include "ScalarSourceModel.h"

#include "ScatterplotPlugin.h"
#include "Application.h"

#include <QHBoxLayout>

using namespace hdps::gui;

ScalarAction::ScalarAction(ScatterplotPlugin* scatterplotPlugin, const QString& title, const float& minimum, const float& maximum, const float& value, const float& defaultValue) :
    PluginAction(scatterplotPlugin, "Scalar"),
    _magnitudeAction(this, title, minimum, maximum, value, defaultValue),
    _sourceAction(scatterplotPlugin, title)
{
    setText(title);

    _scatterplotPlugin->getWidget().addAction(&_sourceAction);

    // Notify others when the source selection changes
    connect(&_sourceAction.getPickerAction(), &OptionAction::currentIndexChanged, this, [this](const std::uint32_t& currentIndex) {
        emit sourceSelectionChanged(currentIndex);
    });

    // Pass-through magnitude updates
    connect(&_magnitudeAction, &DecimalAction::valueChanged, this, [this](const float& value) {
        emit magnitudeChanged(value);
    });

    // Pass-through scalar range updates
    connect(&_sourceAction, &ScalarSourceAction::scalarRangeChanged, this, [this](const float& minimum, const float& maximum) {
        emit scalarRangeChanged(minimum, maximum);
    });

    // Pass-through scalar offset updates
    connect(&_sourceAction.getOffsetAction(), &DecimalAction::valueChanged, this, [this](const float& value) {
        emit offsetChanged(value);
    });
}

void ScalarAction::addDataset(const Dataset<DatasetImpl>& dataset)
{
    // Get reference to the point size source model
    auto& sourceModel = _sourceAction.getModel();

    // Add dataset to the list of candidate datasets
    sourceModel.addDataset(dataset);

    // Connect to the data changed signal so that we can update the scatter plot point size appropriately
    connect(&sourceModel.getDatasets().last(), &Dataset<DatasetImpl>::dataChanged, this, [this, dataset]() {

        // Get smart pointer to current dataset
        const auto currentDataset = getCurrentDataset();

        // Only proceed if we have a valid point size dataset
        if (!currentDataset.isValid())
            return;

        // To do: schedule more efficient
        _sourceAction.updateScalarRange();

        // Update scatter plot widget point size if the dataset matches
        if (currentDataset == dataset)
            emit sourceDataChanged(dataset);
    });

    // Connect to the data changed signal so that we can update the scatter plot point size appropriately
    connect(&_magnitudeAction, &DecimalAction::valueChanged, this, [this, dataset](const float& value) {
        emit magnitudeChanged(value);
    });
}

void ScalarAction::removeAllDatasets()
{
    _sourceAction.getModel().removeAllDatasets();
}

Dataset<DatasetImpl> ScalarAction::getCurrentDataset()
{
    // Get reference to the scalar source model
    auto& scalarSourceModel = _sourceAction.getModel();

    // Get current scalar source index
    const auto currentSourceIndex = _sourceAction.getPickerAction().getCurrentIndex();

    // Only proceed if we have a valid point size dataset row index
    if (currentSourceIndex < ScalarSourceModel::DefaultRow::DatasetStart)
        return Dataset<DatasetImpl>();

    return scalarSourceModel.getDataset(currentSourceIndex);
}

void ScalarAction::setCurrentDataset(const Dataset<DatasetImpl>& dataset)
{
    // Obtain row index of the dataset
    const auto datasetRowIndex = _sourceAction.getModel().rowIndex(dataset);

    // Set picker current index if the dataset was found
    if (datasetRowIndex >= 0)
        _sourceAction.getPickerAction().setCurrentIndex(datasetRowIndex);
}

void ScalarAction::setCurrentSourceIndex(bool sourceIndex)
{
    _sourceAction.getPickerAction().setCurrentIndex(sourceIndex);
}

bool ScalarAction::isSourceConstant() const
{
    return _sourceAction.getPickerAction().getCurrentIndex() == ScalarSourceModel::DefaultRow::Constant;
}

bool ScalarAction::isSourceSelection() const
{
    return _sourceAction.getPickerAction().getCurrentIndex() == ScalarSourceModel::DefaultRow::Selection;
}

bool ScalarAction::isSourceDataset() const
{
    return _sourceAction.getPickerAction().getCurrentIndex() >= ScalarSourceModel::DefaultRow::DatasetStart;
}

ScalarAction::Widget::Widget(QWidget* parent, ScalarAction* scalarAction) :
    WidgetActionWidget(parent, scalarAction)
{
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);

    auto layout = new QHBoxLayout();

    layout->setContentsMargins(0, 0, 0, 0);

    // Create action widgets
    auto magnitudeWidget        = scalarAction->getMagnitudeAction().createWidget(this);
    auto sourceWidget           = scalarAction->getSourceAction().createCollapsedWidget(this);

    // And add them to the layout
    layout->addWidget(magnitudeWidget);
    layout->addWidget(sourceWidget);

    // Adjust size of the combo boxes to the contents
    //pointSizeByWidget->findChild<QComboBox*>("ComboBox")->setSizeAdjustPolicy(QComboBox::AdjustToContents);

    setLayout(layout);

    connect(&scalarAction->getSourceAction().getPickerAction(), &OptionAction::currentIndexChanged, this, [this, layout]() {
        layout->update();
    });
}
