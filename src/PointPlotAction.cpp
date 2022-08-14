#include "PointPlotAction.h"
#include "ScalarSourceAction.h"

#include "Application.h"

#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"
#include "DataHierarchyItem.h"

using namespace gui;

PointPlotAction::PointPlotAction(ScatterplotPlugin* scatterplotPlugin) :
    PluginAction(scatterplotPlugin, "Point"),
    _sizeAction(scatterplotPlugin, "Point size", 0.0, 100.0, DEFAULT_POINT_SIZE, DEFAULT_POINT_SIZE),
    _opacityAction(scatterplotPlugin, "Point opacity", 0.0, 100.0, DEFAULT_POINT_OPACITY, DEFAULT_POINT_OPACITY),
    _pointSizeScalars(),
    _pointOpacityScalars(),
    _focusSelection(this, "Focus selection"),
    _lastOpacitySourceIndex(-1)
{
    _scatterplotPlugin->getWidget().addAction(&_sizeAction);
    _scatterplotPlugin->getWidget().addAction(&_opacityAction);

    _sizeAction.getMagnitudeAction().setSuffix("px");
    _opacityAction.getMagnitudeAction().setSuffix("%");

    _opacityAction.getSourceAction().getOffsetAction().setSuffix("%");

    _focusSelection.setToolTip("Put focus on selected points by modulating the point opacity");

    // Update size by action when the position dataset changes
    connect(&_scatterplotPlugin->getPositionDataset(), &Dataset<Points>::changed, this, [this]() {

        // Get reference to position dataset
        const auto positionDataset = _scatterplotPlugin->getPositionDataset();

        // Do not update if no position dataset is loaded
        if (!positionDataset.isValid())
            return;

        // Remove all datasets from the models
        _sizeAction.removeAllDatasets();
        _opacityAction.removeAllDatasets();

        // Add the position dataset
        _sizeAction.addDataset(positionDataset);
        _opacityAction.addDataset(positionDataset);

        // Get smart pointer to position source dataset
        const auto positionSourceDataset = _scatterplotPlugin->getPositionSourceDataset();

        // Add source position dataset (if position dataset is derived)
        if (positionSourceDataset.isValid()) {

            // Add the position dataset
            _sizeAction.addDataset(positionSourceDataset);
            _opacityAction.addDataset(positionSourceDataset);
        }

        // Update the color by action
        updateDefaultDatasets();

        // Reset the point size and opacity scalars
        updateScatterPlotWidgetPointSizeScalars();
        updateScatterPlotWidgetPointOpacityScalars();

        // Reset
        _sizeAction.getSourceAction().getPickerAction().setCurrentIndex(0);
        _opacityAction.getSourceAction().getPickerAction().setCurrentIndex(0);
    });

    // Update default datasets when a child is added to or removed from the position dataset
    connect(&_scatterplotPlugin->getPositionDataset(), &Dataset<Points>::dataChildAdded, this, &PointPlotAction::updateDefaultDatasets);
    connect(&_scatterplotPlugin->getPositionDataset(), &Dataset<Points>::dataChildRemoved, this, &PointPlotAction::updateDefaultDatasets);

    // Update scatter plot widget point size scalars
    connect(&_sizeAction, &ScalarAction::magnitudeChanged, this, &PointPlotAction::updateScatterPlotWidgetPointSizeScalars);
    connect(&_sizeAction, &ScalarAction::offsetChanged, this, &PointPlotAction::updateScatterPlotWidgetPointSizeScalars);
    connect(&_sizeAction, &ScalarAction::sourceSelectionChanged, this, &PointPlotAction::updateScatterPlotWidgetPointSizeScalars);
    connect(&_sizeAction, &ScalarAction::sourceDataChanged, this, &PointPlotAction::updateScatterPlotWidgetPointSizeScalars);
    connect(&_sizeAction, &ScalarAction::scalarRangeChanged, this, &PointPlotAction::updateScatterPlotWidgetPointSizeScalars);

    // Update scatter plot widget point opacity scalars
    connect(&_opacityAction, &ScalarAction::magnitudeChanged, this, &PointPlotAction::updateScatterPlotWidgetPointOpacityScalars);
    connect(&_opacityAction, &ScalarAction::offsetChanged, this, &PointPlotAction::updateScatterPlotWidgetPointOpacityScalars);
    connect(&_opacityAction, &ScalarAction::sourceSelectionChanged, this, &PointPlotAction::updateScatterPlotWidgetPointOpacityScalars);
    connect(&_opacityAction, &ScalarAction::sourceDataChanged, this, &PointPlotAction::updateScatterPlotWidgetPointOpacityScalars);
    connect(&_opacityAction, &ScalarAction::scalarRangeChanged, this, &PointPlotAction::updateScatterPlotWidgetPointOpacityScalars);

    // Update the point size and opacity scalars when the selection of the position dataset changes
    connect(&_scatterplotPlugin->getPositionDataset(), &Dataset<Points>::dataSelectionChanged, this, &PointPlotAction::updateScatterPlotWidgetPointSizeScalars);
    connect(&_scatterplotPlugin->getPositionDataset(), &Dataset<Points>::dataSelectionChanged, this, &PointPlotAction::updateScatterPlotWidgetPointOpacityScalars);

    // For convenience, set the offset to double the magnitude in case of a selection source
    connect(&_sizeAction, &ScalarAction::sourceSelectionChanged, this, [this](const std::uint32_t& sourceSelectionIndex) {
        switch (sourceSelectionIndex)
        {
            case ScalarSourceModel::DefaultRow::Constant:
            {
                _sizeAction.getSourceAction().getOffsetAction().setValue(0.0f);

                break;
            }

            case ScalarSourceModel::DefaultRow::Selection:
            {
                _sizeAction.getSourceAction().getOffsetAction().setValue(_sizeAction.getMagnitudeAction().getValue());

                break;
            }

            case ScalarSourceModel::DefaultRow::DatasetStart:
                break;

            default:
                break;
        }
    });

    // Set preferred point opacity magnitude and offset depending on the scalar source type
    connect(&_opacityAction, &ScalarAction::sourceSelectionChanged, this, [this](const std::uint32_t& sourceSelectionIndex) {
        switch (sourceSelectionIndex)
        {
            case ScalarSourceModel::DefaultRow::Constant:
            {
                _opacityAction.getMagnitudeAction().setValue(50.0f);
                _opacityAction.getSourceAction().getOffsetAction().setValue(0.0f);
                _focusSelection.setChecked(false);

                // Cache the last opacity source index
                //_lastOpacitySourceIndex = sourceSelectionIndex;

                break;
            }

            case ScalarSourceModel::DefaultRow::Selection:
            {
                _opacityAction.getMagnitudeAction().setValue(10.0f);
                _opacityAction.getSourceAction().getOffsetAction().setValue(90.0f);
                _focusSelection.setChecked(true);

                break;
            }

            case ScalarSourceModel::DefaultRow::DatasetStart:
            {
                _focusSelection.setChecked(false);

                // Cache the last opacity source index
                //_lastOpacitySourceIndex = sourceSelectionIndex;

                break;
            }

            default:
                break;
        }
    });

    // Set the current source index to selection when the focus selection is toggled on
    connect(&_focusSelection, &ToggleAction::toggled, this, [this](const bool& toggled) {

        // Get index of the opacity source
        const auto opacitySourceIndex  = _opacityAction.getSourceAction().getPickerAction().getCurrentIndex();

        if (toggled) {

            if (opacitySourceIndex != ScalarSourceModel::DefaultRow::Selection) {
                _opacityAction.setCurrentSourceIndex(ScalarSourceModel::DefaultRow::Selection);
            }

            // Cache the last opacity source index
            _lastOpacitySourceIndex = opacitySourceIndex;
        }
        else {
            if (_lastOpacitySourceIndex != ScalarSourceModel::DefaultRow::Selection)
                _opacityAction.setCurrentSourceIndex(_lastOpacitySourceIndex);
        }
    });
}

QMenu* PointPlotAction::getContextMenu()
{
    auto menu = new QMenu("Plot settings");

    const auto renderMode = getScatterplotWidget().getRenderMode();

    const auto addActionToMenu = [menu](QAction* action) {
        auto actionMenu = new QMenu(action->text());

        actionMenu->addAction(action);

        menu->addMenu(actionMenu);
    };

    addActionToMenu(&_sizeAction);
    addActionToMenu(&_opacityAction);

    return menu;
}

void PointPlotAction::addPointSizeDataset(const Dataset<DatasetImpl>& pointSizeDataset)
{
    // Only proceed if we have a valid point size dataset
    if (!pointSizeDataset.isValid())
        return;

    // Add the dataset to the size action
    _sizeAction.addDataset(pointSizeDataset);
}

void PointPlotAction::addPointOpacityDataset(const Dataset<DatasetImpl>& pointOpacityDataset)
{
    // Only proceed if we have a valid point opacity dataset
    if (!pointOpacityDataset.isValid())
        return;

    // Add the dataset to the opacity action
    _opacityAction.addDataset(pointOpacityDataset);
}

void PointPlotAction::updateDefaultDatasets()
{
    // Get smart pointer to the position dataset
    auto positionDataset = Dataset<Points>(_scatterplotPlugin->getPositionDataset());

    // Only proceed if the position dataset is loaded
    if (!positionDataset.isValid())
        return;

    // Get child data hierarchy items of the position dataset
    const auto children = positionDataset->getDataHierarchyItem().getChildren();

    // Loop over all children and possibly add them to the datasets vector
    for (auto child : children) {

        // Get smart pointer to child dataset
        const auto childDataset = child->getDataset();

        // Get the data type
        const auto dataType = childDataset->getDataType();

        // Add if points/clusters and not derived
        if (dataType != PointType)
            continue;

        // Convert child dataset to points smart pointer
        auto points = Dataset<Points>(childDataset);

        // Add datasets
        _sizeAction.addDataset(points);
        _opacityAction.addDataset(points);
    }
}

void PointPlotAction::updateScatterPlotWidgetPointSizeScalars()
{
    // Number of points
    const auto numberOfPoints = _scatterplotPlugin->getPositionDataset()->getNumPoints();

    // Resize to number of points if needed
    if (numberOfPoints != _pointSizeScalars.size())
        _pointSizeScalars.resize(numberOfPoints);

    // Fill with ones for constant point size
    std::fill(_pointSizeScalars.begin(), _pointSizeScalars.end(), _sizeAction.getMagnitudeAction().getValue());

    // Modulate point size by selection
    if (_sizeAction.isSourceSelection()) {

        // Get smart pointer to current position dataset
        auto positionDataset = _scatterplotPlugin->getPositionDataset();

        // Default point size for all
        std::fill(_pointSizeScalars.begin(), _pointSizeScalars.end(), _sizeAction.getMagnitudeAction().getValue());

        // Establish point size of selected points
        const auto pointSizeSelectedPoints = _sizeAction.getMagnitudeAction().getValue() + _sizeAction.getSourceAction().getOffsetAction().getValue();

        std::vector<bool> selected;

        // Get selected local indices from position dataset
        positionDataset->selectedLocalIndices(positionDataset->getSelection<Points>()->indices, selected);

        for (int i = 0; i < selected.size(); i++) {
            if (!selected[i])
                continue;

            // Selected point size for selected points
            _pointSizeScalars[i] = pointSizeSelectedPoints;
        }
    }

    // Modulate point size by dataset
    if (_sizeAction.isSourceDataset()) {

        // Get current point size source dataset
        auto pointSizeSourceDataset = Dataset<Points>(_sizeAction.getCurrentDataset());

        // Only populate scalars from dataset if the number of points in the source and target dataset match and we have a valid input dataset
        if (pointSizeSourceDataset.isValid() && pointSizeSourceDataset->getNumPoints() == _scatterplotPlugin->getPositionDataset()->getNumPoints())
        {
            // Visit the points dataset to get access to the point values
            pointSizeSourceDataset->visitData([this, pointSizeSourceDataset, numberOfPoints](auto pointData) {

                // Get current dimension index
                const auto currentDimensionIndex = _sizeAction.getSourceAction().getDimensionPickerAction().getCurrentDimensionIndex();

                // Get range for selected dimension
                const auto rangeMin     = _sizeAction.getSourceAction().getRangeAction().getMinimum();
                const auto rangeMax     = _sizeAction.getSourceAction().getRangeAction().getMaximum();
                const auto rangeLength  = rangeMax - rangeMin;

                // Prevent zero division in normalization
                if (rangeLength > 0) {

                    // Loop over all points and compute the point size scalar
                    for (std::uint32_t pointIndex = 0; pointIndex < numberOfPoints; pointIndex++) {

                        // Get point value for dimension
                        auto pointValue = static_cast<float>(pointData[pointIndex][currentDimensionIndex]);

                        // Clamp the point value to the range
                        const auto pointValueClamped = std::max(rangeMin, std::min(rangeMax, pointValue));

                        // Compute normalized point value
                        const auto pointValueNormalized = (pointValueClamped - rangeMin) / rangeLength;

                        // Compute normalized point size scalar
                        _pointSizeScalars[pointIndex] = _sizeAction.getSourceAction().getOffsetAction().getValue() + (pointValueNormalized * _sizeAction.getMagnitudeAction().getValue());
                    }
                }
                else {

                    // Zero division since rangeMin == rangeMax, so reset all point size scalars to constant value
                    std::fill(_pointSizeScalars.begin(), _pointSizeScalars.end(), _sizeAction.getSourceAction().getOffsetAction().getValue() + (rangeMin * _sizeAction.getMagnitudeAction().getValue()));
                }
            });
        }
    }

    // Set scatter plot point size scalars
    _scatterplotPlugin->getScatterplotWidget().setPointSizeScalars(_pointSizeScalars);
}

void PointPlotAction::updateScatterPlotWidgetPointOpacityScalars()
{
    // Enable the opacity magnitude action in constant mode
    //_opacityAction.getMagnitudeAction().setEnabled(_opacityAction.isConstant());

    // Number of points
    const auto numberOfPoints = _scatterplotPlugin->getPositionDataset()->getNumPoints();

    // Resize to number of points
    if (numberOfPoints != _pointOpacityScalars.size())
        _pointOpacityScalars.resize(numberOfPoints);

    // Establish opacity magnitude
    const auto opacityMagnitude = 0.01f * _opacityAction.getMagnitudeAction().getValue();

    // Fill with ones for constant point opacity
    std::fill(_pointOpacityScalars.begin(), _pointOpacityScalars.end(), opacityMagnitude);

    // Modulate point opacity by point selection
    if (_opacityAction.isSourceSelection()) {

        // Get smart pointer to current position dataset
        auto positionDataset = _scatterplotPlugin->getPositionDataset();

        // Get smart pointer to the position selection dataset
        auto selectionSet = positionDataset->getSelection<Points>();

        // Default point opacity for all
        std::fill(_pointOpacityScalars.begin(), _pointOpacityScalars.end(), 0.01f * _opacityAction.getMagnitudeAction().getValue());

        // Establish opacity offset
        const auto opacityOffset = 0.01f * _opacityAction.getSourceAction().getOffsetAction().getValue();

        // Establish point opacity of selected points
        const auto pointOpacitySelectedPoints = std::min(1.0f, opacityMagnitude + opacityOffset);

        // And selected point opacity for selected points
        for (const auto& selectionIndex : selectionSet->indices)
            _pointOpacityScalars[selectionIndex] = pointOpacitySelectedPoints;
    }

    // Modulate point opacity by dataset
    if (_opacityAction.isSourceDataset()) {

        // Get current point opacity source dataset
        auto pointOpacitySourceDataset = Dataset<Points>(_opacityAction.getCurrentDataset());

        // Only populate scalars from dataset if the number of points in the source and target dataset match and we have a valid input dataset
        if (pointOpacitySourceDataset.isValid() && pointOpacitySourceDataset->getNumPoints() == _scatterplotPlugin->getPositionDataset()->getNumPoints())
        {
            // Visit the points dataset to get access to the point values
            pointOpacitySourceDataset->visitData([this, pointOpacitySourceDataset, numberOfPoints, opacityMagnitude](auto pointData) {

                // Get current dimension index
                const auto currentDimensionIndex = _opacityAction.getSourceAction().getDimensionPickerAction().getCurrentDimensionIndex();

                // Get opacity offset
                const auto opacityOffset = 0.01f * _opacityAction.getSourceAction().getOffsetAction().getValue();

                // Get range for selected dimension
                const auto rangeMin     = _opacityAction.getSourceAction().getRangeAction().getMinimum();
                const auto rangeMax     = _opacityAction.getSourceAction().getRangeAction().getMaximum();
                const auto rangeLength  = rangeMax - rangeMin;

                // Prevent zero division in normalization
                if (rangeLength > 0) {

                    // Loop over all points and compute the point size opacity
                    for (std::uint32_t pointIndex = 0; pointIndex < numberOfPoints; pointIndex++) {

                        // Get point value for dimension
                        auto pointValue = static_cast<float>(pointData[pointIndex][currentDimensionIndex]);

                        // Clamp the point value to the range
                        const auto pointValueClamped = std::max(rangeMin, std::min(rangeMax, pointValue));

                        // Compute normalized point value
                        const auto pointValueNormalized = (pointValueClamped - rangeMin) / rangeLength;

                        // Compute normalized point opacity scalar
                        if (opacityOffset == 1.0f)
                            _pointOpacityScalars[pointIndex] = 1.0f;
                        else
                            _pointOpacityScalars[pointIndex] = opacityMagnitude * (opacityOffset + (pointValueNormalized / (1.0f - opacityOffset)));
                    }
                }
                else {

                    // Get reference to range action
                    auto& rangeAction = _opacityAction.getSourceAction().getRangeAction();

                    // Handle zero division
                    if (rangeAction.getRangeMinAction().getValue() == rangeAction.getRangeMaxAction().getValue())
                        std::fill(_pointOpacityScalars.begin(), _pointOpacityScalars.end(), 0.0f);
                    else
                        std::fill(_pointOpacityScalars.begin(), _pointOpacityScalars.end(), 1.0f);
                }
            });
        }
    }

    // Set scatter plot point size scalars
    _scatterplotPlugin->getScatterplotWidget().setPointOpacityScalars(_pointOpacityScalars);
}

PointPlotAction::Widget::Widget(QWidget* parent, PointPlotAction* pointPlotAction, const std::int32_t& widgetFlags) :
    WidgetActionWidget(parent, pointPlotAction, widgetFlags)
{
    setToolTip("Point plot settings");

    // Add widgets
    if (widgetFlags & PopupLayout) {
        auto layout = new QGridLayout();

        layout->setContentsMargins(0, 0, 0, 0);

        layout->addWidget(pointPlotAction->getSizeAction().createLabelWidget(this), 0, 0);
        layout->addWidget(pointPlotAction->getSizeAction().createWidget(this), 0, 1);

        layout->addWidget(pointPlotAction->getOpacityAction().createLabelWidget(this), 1, 0);
        layout->addWidget(pointPlotAction->getOpacityAction().createWidget(this), 1, 1);

        setLayout(layout);
    }
    else {
        auto layout = new QHBoxLayout();

        layout->setContentsMargins(0, 0, 0, 0);

        layout->addWidget(pointPlotAction->getSizeAction().createLabelWidget(this));
        layout->addWidget(pointPlotAction->getSizeAction().createWidget(this));
        layout->addWidget(pointPlotAction->getOpacityAction().createLabelWidget(this));
        layout->addWidget(pointPlotAction->getOpacityAction().createWidget(this));

        setLayout(layout);
    }
}
