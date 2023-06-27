#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"
#include "ProjectionView.h"
#include "DataHierarchyItem.h"
#include "Application.h"
#include "actions/GroupsAction.h"

#include "util/PixelSelectionTool.h"
#include "PointData/DimensionsPickerAction.h"

#include "PointData/PointData.h"
#include "ColorData/ColorData.h"
#include "ClusterData/ClusterData.h"

#include "graphics/Vector2f.h"
#include "graphics/Vector3f.h"
#include "widgets/DropWidget.h"

#include <Eigen/Dense>
#include "Compute/LocalDimensionality.h"
#include "Compute/RandomWalks.h"
#include "Compute/DataTransformations.h"
#include "Timer.h"
#include "Types.h"

#include <QtCore>
#include <QApplication>
#include <QDebug>
#include <QMenu>
#include <QAction>
#include <QMetaType>
#include <QVector>
#include <QFileDialog>

#include <algorithm>
#include <functional>
#include <limits>
#include <set>
#include <vector>
#include <iostream>
#include <random>
#include <chrono>

Q_PLUGIN_METADATA(IID "nl.biovault.GradientExplorerPlugin")

using namespace hdps;
using namespace hdps::util;

namespace
{
    void computeDirection(DataMatrix& dataMatrix, DataMatrix& projMatrix, KnnGraph& knnGraph, int numSteps, std::vector<Vector2f>& directions)
    {
        for (int p = 0; p < dataMatrix.rows(); p++)
        {
            FloodFill directionsFloodFill(numSteps);
            directionsFloodFill.compute(knnGraph, p);

            int depth = 3;
            bigint numNodes = 0;
            for (int i = 0; i < depth; i++)
            {
                numNodes += directionsFloodFill.getWaves()[i].size();
            }

            Eigen::MatrixXf nodeLocations(numNodes, 2);
            int n = 0;
            for (int i = 0; i < depth; i++)
            {
                for (int j = 0; j < directionsFloodFill.getWaves()[i].size(); j++)
                {
                    int index = directionsFloodFill.getWaves()[i][j];

                    nodeLocations(n, 0) = projMatrix(index, 0);
                    nodeLocations(n, 1) = projMatrix(index, 1);
                    n++;
                }
            }

            //
            // Mean centering data.
            Eigen::MatrixXf centered = nodeLocations.rowwise() - nodeLocations.row(0);
            // Compute the covariance matrix.
            Eigen::MatrixXf cov = centered.adjoint() * centered;
            cov = cov / (nodeLocations.rows() - 1);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(cov);
            // Normalize eigenvalues to make them represent percentages.
            Eigen::MatrixXf evecs = eig.eigenvectors();
            // Get the two major eigenvectors and omit the others.
            Vector2f majorEigenVector;
            if (eig.eigenvalues()(0) > eig.eigenvalues()(1))
                majorEigenVector.set(evecs(0, 0), evecs(1, 0));
            else
                majorEigenVector.set(evecs(0, 1), evecs(1, 1));

            directions.push_back(Vector2f(projMatrix(p, 0), projMatrix(p, 1)));
            directions.push_back(majorEigenVector);
        }
    }
}

ScatterplotPlugin::ScatterplotPlugin(const PluginFactory* factory) :
    ViewPlugin(factory),
    _positionDataset(),
    _positionSourceDataset(),
    _positions(),
    _numPoints(0),
    _scatterPlotWidget(new ScatterplotWidget()),
    _projectionViews(2, nullptr),
    _selectedView(),
    _dropWidget(nullptr),
    _settingsAction(this),
    _gradientGraph(new GradientGraph()),
    _selectedDimension(-1),
    _floodFill(10),
    _filterType(filters::FilterType::SPATIAL_PEAK),
    _overlayType(OverlayType::NONE),
    _colorMapAction(this, "Color map", util::ColorMap::Type::OneDimensional, "RdYlBu", "RdYlBu"),
    _graphTimer(new QTimer(this)),
    _filterLabel(nullptr)
{
    setObjectName("GradientExplorer");

    _dropWidget = new DropWidget(_scatterPlotWidget);

    for (int i = 0; i < _projectionViews.size(); i++)
    {
        _projectionViews[i] = new ProjectionView();
    }
    _selectedView = new ProjectionView();

    // Connect signals from views
    connect(_projectionViews[0], &ProjectionView::viewSelected, this, [this]() { _selectedViewIndex = 1; updateViews(); });
    connect(_projectionViews[1], &ProjectionView::viewSelected, this, [this]() { _selectedViewIndex = 2; updateViews(); });
    connect(_selectedView, &ProjectionView::viewSelected, this, [this]() { _selectedViewIndex = 3; updateViews(); });

    connect(_scatterPlotWidget, &ScatterplotWidget::initialized, this, [this]()
    {
        getScatterplotWidget().setColorMap(_colorMapAction.getColorMapImage().mirrored(false, true));
        getScatterplotWidget().setScalarEffect(PointEffect::Color);
    });
    connect(_projectionViews[0], &ProjectionView::initialized, this, [this]() {_projectionViews[0]->setColorMap(_colorMapAction.getColorMapImage().mirrored(false, true)); });
    connect(_projectionViews[1], &ProjectionView::initialized, this, [this]() {_projectionViews[1]->setColorMap(_colorMapAction.getColorMapImage().mirrored(false, true)); });
    connect(_selectedView, &ProjectionView::initialized, this, [this]() {_selectedView->setColorMap(_colorMapAction.getColorMapImage().mirrored(false, true)); });

    connect(_gradientGraph, &GradientGraph::lineClicked, this, &ScatterplotPlugin::onLineClicked);
    _graphTimer->setSingleShot(true);
    connect(_graphTimer, &QTimer::timeout, this, &ScatterplotPlugin::computeGraphs);

    getWidget().setFocusPolicy(Qt::ClickFocus);

    connect(_scatterPlotWidget, &ScatterplotWidget::customContextMenuRequested, this, [this](const QPoint& point) {
        if (!_positionDataset.isValid())
            return;

        auto contextMenu = _settingsAction.getContextMenu();
        
        contextMenu->addSeparator();

        _positionDataset->populateContextMenu(contextMenu);

        contextMenu->exec(getWidget().mapToGlobal(point));
    });

    _dropWidget->setDropIndicatorWidget(new DropWidget::DropIndicatorWidget(&getWidget(), "No data loaded", "Drag an item from the data hierarchy and drop it here to visualize data..."));
    _dropWidget->initialize([this](const QMimeData* mimeData) -> DropWidget::DropRegions {
        DropWidget::DropRegions dropRegions;

        const auto mimeText = mimeData->text();
        const auto tokens   = mimeText.split("\n");

        if (tokens.count() == 1)
            return dropRegions;

        const auto datasetGuiName       = tokens[0];
        const auto datasetId            = tokens[1];
        const auto dataType             = DataType(tokens[2]);
        const auto dataTypes            = DataTypes({ PointType , ColorType, ClusterType });

        // Check if the data type can be dropped
        if (!dataTypes.contains(dataType))
            dropRegions << new DropWidget::DropRegion(this, "Incompatible data", "This type of data is not supported", "exclamation-circle", false);

        // Points dataset is about to be dropped
        if (dataType == PointType) {

            // Get points dataset from the core
            auto candidateDataset = _core->requestDataset<Points>(datasetId);

            // Establish drop region description
            const auto description = QString("Visualize %1 as points or density/contour map").arg(datasetGuiName);

            if (!_positionDataset.isValid()) {

                // Load as point positions when no dataset is currently loaded
                dropRegions << new DropWidget::DropRegion(this, "Point position", description, "map-marker-alt", true, [this, candidateDataset]() {
                    _positionDataset = candidateDataset;
                    positionDatasetChanged();
                });
            }
            else {
                if (_positionDataset != candidateDataset && candidateDataset->getNumDimensions() >= 2) {

                    // The number of points is equal, so offer the option to replace the existing points dataset
                    dropRegions << new DropWidget::DropRegion(this, "Point position", description, "map-marker-alt", true, [this, candidateDataset]() {
                        _positionDataset = candidateDataset;
                        positionDatasetChanged();
                    });
                }

                if (candidateDataset->getNumPoints() == _positionDataset->getNumPoints()) {

                    //// The number of points is equal, so offer the option to use the points dataset as source for points colors
                    //dropRegions << new DropWidget::DropRegion(this, "Point color", QString("Colorize %1 points with %2").arg(_positionDataset->getGuiName(), candidateDataset->getGuiName()), "palette", true, [this, candidateDataset]() {
                    //    _settingsAction.getColoringAction().addColorDataset(candidateDataset);
                    //    _settingsAction.getColoringAction().setCurrentColorDataset(candidateDataset);
                    //});

                    // The number of points is equal, so offer the option to use the points dataset as source for points size
                    dropRegions << new DropWidget::DropRegion(this, "Point size", QString("Size %1 points with %2").arg(_positionDataset->getGuiName(), candidateDataset->getGuiName()), "ruler-horizontal", true, [this, candidateDataset]() {
                        _settingsAction.getPlotAction().getPointPlotAction().addPointSizeDataset(candidateDataset);
                        _settingsAction.getPlotAction().getPointPlotAction().getSizeAction().setCurrentDataset(candidateDataset);
                    });

                    // The number of points is equal, so offer the option to use the points dataset as source for points opacity
                    dropRegions << new DropWidget::DropRegion(this, "Point opacity", QString("Set %1 points opacity with %2").arg(_positionDataset->getGuiName(), candidateDataset->getGuiName()), "brush", true, [this, candidateDataset]() {
                        _settingsAction.getPlotAction().getPointPlotAction().addPointOpacityDataset(candidateDataset);
                        _settingsAction.getPlotAction().getPointPlotAction().getOpacityAction().setCurrentDataset(candidateDataset);
                    });
                }
            }
        }

        return dropRegions;
    });

    _scatterPlotWidget->installEventFilter(this);
}

ScatterplotPlugin::~ScatterplotPlugin()
{
}

void ScatterplotPlugin::init()
{
    auto layout = new QVBoxLayout();
    auto gradientViewLayout = new QVBoxLayout();
    auto dimensionViewsLayout = new QHBoxLayout();

    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    layout->addWidget(_settingsAction.createWidget(&getWidget()));

    _filterLabel = new QLabel();
    _filterLabel->setText("Spatial Peak Ranking");
    QFont font = _filterLabel->font();
    font.setPointSize(font.pointSize() * 2);
    _filterLabel->setFont(font);

    gradientViewLayout->setContentsMargins(6, 0, 6, 0);
    gradientViewLayout->addWidget(_filterLabel);
    gradientViewLayout->addWidget(_projectionViews[0], 50);
    gradientViewLayout->addWidget(_projectionViews[1], 50);

    // Dimension selection
    QLabel* dimensionSelectionLabel = new QLabel("Dimension Selection");
    dimensionSelectionLabel->setFont(font);
    gradientViewLayout->addWidget(dimensionSelectionLabel);
    gradientViewLayout->addWidget(_selectedView, 50);

    // Expression graph
    QLabel* sortedExpressionGraphLabel = new QLabel("Sorted Expression Graph");
    sortedExpressionGraphLabel->setFont(font);
    gradientViewLayout->addWidget(sortedExpressionGraphLabel);
    gradientViewLayout->addWidget(_gradientGraph, 70);

    auto leftPanel = new QVBoxLayout();
    leftPanel->addWidget(_scatterPlotWidget, 90);

    auto centralPanel = new QHBoxLayout();
    centralPanel->addLayout(leftPanel, 80);
    centralPanel->addLayout(gradientViewLayout, 20);

    layout->addLayout(centralPanel, 100);

    auto bottomToolbarWidget = new QWidget();
    auto bottomToolbarLayout = new QHBoxLayout();

    bottomToolbarWidget->setAutoFillBackground(true);
    bottomToolbarWidget->setLayout(bottomToolbarLayout);

    bottomToolbarLayout->setContentsMargins(0, 0, 0, 0);
    //bottomToolbarLayout->addWidget(_settingsAction.getColoringAction().getColorMapAction().createLabelWidget(&getWidget()));
    //bottomToolbarLayout->addWidget(_settingsAction.getColoringAction().getColorMapAction().createWidget(&getWidget()));
    bottomToolbarLayout->addWidget(_settingsAction.getPlotAction().getPointPlotAction().getFocusSelection().createWidget(&getWidget()));
    bottomToolbarLayout->addStretch(1);
    bottomToolbarLayout->addWidget(_settingsAction.getExportAction().createWidget(&getWidget()));
    bottomToolbarLayout->addWidget(_settingsAction.getMiscellaneousAction().createCollapsedWidget(&getWidget()));

    layout->addWidget(bottomToolbarWidget, 1);

    getWidget().setLayout(layout);

    // Update the data when the scatter plot widget is initialized
    connect(_scatterPlotWidget, &ScatterplotWidget::initialized, this, &ScatterplotPlugin::updateData);

    //_eventListener.setEventCore(Application::core());
    _eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DataSelectionChanged));
    _eventListener.registerDataEventByType(PointType, std::bind(&ScatterplotPlugin::onDataEvent, this, std::placeholders::_1));

    // Load points when the pointer to the position dataset changes
    //connect(&_positionDataset, &Dataset<Points>::changed, this, &ScatterplotPlugin::positionDatasetChanged);

    // Update points when the position dataset data changes
    //connect(&_positionDataset, &Dataset<Points>::dataChanged, this, &ScatterplotPlugin::updateData);

    // Update point selection when the position dataset data changes
    connect(&_positionDataset, &Dataset<Points>::dataSelectionChanged, this, &ScatterplotPlugin::updateSelection);

    // Update the window title when the GUI name of the position dataset changes
    connect(&_positionDataset, &Dataset<Points>::dataGuiNameChanged, this, &ScatterplotPlugin::updateWindowTitle);

    // Do an initial update of the window title
    updateWindowTitle();

    _floodScalars = _core->addDataset<Points>("Points", "Flood Nodes");

    events().notifyDatasetAdded(_floodScalars);
}

void ScatterplotPlugin::updateWindowTitle()
{
    if (!_positionDataset.isValid())
        getWidget().setWindowTitle(getGuiName());
    else
        getWidget().setWindowTitle(QString("%1: %2").arg(getGuiName(), _positionDataset->getDataHierarchyItem().getFullPathName()));
}

void ScatterplotPlugin::updateColorMapActionScalarRange()
{
    // Get the color map range from the scatter plot widget
    const auto colorMapRange = getScatterplotWidget().getColorMapRange();
    const auto colorMapRangeMin = colorMapRange.x;
    const auto colorMapRangeMax = colorMapRange.y;

    // Get reference to color map range action
    const auto& colorMapRangeAction = _colorMapAction.getRangeAction(ColorMapAction::Axis::X);

    // Initialize the color map range action with the color map range from the scatter plot 
    _colorMapAction.getDataRangeAction(ColorMapAction::Axis::X).setRange({ colorMapRangeMin, colorMapRangeMax });
}

void ScatterplotPlugin::loadData(const Datasets& datasets)
{
    // Exit if there is nothing to load
    if (datasets.isEmpty())
        return;

    // Load the first dataset
    _positionDataset = datasets.first();

    // And set the coloring mode to constant
    //_settingsAction.getColoringAction().getColorByAction().setCurrentIndex(1);
}

void ScatterplotPlugin::onDataEvent(hdps::DataEvent* dataEvent)
{
    if (dataEvent->getType() == EventType::DataSelectionChanged)
    {
        if (dataEvent->getDataset() == _positionDataset)
        {
            if (_positionDataset->isDerivedData())
            {
                onPointSelection();
            }
        }
    }
}

void ScatterplotPlugin::positionDatasetChanged()
{
    // Only proceed if we have a valid position dataset
    if (!_positionDataset.isValid())
        return;

    // Reset dataset references
    _positionSourceDataset.reset();

    // Set position source dataset reference when the position dataset is derived
    if (_positionDataset->isDerivedData())
        _positionSourceDataset = _positionDataset->getSourceDataset<Points>();

    // Do not show the drop indicator if there is a valid point positions dataset
    _dropWidget->setShowDropIndicator(!_positionDataset.isValid());

    _dataInitialized = true;

    // Update positions data
    updateData();

    // Update the window title to reflect the position dataset change
    updateWindowTitle();

    computeStaticData();
}

std::uint32_t ScatterplotPlugin::getNumberOfPoints() const
{
    if (!_positionDataset.isValid())
        return 0;

    return _positionDataset->getNumPoints();
}

void ScatterplotPlugin::createSubset(const bool& fromSourceData /*= false*/, const QString& name /*= ""*/)
{
    auto subsetPoints = fromSourceData ? _positionDataset->getSourceDataset<Points>() : _positionDataset;

    // Create the subset
    auto subset = subsetPoints->createSubsetFromSelection(_positionDataset->getGuiName(), _positionDataset);

    // Notify others that the subset was added
    events().notifyDatasetAdded(subset);

    // And select the subset
    subset->getDataHierarchyItem().select();
}

void ScatterplotPlugin::computeStaticData()
{
    Timer timer;
    timer.start();

    std::cout << "Start conversion" << std::endl;
    convertToEigenMatrix(_positionDataset, _positionSourceDataset, _dataMatrix);
    convertToEigenMatrixProjection(_positionDataset, _fullProjMatrix);

    int xDim = _settingsAction.getPositionAction().getDimensionX();
    int yDim = _settingsAction.getPositionAction().getDimensionY();
    _projMatrix = _fullProjMatrix(Eigen::all, std::vector<int> { xDim, yDim });

    // Set mask to include all points
    //_mask.resize(_dataMatrix.rows());
    //std::iota(_mask.begin(), _mask.end(), 0);

    std::cout << "Number of enabled dimensions in the dataset : " << _dataMatrix.cols() << std::endl;

    timer.mark("Data preparation");
    qDebug() << "Standardize..";
    // Standardize
    standardizeData(_dataMatrix, _variances);
    qDebug() << "Normalize..";
    // Compute normalized data
    normalizeData(_dataMatrix, _normalizedData);

    _bins.resize(_dataMatrix.cols(), std::vector<int>(30));

    Bounds bounds = _scatterPlotWidget->getBounds();
    _projectionSize = bounds.getWidth() > bounds.getHeight() ? bounds.getWidth() : bounds.getHeight();
    std::cout << "Projection size: " << _projectionSize << std::endl;

    timer.mark("Data transformations");

    if (_computeOnLoad)
    {
        createKnnIndex();
        timer.mark("Computing KNN index");
        computeKnnGraph();
    }

    timer.mark("Computing KNN graph");

    _largeKnnGraph.writeToFile();

    //_localSpatialDimensionality.clear();
    //_localHighDimensionality.clear();

    //timer.mark("Local dimensionality");

    //compute::computeHDLocalDimensionality(_dataMatrix, _largeKnnGraph, _localHighDimensionality);
    //getScatterplotWidget().setColorMap(_colorMapAction.getColorMapImage());
    //getScatterplotWidget().setColorMapRange(0, 1);
    //getScatterplotWidget().setScalars(_localHighDimensionality);

    //Dataset<Points> localDims = _core->addDataset("Points", "Dimensionality");

    //localDims->setData(_localHighDimensionality.data(), _localHighDimensionality.size(), 1);

    //_core->notifyDatasetAdded(localDims);

    //computeDirection(_dataMatrix, _projMatrix, _knnGraph, _directions);
    //getScatterplotWidget().setDirections(_directions);

    timer.mark("Directions");

    qDebug() << _positionSourceDataset->getDimensionNames().size() << _positionSourceDataset->getNumDimensions();

    // Get enabled dimension names
    const auto& dimNames = _positionSourceDataset->getDimensionNames();
    auto enabledDimensions = _positionSourceDataset->getDimensionsPickerAction().getEnabledDimensions();

    _enabledDimNames.clear();
    for (int i = 0; i < enabledDimensions.size(); i++)
    {
        if (enabledDimensions[i])
            _enabledDimNames.push_back(dimNames[i]);
    }

    // Set up chart
    _gradientGraph->setNumDimensions((bigint)_enabledDimNames.size());

    timer.finish("Graph init");
}

void ScatterplotPlugin::onPointSelection()
{
    Timer timer;

    if (!_positionDataset.isValid() || !_positionSourceDataset.isValid() || !_dataInitialized)
        return;

    hdps::Dataset<Points> selection = _positionSourceDataset->getSelection();

    int numDimensions = _dataMatrix.cols();

    if (selection->indices.size() > 0)
    {
timer.start();
        //_selectedPoint = selection->indices[0];
        Vector2f center = _positions[_globalSelectedPoint];

        getScatterplotWidget().setCurrentPosition(center);
        getProjectionViews()[0]->setCurrentPosition(center);
        getProjectionViews()[1]->setCurrentPosition(center);
        _selectedView->setCurrentPosition(center);
        getScatterplotWidget().setFilterRadii(Vector2f(_spatialPeakFilter.getInnerFilterRadius() * _projectionSize, _spatialPeakFilter.getOuterFilterRadius() * _projectionSize));

        KnnGraph& knnGraph = _mask.empty() ? _knnGraph : _maskedKnnGraph;
        DataMatrix& dataMatrix = _mask.empty() ? _dataMatrix : _maskedDataMatrix;
        DataMatrix& projMatrix = _mask.empty() ? _projMatrix : _maskedProjMatrix;

        /////////////////////
        // Gradient picker //
        /////////////////////
        std::vector<int> dimRanking;
        switch (_filterType)
        {
        case filters::FilterType::SPATIAL_PEAK:
        {
            if (_settingsAction.getFilterAction().getRestrictToFloodAction().isChecked())
                _spatialPeakFilter.computeDimensionRanking(_selectedPoint, dataMatrix, _variances, projMatrix, _projectionSize, dimRanking, _floodFill.getAllNodes());
            else
                _spatialPeakFilter.computeDimensionRanking(_selectedPoint, dataMatrix, _variances, projMatrix, _projectionSize, dimRanking);
            break;
        }
        case filters::FilterType::HD_PEAK:
        {
            _hdFloodPeakFilter.computeDimensionRanking(_selectedPoint, dataMatrix, _variances, _floodFill, dimRanking);
            break;
        }
        }
timer.mark("Ranking");

        // Set appropriate coloring of gradient view, FIXME use colormap later
        for (int pi = 0; pi < _projectionViews.size(); pi++)
        {
            const auto dimValues = _dataMatrix(Eigen::all, dimRanking[pi]);
            std::vector<float> dimV(dimValues.data(), dimValues.data() + dimValues.size());
            _projectionViews[pi]->setShownDimension(dimRanking[pi]);
            _projectionViews[pi]->setScalars(dimV, _globalSelectedPoint);
            _projectionViews[pi]->setProjectionName(_enabledDimNames[dimRanking[pi]]);
        }
        // Set selected gradient view
        if (_selectedDimension >= 0)
        {
            qDebug() << "SEL DIM:" << _selectedDimension;
            const auto dimValues = _dataMatrix(Eigen::all, _selectedDimension);
            std::vector<float> dimV(dimValues.data(), dimValues.data() + dimValues.size());
            _selectedView->setShownDimension(_selectedDimension);
            _selectedView->setScalars(dimV, _globalSelectedPoint);
            _selectedView->setProjectionName(_enabledDimNames[_selectedDimension]);
        }

        _gradientGraph->setTopDimensions(dimRanking[0], dimRanking[1]);

timer.mark("Filter");

        //////////////////
        // Do floodfill //
        //////////////////

        if (_graphAvailable)
            _floodFill.compute(knnGraph, _selectedPoint);

        timer.mark("Floodfill");

        /////////////////////
        // Coloring        //
        /////////////////////
        std::vector<float> colorScalars(_dataMatrix.rows(), 0);

        if (_graphAvailable)
        {
            switch (_overlayType)
            {
            case OverlayType::NONE:
            {
                _scatterPlotWidget->setColoredBy("Colored by - Flood fill step");
                for (int i = 0; i < _floodFill.getNumWaves(); i++)
                {
                    for (int j = 0; j < _floodFill.getWaves()[i].size(); j++)
                    {
                        int index = _floodFill.getWaves()[i][j];
                        colorScalars[_mask.empty() ? index : _mask[index]] = 1 - (1.0f / _floodFill.getNumWaves()) * i;
                    }
                }

                _scatterPlotWidget->setColorMapRange(0, 1);
                break;
            }
            case OverlayType::DIM_VALUES:
            {
                _scatterPlotWidget->setColoredBy("Colored by - Dim: " + _enabledDimNames[dimRanking[0]]);
                for (int i = 0; i < _floodFill.getTotalNumNodes(); i++)
                {
                    int node = _floodFill.getAllNodes()[i];
                    int index = _mask.empty() ? node : _mask[node];
                    colorScalars[index] = _normalizedData[dimRanking[0]][index];
                }
                break;
            }
            case OverlayType::LOCAL_DIMENSIONALITY:
            {
                _scatterPlotWidget->setColoredBy("Colored by - Local Dimensionality");
                if (_localHighDimensionality.empty()) break;

                for (int i = 0; i < _floodFill.getTotalNumNodes(); i++)
                {
                    int node = _floodFill.getAllNodes()[i];
                    int index = _mask.empty() ? node : _mask[node];
                    colorScalars[index] = _localHighDimensionality[node];
                }
                break;
            }
            case OverlayType::DIRECTIONS:
            {
                std::vector<float> scalars(dataMatrix.rows(), 0);
                std::vector<Vector2f> directions(_floodFill.getTotalNumNodes() * 2);

                for (int i = 0; i < _floodFill.getTotalNumNodes(); i++)
                {
                    float idx = _floodFill.getAllNodes()[i];
                    directions[i * 2 + 0] = _directions[idx * 2 + 0];
                    directions[i * 2 + 1] = _directions[idx * 2 + 1];
                    scalars[idx] = 0.5f;
                }
                getScatterplotWidget().setDirections(directions);

                break;
            }
            }
        }
        
        getScatterplotWidget().setScalars(colorScalars);

        // Store scalars in floodfill dataset
        _floodScalars->setData<float>(colorScalars.data(), colorScalars.size(), 1);
        events().notifyDatasetChanged(_floodScalars);

        timer.mark("Overlay");

        /////////////////////
        // Graphs          //
        /////////////////////

        // Start a timer to compute the graphs in 100ms, if the timer is restarted before graphs are not computed
        _graphTimer->start(100);

        timer.finish("Graphs");
    }
}

void ScatterplotPlugin::updateData()
{
    // Check if the scatter plot is initialized, if not, don't do anything
    if (!_scatterPlotWidget->isInitialized() || !_dataInitialized)
        return;
    
    // If no dataset has been selected, don't do anything
    if (_positionDataset.isValid()) {

        // Get the selected dimensions to use as X and Y dimension in the plot
        const dint xDim = _settingsAction.getPositionAction().getDimensionX();
        const dint yDim = _settingsAction.getPositionAction().getDimensionY();

        // If one of the dimensions was not set, do not draw anything
        if (xDim < 0 || yDim < 0)
            return;

        // Determine number of points depending on if its a full dataset or a subset
        _numPoints = _positionDataset->getNumPoints();

        // Extract 2-dimensional points from the data set based on the selected dimensions
        calculatePositions(*_positionDataset);

        // Pass the 2D points to the scatter plot widget
        _scatterPlotWidget->setData(&_positions);
        for (int i = 0; i < (int) _projectionViews.size(); i++)
            _projectionViews[i]->setData(&_positions);
        _selectedView->setData(&_positions);

        updateSelection();

        // Update proj matrix
        _projMatrix = _fullProjMatrix(Eigen::all, std::vector<dint> { xDim, yDim });

        Bounds bounds = _scatterPlotWidget->getBounds();
        _projectionSize = bounds.getWidth() > bounds.getHeight() ? bounds.getWidth() : bounds.getHeight();
        std::cout << "Projection size: " << _projectionSize << std::endl;

        //std::vector<Vector2f> directions;
        //computeDirection(_dataMatrix, _projMatrix, _knnGraph, directions);
        //getScatterplotWidget().setDirections(directions);
    }
    else {
        _positions.clear();
        _scatterPlotWidget->setData(&_positions);
    }
}

void ScatterplotPlugin::calculatePositions(const Points& points)
{
    points.extractDataForDimensions(_positions, _settingsAction.getPositionAction().getDimensionX(), _settingsAction.getPositionAction().getDimensionY());
}

void ScatterplotPlugin::updateSelection()
{
    if (!_positionDataset.isValid())
        return;

    auto selection = _positionDataset->getSelection<Points>();

    std::vector<bool> selected;
    std::vector<char> highlights;

    _positionDataset->selectedLocalIndices(selection->indices, selected);

    highlights.resize(_positionDataset->getNumPoints(), 0);

    for (int i = 0; i < selected.size(); i++)
        highlights[i] = selected[i] ? 1 : 0;

    _scatterPlotWidget->setHighlights(highlights, static_cast<std::int32_t>(selection->indices.size()));
}

void ScatterplotPlugin::updateViews()
{
    _projectionViews[0]->selectView(false);
    _projectionViews[1]->selectView(false);
    _selectedView->selectView(false);

    ProjectionView* selectedView = nullptr;
    if (_selectedViewIndex == 1) selectedView = _projectionViews[0];
    else if (_selectedViewIndex == 2) selectedView = _projectionViews[1];
    else if (_selectedViewIndex == 3) selectedView = _selectedView;

    if (selectedView != nullptr)
    {
        selectedView->selectView(true);

        int selectedDimension = selectedView->getShownDimension();
        const auto dimValues = _dataMatrix(Eigen::all, selectedDimension);
        std::vector<float> dimV(dimValues.data(), dimValues.data() + dimValues.size());
        getScatterplotWidget().setScalars(dimV);
        getScatterplotWidget().setProjectionName("Dimension View: " + _enabledDimNames[selectedDimension]);
        getScatterplotWidget().setColoredBy("");
        //setProjectionName(_enabledDimNames[_selectedDimension]);
    }
    else
    {
        getScatterplotWidget().setProjectionName("Floodfill View");
    }
}

bool ScatterplotPlugin::eventFilter(QObject* target, QEvent* event)
{
    auto shouldPaint = false;

    switch (event->type())
    {
    case QEvent::Resize:
    {
        const auto resizeEvent = static_cast<QResizeEvent*>(event);

        break;
    }

    case QEvent::MouseButtonPress:
    {
        auto mouseEvent = static_cast<QMouseEvent*>(event);

        qDebug() << "Mouse button press";

        _mousePressed = true;

        _selectedViewIndex = 0;
        updateViews();

        break;
    }

    case QEvent::MouseButtonRelease:
    {
        auto mouseEvent = static_cast<QMouseEvent*>(event);

        _mousePressed = false;

        break;
    }

    case QEvent::MouseMove:
    {
        auto mouseEvent = static_cast<QMouseEvent*>(event);

        if (!_mousePressed)
            break;

        _mousePos = QPoint(mouseEvent->position().x(), mouseEvent->position().y());

        hdps::Bounds bounds = _scatterPlotWidget->getBounds();

        if (!_positionDataset.isValid())
            return QObject::eventFilter(target, event);

        // Define function finding out which points are selected or not
        auto findSelectedPoint = [this, bounds](std::vector<int> mask) -> int
        {
            const auto w = _scatterPlotWidget->width();
            const auto h = _scatterPlotWidget->height();
            const auto size = w < h ? w : h;

            int closestIndex = 0;
            float minDist = std::numeric_limits<float>::max();

            for (int i = 0; i < mask.size(); i++)
            //for (const int maskIndex : mask)
            {
                int maskIndex = mask[i];
                const Vector2f& position = _positions[maskIndex];
                const auto pointUV = Vector2f((position.x - bounds.getLeft()) / bounds.getWidth(), (bounds.getTop() - position.y) / bounds.getHeight());
                const auto uvOffset = Vector2f((w - size) / 2.0f, (h - size) / 2.0f);
                const auto uv = uvOffset + Vector2f(pointUV.x * size, pointUV.y * size);

                Vector2f diff = uv - Vector2f(_mousePos.x(), _mousePos.y());
                float sqrDist = diff.x * diff.x + diff.y * diff.y;
                if (sqrDist < minDist)
                {
                    closestIndex = i;
                    minDist = sqrDist;
                }
            }
            return closestIndex;
        };

        // Loop over either all points or only the masked points and establish whether they are selected or not
        std::vector<int> full(_dataMatrix.rows());
        std::iota(full.begin(), full.end(), 0);
        _selectedPoint = findSelectedPoint(_mask.empty() ? full: _mask);
        _globalSelectedPoint = _mask.empty() ? _selectedPoint : _mask[_selectedPoint];
        int selectedPoint = _globalSelectedPoint;

        // If we're looking at a subset of the projection, apply subset indirection to selected index
        if (!_positionDataset->isFull())
            selectedPoint = _positionDataset->indices[selectedPoint];

        // Create vector for target selection indices
        std::vector<std::uint32_t> targetSelectionIndices;
        targetSelectionIndices.push_back(selectedPoint);

        // Apply the selection indices
        _positionDataset->setSelectionIndices(targetSelectionIndices);

        // Notify others that the selection changed
        events().notifyDatasetSelectionChanged(_positionDataset);

        //_explanationWidget->update();

        break;
    }

    default:
        break;
    }

    return QObject::eventFilter(target, event);
}

/******************************************************************************
 * Graphs
 ******************************************************************************/

void ScatterplotPlugin::computeGraphs()
{
    // Binning
    int numBins = (int)_bins.size();
    int binSteps = (int)_bins[0].size();

    for (int d = 0; d < numBins; d++)
        std::fill(_bins[d].begin(), _bins[d].end(), 0);

#pragma omp parallel for
    for (int d = 0; d < numBins; d++)
    {
        int* const bins_d = &_bins[d][0];
        float* const norm_d = &_normalizedData[d][0];

        for (bigint i = 0; i < _floodFill.getTotalNumNodes(); i++)
        {
            const float& f = norm_d[_floodFill.getAllNodes()[i]];
            bins_d[(bigint)(f * binSteps)]++;
        }
    }
    qDebug() << "Graphs computed";
    _gradientGraph->setBins(_bins);
}

void ScatterplotPlugin::onLineClicked(dint dim)
{
    qDebug() << "Dim: " << dim;
    _selectedDimension = dim;

    //const auto dimValues = _dataMatrix(Eigen::all, _selectedDimension);

    //float minV = *std::min_element(dimValues.begin(), dimValues.end());
    //float maxV = *std::max_element(dimValues.begin(), dimValues.end());

    //std::vector<Vector3f> colors(_positions.size());
    //for (int i = 0; i < dimValues.size(); i++)
    //{
    //    float dimValue = dimValues[i] / (maxV - minV);

    //    colors[i] = (i == _selectedPoint) ? Vector3f(1, 0, 0) : Vector3f(1 - dimValue, 1 - dimValue, 1 - dimValue);
    //}
    //_selectedView->setColors(colors);

    _selectedView->setProjectionName(_enabledDimNames[_selectedDimension]);

    onPointSelection();
    updateViews();
}

/******************************************************************************
 * Import / Export
 ******************************************************************************/

void ScatterplotPlugin::exportRankings()
{
    std::vector<std::vector<int>> perPointDimRankings(_dataMatrix.rows());
    for (int i = 0; i < _dataMatrix.rows(); i++)
    {
        FloodFill exportFloodFill(_floodFill.getNumWaves());
        exportFloodFill.compute(_knnGraph, i);

        switch (_filterType)
        {
        case filters::FilterType::SPATIAL_PEAK:
            if (_settingsAction.getFilterAction().getRestrictToFloodAction().isChecked())
                _spatialPeakFilter.computeDimensionRanking(i, _dataMatrix, _variances, _projMatrix, _projectionSize, perPointDimRankings[i], exportFloodFill.getAllNodes());
            else
                _spatialPeakFilter.computeDimensionRanking(i, _dataMatrix, _variances, _projMatrix, _projectionSize, perPointDimRankings[i]);
            break;
        case filters::FilterType::HD_PEAK:
            _hdFloodPeakFilter.computeDimensionRanking(i, _dataMatrix, _variances, exportFloodFill, perPointDimRankings[i]);
            break;
        }
    }

    writeDimensionRanking(perPointDimRankings, _enabledDimNames);
}

/******************************************************************************
 * Flooding
 ******************************************************************************/

void ScatterplotPlugin::createKnnIndex()
{
    qDebug() << "Creating index";
    if (_dataMatrix.cols() <= 200)
        _knnIndex.create(_dataMatrix.cols(), knn::Metric::MANHATTAN);
    else
        _knnIndex.create(_dataMatrix.cols(), knn::Metric::EUCLIDEAN);
    qDebug() << "Adding data";
    _knnIndex.addData(_dataMatrix);
    qDebug() << "Done creating index";
}

void ScatterplotPlugin::computeKnnGraph()
{
    qDebug() << "Building KNN Graph..";
    if (!_preloadedKnnGraph)
        _largeKnnGraph.build(_dataMatrix, _knnIndex, 30);

    if (_useSharedDistances)
    {
        _sourceKnnGraph.build(_dataMatrix, _knnIndex, 100);
        _knnGraph.build(_sourceKnnGraph, 10, true);
    }
    else
        _knnGraph.build(_largeKnnGraph, 10);

    _graphAvailable = true;
    qDebug() << "Done building KNN Graph! Ready for flood-fill.";
}

void ScatterplotPlugin::exportFloodnodes()
{
    std::vector<std::vector<int>> perPointFloodNodes(_dataMatrix.rows());
    for (int p = 0; p < _dataMatrix.rows(); p++)
    {
        FloodFill exportFloodFill(_floodFill.getNumWaves());
        exportFloodFill.compute(_knnGraph, p);

        // Store all flood nodes together
        perPointFloodNodes[p].resize(exportFloodFill.getTotalNumNodes() + exportFloodFill.getNumWaves());
        int n = 0;
        for (int i = 0; i < exportFloodFill.getNumWaves(); i++)
        {
            perPointFloodNodes[p][n++] = -1;
            for (int j = 0; j < exportFloodFill.getWaves()[i].size(); j++)
            {
                perPointFloodNodes[p][n++] = exportFloodFill.getWaves()[i][j];
            }
        }
    }

    writeFloodNodes(perPointFloodNodes);
}

void ScatterplotPlugin::importKnnGraph()
{
    QString fileName = QFileDialog::getOpenFileName(&getWidget(),
        tr("Open Knn Graph"), "", tr("KNN Files (*.knn)"));

    _largeKnnGraph.readFromFile(fileName);
    _knnGraph.build(_largeKnnGraph, 10);

    _graphAvailable = true;
}

/******************************************************************************
 * Serialization
 ******************************************************************************/

void ScatterplotPlugin::fromVariantMap(const QVariantMap& variantMap)
{
    ViewPlugin::fromVariantMap(variantMap);

    variantMapMustContain(variantMap, "Settings");

    _settingsAction.fromVariantMap(variantMap["Settings"].toMap());

    _overlayType = static_cast<OverlayType>(variantMap["OverlayType"].toInt());

    int numPoints = static_cast<size_t>(variantMap["numPoints"].toInt());
    int numNeighbours = static_cast<size_t>(variantMap["numNeighbours"].toInt());
    const auto qneighbours = variantMap["largeKnnGraph"].toMap();

    std::vector<int> linearNeighbours(numPoints * numNeighbours);

    populateDataBufferFromVariantMap(qneighbours, (char*)linearNeighbours.data());

    std::vector<std::vector<int>> neighbours(numPoints, std::vector<int>(numNeighbours));
    int c = 0;
    for (int i = 0; i < neighbours.size(); i++)
    {
        for (int j = 0; j < neighbours[i].size(); j++)
        {
            neighbours[c / numNeighbours][c % numNeighbours] = linearNeighbours[c];
            c++;
        }
    }

    _largeKnnGraph._neighbours = neighbours;
    _largeKnnGraph._numNeighbours = numNeighbours;

    _preloadedKnnGraph = true;

    positionDatasetChanged();
}

QVariantMap ScatterplotPlugin::toVariantMap() const
{
    QVariantMap variantMap = ViewPlugin::toVariantMap();

    _settingsAction.insertIntoVariantMap(variantMap);

    variantMap.insert("OverlayType", static_cast<int>(_overlayType));

    const std::vector<std::vector<int>>& neighbours = _largeKnnGraph.getNeighbours();

    // Linearize data
    std::vector<int> linearNeighbours(neighbours.size() * neighbours[0].size());
    int c = 0;
    for (int i = 0; i < neighbours.size(); i++)
    {
        const std::vector<int>& n = neighbours[i];
        for (int j = 0; j < neighbours[i].size(); j++)
            linearNeighbours[c++] = neighbours[i][j];
    }

    QVariantMap qneighbours = rawDataToVariantMap((char*)linearNeighbours.data(), linearNeighbours.size() * sizeof(std::int32_t), true);

    variantMap.insert("largeKnnGraph", qneighbours);
    variantMap.insert("numPoints", neighbours.size());
    variantMap.insert("numNeighbours", neighbours[0].size());

    return variantMap;
}


/******************************************************************************
 * Mask
 ******************************************************************************/

bool ScatterplotPlugin::hasMaskApplied()
{
    return !_mask.empty();
}

void ScatterplotPlugin::clearMask()
{
    _mask.clear();

    // Set point opacity
    std::vector<float> opacityScalars(_dataMatrix.rows(), 1.0f);
    getScatterplotWidget().setPointOpacityScalars(opacityScalars);
}

void ScatterplotPlugin::useSelectionAsMask()
{
    // Get current selection
    // Compute the indices that are selected in this local dataset
    std::vector<uint32_t> localSelectionIndices;
    _positionDataset->getLocalSelectionIndices(localSelectionIndices);

    _mask.assign(localSelectionIndices.begin(), localSelectionIndices.end());

    // Set point opacity
    std::vector<float> opacityScalars(_dataMatrix.rows(), 0.2f);
    for (const int maskIndex : _mask)
        opacityScalars[maskIndex] = 1.0f;
    getScatterplotWidget().setPointOpacityScalars(opacityScalars);

    _maskedDataMatrix = _dataMatrix(_mask, Eigen::all);
    if (_maskedDataMatrix.rows() < 5000)
        _maskedKnnIndex.create(_maskedDataMatrix.cols(), knn::Metric::MANHATTAN);
    else
        _maskedKnnIndex.create(_maskedDataMatrix.cols(), knn::Metric::EUCLIDEAN);
    _maskedKnnIndex.addData(_maskedDataMatrix);

    _maskedProjMatrix = _projMatrix(_mask, Eigen::all);
    //_largeKnnGraph.build(_maskedDataMatrix, _maskedKnnIndex, 30);

    if (_maskedDataMatrix.rows() < 5000)
    {
        _maskedSourceKnnGraph.build(_maskedDataMatrix, _maskedKnnIndex, 100);
        _maskedKnnGraph.build(_maskedSourceKnnGraph, 10, true);
    }
    else
        _maskedKnnGraph.build(_maskedDataMatrix, _maskedKnnIndex, 10);
}

/******************************************************************************
 * Factory
 ******************************************************************************/

QIcon ScatterplotPluginFactory::getIcon(const QColor& color /*= Qt::black*/) const
{
    return Application::getIconFont("FontAwesome").getIcon("braille", color);
}

ViewPlugin* ScatterplotPluginFactory::produce()
{
    return new ScatterplotPlugin(this);
}

PluginTriggerActions ScatterplotPluginFactory::getPluginTriggerActions(const hdps::Datasets& datasets) const
{
    PluginTriggerActions pluginTriggerActions;

    const auto getInstance = [this]() -> ScatterplotPlugin* {
        return dynamic_cast<ScatterplotPlugin*>(Application::core()->getPluginManager().requestViewPlugin(getKind()));
    };

    const auto numberOfDatasets = datasets.count();

    return pluginTriggerActions;
}
