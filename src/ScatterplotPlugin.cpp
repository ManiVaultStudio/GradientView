#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"
#include "ProjectionView.h"
#include "DataHierarchyItem.h"
#include "Application.h"
#include "actions/GroupsAction.h"

#include "util/PixelSelectionTool.h"
#include "DimensionsPickerAction.h"

#include "PointData.h"
#include "ClusterData.h"
#include "ColorData.h"

#include "graphics/Vector2f.h"
#include "graphics/Vector3f.h"
#include "widgets/DropWidget.h"

#include <Eigen/Dense>
#include "LocalDimensionality.h"
#include "RandomWalks.h"
#include "Timer.h"

#include <QtCore>
#include <QApplication>
#include <QDebug>
#include <QMenu>
#include <QAction>
#include <QMetaType>

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
    void convertToEigenMatrix(hdps::Dataset<Points> dataset, DataMatrix& dataMatrix)
    {
        // Compute num points
        std::vector<bool> enabledDims = dataset->getDimensionsPickerAction().getEnabledDimensions();
        int numPoints = dataset->getNumPoints();
        int numDimensions = dataset->getNumDimensions();
        int numEnabledDims =  std::count(enabledDims.begin(), enabledDims.end(), true);

        dataMatrix.resize(numPoints, numEnabledDims);

        if (dataset->isFull())
        {
            int col = 0;
            for (int d = 0; d < numDimensions; d++)
            {
                if (!enabledDims[d]) continue;

                std::vector<float> dimData;
                dataset->extractDataForDimension(dimData, d);
                for (int i = 0; i < numPoints; i++)
                    dataMatrix(i, col) = dimData[i];

                col++;
            }
        }
    }

    void indicesToVectors(const std::vector<int>& indices, std::vector<Vector2f>& vec)
    {

    }

    void computeDirection(DataMatrix& dataMatrix, DataMatrix& projMatrix, KnnGraph& knnGraph, int numSteps, std::vector<Vector2f>& directions)
    {
        for (int p = 0; p < dataMatrix.rows(); p++)
        {
            std::vector<std::vector<int>> floodFill;
            compute::doFloodFill(dataMatrix, projMatrix, knnGraph, p, numSteps, floodFill);

            int depth = 3;
            int numNodes = 0;
            for (int i = 0; i < depth; i++)
            {
                numNodes += floodFill[i].size();
            }

            Eigen::MatrixXf nodeLocations(numNodes, 2);
            int n = 0;
            for (int i = 0; i < depth; i++)
            {
                for (int j = 0; j < floodFill[i].size(); j++)
                {
                    int index = floodFill[i][j];

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
    _projectionViews(2),
    _selectedView(),
    _dropWidget(nullptr),
    _settingsAction(this),
    _gradientGraph(new GradientGraph()),
    _selectedDimension(-1),
    _numFloodSteps(10),
    _filterType(filters::FilterType::SPATIAL_PEAK),
    _overlayType(OverlayType::NONE)
{
    setObjectName("GradientExplorer");

    _dropWidget = new DropWidget(_scatterPlotWidget);

    for (int i = 0; i < _projectionViews.size(); i++)
    {
        _projectionViews[i] = new ProjectionView();
    }
    _selectedView = new ProjectionView();

    connect(_gradientGraph, &GradientGraph::lineClicked, this, &ScatterplotPlugin::onLineClicked);

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
                });
            }
            else {
                if (_positionDataset != candidateDataset && candidateDataset->getNumDimensions() >= 2) {

                    // The number of points is equal, so offer the option to replace the existing points dataset
                    dropRegions << new DropWidget::DropRegion(this, "Point position", description, "map-marker-alt", true, [this, candidateDataset]() {
                        _positionDataset = candidateDataset;
                    });
                }

                if (candidateDataset->getNumPoints() == _positionDataset->getNumPoints()) {

                    // The number of points is equal, so offer the option to use the points dataset as source for points colors
                    dropRegions << new DropWidget::DropRegion(this, "Point color", QString("Colorize %1 points with %2").arg(_positionDataset->getGuiName(), candidateDataset->getGuiName()), "palette", true, [this, candidateDataset]() {
                        _settingsAction.getColoringAction().addColorDataset(candidateDataset);
                        _settingsAction.getColoringAction().setCurrentColorDataset(candidateDataset);
                    });

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

        // Cluster dataset is about to be dropped
        if (dataType == ClusterType) {

            // Get clusters dataset from the core
            auto candidateDataset  = _core->requestDataset<Clusters>(datasetId);
            
            // Establish drop region description
            const auto description  = QString("Color points by %1").arg(candidateDataset->getGuiName());

            // Only allow user to color by clusters when there is a positions dataset loaded
            if (_positionDataset.isValid()) {

                if (_settingsAction.getColoringAction().hasColorDataset(candidateDataset)) {

                    // The clusters dataset is already loaded
                    dropRegions << new DropWidget::DropRegion(this, "Color", description, "palette", true, [this, candidateDataset]() {
                        _settingsAction.getColoringAction().setCurrentColorDataset(candidateDataset);
                    });
                }
                else {

                    // Use the clusters set for points color
                    dropRegions << new DropWidget::DropRegion(this, "Color", description, "palette", true, [this, candidateDataset]() {
                        _settingsAction.getColoringAction().addColorDataset(candidateDataset);
                        _settingsAction.getColoringAction().setCurrentColorDataset(candidateDataset);
                    });
                }
            }
            else {

                // Only allow user to color by clusters when there is a positions dataset loaded
                dropRegions << new DropWidget::DropRegion(this, "No points data loaded", "Clusters can only be visualized in concert with points data", "exclamation-circle", false);
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
    auto plotLayout = new QHBoxLayout();
    auto gradientViewLayout = new QVBoxLayout();
    auto dimensionViewsLayout = new QHBoxLayout();

    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    layout->addWidget(_settingsAction.createWidget(&getWidget()));

    auto groupsAction = new GroupsAction(&getWidget());

    QLabel* filterLabel = new QLabel();
    filterLabel->setText("Spatial Peak Ranking");
    QFont font = filterLabel->font();
    font.setPointSize(font.pointSize() * 2);
    filterLabel->setFont(font);

    gradientViewLayout->addWidget(filterLabel);

    dimensionViewsLayout->addWidget(_projectionViews[0], 50);
    dimensionViewsLayout->addWidget(_projectionViews[1], 50);
    gradientViewLayout->addLayout(dimensionViewsLayout, 50);

    // Filter layout
    QHBoxLayout* filterLayout = new QHBoxLayout();
    QPushButton* spatialPeakFilter = new QPushButton("Spatial Peak Filter");
    QPushButton* hdPeakFilter = new QPushButton("HD Peak Filter");

    connect(spatialPeakFilter, &QPushButton::pressed, this, [this, filterLabel]() {
        setFilterType(filters::FilterType::SPATIAL_PEAK);
        filterLabel->setText("Spatial Peak Ranking");
    });
    connect(hdPeakFilter, &QPushButton::pressed, this, [this, filterLabel]() {
        setFilterType(filters::FilterType::HD_PEAK);
        filterLabel->setText("HD Peak Ranking");
    });

    //gradientViewLayout->addWidget(filtersLabel);
    filterLayout->addWidget(spatialPeakFilter);
    filterLayout->addWidget(hdPeakFilter);
    gradientViewLayout->addLayout(filterLayout);

    // Filter parameters
    auto filterGroupAction = new GroupAction(&getWidget());

    auto innerFilterRadius = new ScalarAction(this, "Spatial Inner Filter Radius", 1, 10, 2.5f, 2.5f);
    auto outerFilterRadius = new ScalarAction(this, "Spatial Outer Filter Radius", 2, 20, 5, 5);
    connect(innerFilterRadius, &ScalarAction::magnitudeChanged, [this](float mag) { _spatialPeakFilter.setInnerFilterRadius(mag * 0.01f); });
    connect(outerFilterRadius, &ScalarAction::magnitudeChanged, [this](float mag) { _spatialPeakFilter.setOuterFilterRadius(mag * 0.01f); });
    *filterGroupAction << *innerFilterRadius;
    *filterGroupAction << *outerFilterRadius;

    auto innerFilterSize = new IntegralAction(this, "HD Inner Filter Size", 1, 10, 5, 5);
    auto outerFilterSize = new IntegralAction(this, "HD Outer Filter Size", 2, 10, 10, 10);
    connect(innerFilterSize, &IntegralAction::valueChanged, [this](int value) { _hdFloodPeakFilter.setInnerFilterSize(value); });
    connect(outerFilterSize, &IntegralAction::valueChanged, [this](int value) { _hdFloodPeakFilter.setOuterFilterSize(value); });
    *filterGroupAction << *innerFilterSize;
    *filterGroupAction << *outerFilterSize;

    gradientViewLayout->addWidget(filterGroupAction->createWidget(&this->getWidget()), Qt::AlignLeft);

    // Dimension selection
    QLabel* dimensionSelectionLabel = new QLabel("Dimension Selection");
    dimensionSelectionLabel->setFont(font);
    gradientViewLayout->addWidget(dimensionSelectionLabel);
    gradientViewLayout->addWidget(_selectedView, 50);
    gradientViewLayout->addWidget(_gradientGraph, 70);

    gradientViewLayout->addWidget(groupsAction->createWidget(&getWidget()), 40);

    // Overlay options
    {
        auto overlayGroupAction = new GroupAction(&getWidget(), true);
        overlayGroupAction->setText("Flood Nodes Overlay");
        overlayGroupAction->setShowLabels(false);
        
        IntegralAction* floodDecimal = new IntegralAction(this, "Flood nodes", 10, 500, 10, 10);
        connect(floodDecimal, &IntegralAction::valueChanged, this, [this](int32_t value)
        {
            _knnGraph.build(_dataMatrix, _knnIndex, value);
        });
        *overlayGroupAction << *floodDecimal;

        IntegralAction* floodStepsAction = new IntegralAction(this, "Flood steps", 3, 20, 10, 10);
        connect(floodStepsAction, &IntegralAction::valueChanged, this, [this](int32_t value)
        {
            _numFloodSteps = value;
            onPointSelection();
        });
        *overlayGroupAction << *floodStepsAction;

        QVector<TriggersAction::Trigger> triggers;
        triggers << TriggersAction::Trigger("Flood Steps", "Color flood points by closeness to seed point in HD space");
        triggers << TriggersAction::Trigger("Top Dimension Values", "Color flood points by values of top ranked dimension");
        triggers << TriggersAction::Trigger("Local Dimensionality", "Color flood points by local intrinsic dimensionality");
        triggers << TriggersAction::Trigger("Directions", "Show major eigenvector directions over flood points");

        TriggersAction* overlayTriggers = new TriggersAction(overlayGroupAction, "Overlay Triggers", triggers);

        connect(overlayTriggers, &TriggersAction::triggered, this, [this](int32_t triggerIndex)
        {
            getScatterplotWidget().showDirections(false);
            switch (triggerIndex)
            {
            case 0: _overlayType = OverlayType::NONE; break;
            case 1: _overlayType = OverlayType::DIM_VALUES; break;
            case 2: _overlayType = OverlayType::LOCAL_DIMENSIONALITY; break;
            case 3: {_overlayType = OverlayType::DIRECTIONS; getScatterplotWidget().showDirections(true); break; }
            }
        });

        groupsAction->addGroupAction(overlayGroupAction);
    }

    // Export options
    {
        auto exportGroupAction = new GroupAction(&getWidget(), false);
        exportGroupAction->setText("Export");
        exportGroupAction->setShowLabels(false);

        QVector<TriggersAction::Trigger> triggers;
        triggers << TriggersAction::Trigger("Export rankings", "");
        triggers << TriggersAction::Trigger("Export flood nodes", "");

        TriggersAction* exportTriggers = new TriggersAction(exportGroupAction, "Export Triggers", triggers);

        connect(exportTriggers, &TriggersAction::triggered, this, [this](int32_t triggerIndex)
        {
            switch (triggerIndex)
            {
            case 0: // Export rankings
            {
                std::vector<std::vector<int>> perPointDimRankings(_dataMatrix.rows());
                for (int i = 0; i < _dataMatrix.rows(); i++)
                {
                    switch (_filterType)
                    {
                    case filters::FilterType::SPATIAL_PEAK:
                        _spatialPeakFilter.computeDimensionRanking(i, _dataMatrix, _variances, _projMatrix, _projectionSize, perPointDimRankings[i]);
                        break;
                    case filters::FilterType::HD_PEAK:
                        std::vector<std::vector<int>> floodFill;
                        compute::doFloodFill(_dataMatrix, _projMatrix, _knnGraph, i, _numFloodSteps, floodFill);
                        _hdFloodPeakFilter.computeDimensionRanking(i, _dataMatrix, _variances, floodFill, perPointDimRankings[i]);
                        break;
                    }
                }

                writeDimensionRanking(perPointDimRankings, _enabledDimNames);
                break;
            }
            case 1: // Export flood nodes
            {
                std::vector<std::vector<int>> perPointFloodNodes(_dataMatrix.rows());
                for (int p = 0; p < _dataMatrix.rows(); p++)
                {
                    std::vector<std::vector<int>> floodFill;
                    compute::doFloodFill(_dataMatrix, _projMatrix, _knnGraph, p, _numFloodSteps, floodFill);

                    int numFloodNodes = 0;
                    for (int i = 0; i < floodFill.size(); i++)
                    {
                        numFloodNodes += floodFill[i].size();
                    }

                    // Store all flood nodes together
                    perPointFloodNodes[p].resize(numFloodNodes);
                    int n = 0;
                    for (int i = 0; i < floodFill.size(); i++)
                    {
                        for (int j = 0; j < floodFill[i].size(); j++)
                        {
                            perPointFloodNodes[p][n++] = floodFill[i][j];
                        }
                    }
                }

                writeFloodNodes(perPointFloodNodes);
                break;
            }
            }
        });

        groupsAction->addGroupAction(exportGroupAction);
    }

    // Map overlay
    {
        auto mapOverlayGroupAction = new GroupAction(&getWidget(), true);
        mapOverlayGroupAction->setText("Map Overlay");
        mapOverlayGroupAction->setShowLabels(false);

        QVector<TriggersAction::Trigger> triggers;
        triggers << TriggersAction::Trigger("Show random walks", "");
        triggers << TriggersAction::Trigger("Show directions", "");
        triggers << TriggersAction::Trigger("Show Spatial Local Dimensionality", "");
        triggers << TriggersAction::Trigger("Show HD Local Dimensionality", "");

        TriggersAction* overlayTriggers = new TriggersAction(mapOverlayGroupAction, "Map Overlay Triggers", triggers);

        connect(overlayTriggers, &TriggersAction::triggered, this, [this](int32_t triggerIndex)
        {
            getScatterplotWidget().showDirections(false);
            switch (triggerIndex)
            {
            case 0: getScatterplotWidget().showRandomWalk(); break;
            case 1: getScatterplotWidget().showDirections(true); break;
            case 2:
            {
                if (_localSpatialDimensionality.empty())
                    compute::computeSpatialLocalDimensionality(_dataMatrix, _projMatrix, _localSpatialDimensionality);
                getScatterplotWidget().setScalars(_localSpatialDimensionality);
                break;
            }
            case 3:
            {
                if (_localHighDimensionality.empty())
                    compute::computeHDLocalDimensionality(_dataMatrix, _largeKnnGraph, _localHighDimensionality);
                getScatterplotWidget().setScalars(_localHighDimensionality);
                break;
            }
            }
        });

        groupsAction->addGroupAction(mapOverlayGroupAction);
    }

    plotLayout->addWidget(_scatterPlotWidget, 60);
    plotLayout->addLayout(gradientViewLayout, 30);
    layout->addLayout(plotLayout, 100);

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

    //// Update the selection when the pixel selection tool selected area changed
    //connect(&_scatterPlotWidget->getPixelSelectionTool(), &PixelSelectionTool::areaChanged, [this]() {
    //    if (_scatterPlotWidget->getPixelSelectionTool().isNotifyDuringSelection())
    //        selectPoints();
    //});

    //// Update the selection when the pixel selection process ended
    //connect(&_scatterPlotWidget->getPixelSelectionTool(), &PixelSelectionTool::ended, [this]() {
    //    if (_scatterPlotWidget->getPixelSelectionTool().isNotifyDuringSelection())
    //        return;

    //    selectPoints();
    //});

    _eventListener.setEventCore(Application::core());
    _eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DataSelectionChanged));
    _eventListener.registerDataEventByType(PointType, std::bind(&ScatterplotPlugin::onDataEvent, this, std::placeholders::_1));

    // Load points when the pointer to the position dataset changes
    connect(&_positionDataset, &Dataset<Points>::changed, this, &ScatterplotPlugin::positionDatasetChanged);

    // Update points when the position dataset data changes
    //connect(&_positionDataset, &Dataset<Points>::dataChanged, this, &ScatterplotPlugin::updateData);

    // Update point selection when the position dataset data changes
    connect(&_positionDataset, &Dataset<Points>::dataSelectionChanged, this, &ScatterplotPlugin::updateSelection);

    // Update the window title when the GUI name of the position dataset changes
    connect(&_positionDataset, &Dataset<Points>::dataGuiNameChanged, this, &ScatterplotPlugin::updateWindowTitle);

    // Do an initial update of the window title
    updateWindowTitle();
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

void ScatterplotPlugin::onPointSelection()
{
    Timer timer;

    hdps::Dataset<Points> selection = _positionSourceDataset->getSelection();

    int numDimensions = _dataMatrix.cols();

    if (selection->indices.size() > 0)
    {
timer.start();
        _selectedPoint = selection->indices[0];
        Vector2f center = _positions[_selectedPoint];

        getScatterplotWidget().setCurrentPosition(center);
        getScatterplotWidget().setFilterRadii(Vector2f(_spatialPeakFilter.getInnerFilterRadius() * _projectionSize, _spatialPeakFilter.getOuterFilterRadius() * _projectionSize));

        //////////////////
        // Do floodfill //
        //////////////////
        std::vector<std::vector<int>> floodFill;

        {
            compute::doFloodFill(_dataMatrix, _projMatrix, _knnGraph, _selectedPoint, _numFloodSteps, floodFill);

            getScatterplotWidget().setColorMap(_settingsAction.getColoringAction().getColorMapAction().getColorMapImage());
            _scatterPlotWidget->setColoringMode(ScatterplotWidget::ColoringMode::Data);
            getScatterplotWidget().setScalarEffect(PointEffect::Color);
        }
timer.mark("Floodfill");

        /////////////////////
        // Gradient picker //
        /////////////////////
        std::vector<int> dimRanking;
        switch (_filterType)
        {
        case filters::FilterType::SPATIAL_PEAK:
            _spatialPeakFilter.computeDimensionRanking(_selectedPoint, _dataMatrix, _variances, _projMatrix, _projectionSize, dimRanking);
            break;
        case filters::FilterType::HD_PEAK:
            _hdFloodPeakFilter.computeDimensionRanking(_selectedPoint, _dataMatrix, _variances, floodFill, dimRanking);
            break;
        }
timer.mark("Ranking");

        // Set appropriate coloring of gradient view, FIXME use colormap later
        for (int pi = 0; pi < _projectionViews.size(); pi++)
        {
            const auto dimValues = _dataMatrix(Eigen::all, dimRanking[pi]);
            _projectionViews[pi]->setScalars(dimValues, _selectedPoint);
            _projectionViews[pi]->setProjectionName(_enabledDimNames[dimRanking[pi]]);
        }
        // Set selected gradient view
        if (_selectedDimension > 0)
        {
            const auto dimValues = _dataMatrix(Eigen::all, _selectedDimension);
            _selectedView->setScalars(dimValues, _selectedPoint);
            _selectedView->setProjectionName(_enabledDimNames[_selectedDimension]);
        }

        _gradientGraph->setTopDimensions(dimRanking[0], dimRanking[1]);
timer.mark("Filter");

        /////////////////////
        // Trace lineage   //
        /////////////////////
        int numFloodNodes = 0;
        for (int i = 0; i < floodFill.size(); i++)
        {
            numFloodNodes += floodFill[i].size();
        }

        // Store all flood nodes together
        std::vector<int> floodNodes(numFloodNodes);
        int n = 0;
        for (int i = 0; i < floodFill.size(); i++)
            for (int j = 0; j < floodFill[i].size(); j++)
                floodNodes[n++] = floodFill[i][j];

        // Store dimension values for every flood node
        //DataMatrix subMatrix = _dataMatrix(floodNodes, Eigen::all);
        //std::vector<std::vector<float>> dimValues(subMatrix.data(), subMatrix.data() + subMatrix.rows() * subMatrix.cols());

//        _dimValues.resize(_dataMatrix.cols(), std::vector<float>(floodNodes.size()));
////#pragma omp parallel for
//        for (int i = 0; i < floodNodes.size(); i++)
//        {
//            int floodNode = floodNodes[i];
//            auto floodValues = _dataMatrix.row(floodNode);
//            for (int d = 0; d < _dataMatrix.cols(); d++)
//            {
//                _dimValues[d][i] = _dataD[d][floodNode];
//            }
//        }

        //for (int d = 0; d < _dimValues.size(); d++)
        //{
        //    sort(_dimValues[d].begin(), _dimValues[d].end());
        //}

        // Binning
#pragma omp parallel for
        for (int d = 0; d < _bins.size(); d++)
        {
            for (int i = 0; i < _bins[d].size(); i++)
            {
                _bins[d][i] = 0;
            }
            for (int i = 0; i < floodNodes.size(); i++)
            {
                const float& f = std::min(0.99999f, _normalizedData[d][floodNodes[i]]);
                _bins[d][(int) (f * 30)] += 1;
            }
        }

        _gradientGraph->setBins(_bins);
timer.mark("Graph");

        //////////////////
        std::vector<std::vector<Vector2f>> linPoints(10, std::vector<Vector2f>());
        for (int w = 0; w < 10; w++)
        {
            std::vector<int> lineage;
            lineage.push_back(_selectedPoint);
            compute::traceLineage(_dataMatrix, floodFill, _positions, _selectedPoint, lineage);

            for (int i = 0; i < lineage.size(); i++)
            {
                linPoints[w].push_back(_positions[lineage[i]]);
                //qDebug() << linPoints[w][i].x << linPoints[w][i].y;
            }
        }
        _scatterPlotWidget->setRandomWalks(linPoints);

        // Coloring
        switch (_overlayType)
        {
        case OverlayType::NONE:
        {
            std::vector<float> ccolors(_dataMatrix.rows(), 0);
            for (int i = 0; i < floodFill.size(); i++)
            {
                for (int j = 0; j < floodFill[i].size(); j++)
                {
                    int index = floodFill[i][j];
                    ccolors[index] = 1 - (1.0f/_numFloodSteps) * i;
                }
            }

            _scatterPlotWidget->setScalars(ccolors);
            _scatterPlotWidget->setColorMapRange(0, 1);
            break;
        }
        case OverlayType::DIM_VALUES:
        {
            std::vector<float> dimScalars(_dataMatrix.rows(), 0);

            for (int i = 0; i < floodNodes.size(); i++)
            {
                dimScalars[floodNodes[i]] = _dataMatrix(floodNodes[i], dimRanking[0]);
            }
            getScatterplotWidget().setScalars(dimScalars);
            break;
        }
        case OverlayType::LOCAL_DIMENSIONALITY:
        {
            std::vector<float> scalars(_dataMatrix.rows(), 0);

            if (_localHighDimensionality.empty()) break;

            for (int i = 0; i < floodNodes.size(); i++)
            {
                scalars[floodNodes[i]] = _localHighDimensionality[floodNodes[i]];
            }
            getScatterplotWidget().setScalars(scalars);
            break;
        }
        case OverlayType::DIRECTIONS:
        {
            std::vector<float> scalars(_dataMatrix.rows(), 0);
            std::vector<Vector2f> directions(floodNodes.size() * 2);

            for (int i = 0; i < floodNodes.size(); i++)
            {
                float idx = floodNodes[i];
                directions[i * 2 + 0] = _directions[idx * 2 + 0];
                directions[i * 2 + 1] = _directions[idx * 2 + 1];
                scalars[idx] = 0.5f;
            }
            getScatterplotWidget().setScalars(scalars);
            getScatterplotWidget().setDirections(directions);

            break;
        }
        }
timer.finish("Overlay");
    }
}

void ScatterplotPlugin::computeStaticData()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Start conversion" << std::endl;
    DataMatrix dataMatrix;
    convertToEigenMatrix(_positionSourceDataset, dataMatrix);
    convertToEigenMatrix(_positionDataset, _fullProjMatrix);

    int xDim = _settingsAction.getPositionAction().getDimensionX();
    int yDim = _settingsAction.getPositionAction().getDimensionY();
    _projMatrix = _fullProjMatrix(Eigen::all, std::vector<int> { xDim, yDim });

    // Subsample based on subset or full data
    if (!_positionDataset->isFull())
        _dataMatrix = dataMatrix(_positionDataset->indices, Eigen::all);
    else
        _dataMatrix = dataMatrix;

    std::cout << "Number of enabled dimensions in the dataset : " << _dataMatrix.cols() << std::endl;
    std::chrono::duration<double> elapsedSub = std::chrono::high_resolution_clock::now() - start;
    start = std::chrono::high_resolution_clock::now();

    // Standardize
    // Compute means
    std::vector<float> means(_dataMatrix.cols());
    for (int d = 0; d < _dataMatrix.cols(); d++)
    {
        means[d] = 0;
        for (int i = 0; i < _dataMatrix.rows(); i++)
        {
            means[d] += _dataMatrix(i, d);
        }
        means[d] /= _dataMatrix.rows();
        std::cout << "Mean " << d << " " << means[d] << std::endl;
    }

    // Compute variances
    _variances.resize(_dataMatrix.cols());
    for (int d = 0; d < _dataMatrix.cols(); d++)
    {
        _variances[d] = 0;
        for (int i = 0; i < _dataMatrix.rows(); i++)
        {
            _variances[d] += (_dataMatrix(i, d) - means[d]) * (_dataMatrix(i, d) - means[d]);
        }
        _variances[d] /= _dataMatrix.rows();
        std::cout << "Variance " << d << " " << _variances[d] << std::endl;
    }
    for (int d = 0; d < _dataMatrix.cols(); d++)
    {
        if (_variances[d] <= 0) continue;
        for (int i = 0; i < _dataMatrix.rows(); i++)
        {
            _dataMatrix(i, d) -= means[d];
            _dataMatrix(i, d) /= sqrt(_variances[d]);
        }
    }

    // Compute normalized data
    _dataD.resize(_dataMatrix.cols(), std::vector<float>(_dataMatrix.rows()));
    _normalizedData.resize(_dataMatrix.cols(), std::vector<float>(_dataMatrix.rows()));
    for (int d = 0; d < _dataMatrix.cols(); d++)
    {
        for (int i = 0; i < _dataMatrix.rows(); i++)
            _dataD[d][i] = _dataMatrix(i, d);

        float minVal = *std::min_element(_dataD[d].begin(), _dataD[d].end());
        float maxVal = *std::max_element(_dataD[d].begin(), _dataD[d].end());
        float range = maxVal - minVal;
        if (range == 0) range = 1;

        for (int i = 0; i < _dataMatrix.rows(); i++)
            _normalizedData[d][i] = (_dataD[d][i] - minVal) / range;
    }
    _bins.resize(_dataMatrix.cols(), std::vector<int>(30));

    Bounds bounds = _scatterPlotWidget->getBounds();
    _projectionSize = bounds.getWidth() > bounds.getHeight() ? bounds.getWidth() : bounds.getHeight();
    std::cout << "Projection size: " << _projectionSize << std::endl;

    _knnIndex.create(_dataMatrix.cols(), knn::Metric::EUCLIDEAN);
    _knnIndex.addData(_dataMatrix);

    _knnGraph.build(_dataMatrix, _knnIndex, 10);
    _largeKnnGraph.build(_dataMatrix, _knnIndex, 10);

    std::chrono::duration<double> elapsedKnn = std::chrono::high_resolution_clock::now() - start;
    start = std::chrono::high_resolution_clock::now();
    


    _localSpatialDimensionality.clear();
    _localHighDimensionality.clear();

    std::chrono::duration<double> elapsedDim = std::chrono::high_resolution_clock::now() - start;
    start = std::chrono::high_resolution_clock::now();

    //computeDirection(_dataMatrix, _projMatrix, _knnGraph, _directions);
    //getScatterplotWidget().setDirections(_directions);

    std::chrono::duration<double> elapsedDir = std::chrono::high_resolution_clock::now() - start;
    start = std::chrono::high_resolution_clock::now();

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
    _gradientGraph->setNumDimensions(enabledDimensions.size());

    std::chrono::duration<double> elapsedGraph = std::chrono::high_resolution_clock::now() - start;
    start = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsedStan = std::chrono::high_resolution_clock::now() - start;
    start = std::chrono::high_resolution_clock::now();

    std::cout << "Elapsed time subsetting: " << elapsedSub.count() << " s\n";
    std::cout << "Elapsed time knn: " << elapsedKnn.count() << " s\n";
    std::cout << "Elapsed time local dimensionality: " << elapsedDim.count() << " s\n";
    std::cout << "Elapsed time directions: " << elapsedDir.count() << " s\n";
    std::cout << "Elapsed time graph setup: " << elapsedGraph.count() << " s\n";
    std::cout << "Elapsed time standardization: " << elapsedStan.count() << " s\n";
}

void ScatterplotPlugin::loadData(const Datasets& datasets)
{
    // Exit if there is nothing to load
    if (datasets.isEmpty())
        return;

    // Load the first dataset
    _positionDataset = datasets.first();

    // And set the coloring mode to constant
    _settingsAction.getColoringAction().getColorByAction().setCurrentIndex(1);
}

void ScatterplotPlugin::createSubset(const bool& fromSourceData /*= false*/, const QString& name /*= ""*/)
{
    auto subsetPoints = fromSourceData ? _positionDataset->getSourceDataset<Points>() : _positionDataset;

    // Create the subset
    auto subset = subsetPoints->createSubsetFromSelection(_positionDataset->getGuiName(), _positionDataset);

    // Notify others that the subset was added
    _core->notifyDatasetAdded(subset);
    
    // And select the subset
    subset->getDataHierarchyItem().select();
}

void ScatterplotPlugin::updateWindowTitle()
{
    if (!_positionDataset.isValid())
        getWidget().setWindowTitle(getGuiName());
    else
        getWidget().setWindowTitle(QString("%1: %2").arg(getGuiName(), _positionDataset->getDataHierarchyItem().getFullPathName()));
}

Dataset<Points>& ScatterplotPlugin::getPositionDataset()
{
    return _positionDataset;
}

Dataset<Points>& ScatterplotPlugin::getPositionSourceDataset()
{
    return _positionSourceDataset;
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

    // Update positions data
    updateData();

    // Update the window title to reflect the position dataset change
    updateWindowTitle();

    computeStaticData();
}

void ScatterplotPlugin::loadColors(const Dataset<Points>& points, const std::uint32_t& dimensionIndex)
{
    // Only proceed with valid points dataset
    if (!points.isValid())
        return;

    // Generate point scalars for color mapping
    std::vector<float> scalars;

    if (_positionDataset->getNumPoints() != _numPoints)
    {
        qWarning("Number of points used for coloring does not match number of points in data, aborting attempt to color plot");
        return;
    }

    // Populate point scalars
    if (dimensionIndex >= 0)
        points->extractDataForDimension(scalars, dimensionIndex);

    // Assign scalars and scalar effect
    std::vector<float> tempScalars(scalars.size(), 0);
    _scatterPlotWidget->setScalars(tempScalars);
    _scatterPlotWidget->setScalarEffect(PointEffect::Color);

    // Render
    getWidget().update();
}

void ScatterplotPlugin::loadColors(const Dataset<Clusters>& clusters)
{
    // Only proceed with valid clusters and position dataset
    if (!clusters.isValid() || !_positionDataset.isValid())
        return;

    // Mapping from local to global indices
    std::vector<std::uint32_t> globalIndices;

    // Get global indices from the position dataset
    _positionDataset->getGlobalIndices(globalIndices);

    // Generate color buffer for global and local colors
    std::vector<Vector3f> globalColors(globalIndices.back() + 1);
    std::vector<Vector3f> localColors(_positions.size());

    // Loop over all clusters and populate global colors
    for (const auto& cluster : clusters->getClusters())
        for (const auto& index : cluster.getIndices())
            globalColors[globalIndices[index]] = Vector3f(cluster.getColor().redF(), cluster.getColor().greenF(), cluster.getColor().blueF());

    std::int32_t localColorIndex = 0;

    // Loop over all global indices and find the corresponding local color
    for (const auto& globalIndex : globalIndices)
        localColors[localColorIndex++] = globalColors[globalIndex];

    // Apply colors to scatter plot widget without modification
    _scatterPlotWidget->setColors(localColors);

    // Render
    getWidget().update();
}

ScatterplotWidget& ScatterplotPlugin::getScatterplotWidget()
{
    return *_scatterPlotWidget;
}

void ScatterplotPlugin::updateData()
{
    // Check if the scatter plot is initialized, if not, don't do anything
    if (!_scatterPlotWidget->isInitialized())
        return;
    
    // If no dataset has been selected, don't do anything
    if (_positionDataset.isValid()) {

        // Get the selected dimensions to use as X and Y dimension in the plot
        const auto xDim = _settingsAction.getPositionAction().getDimensionX();
        const auto yDim = _settingsAction.getPositionAction().getDimensionY();

        // If one of the dimensions was not set, do not draw anything
        if (xDim < 0 || yDim < 0)
            return;

        // Determine number of points depending on if its a full dataset or a subset
        _numPoints = _positionDataset->getNumPoints();

        // Extract 2-dimensional points from the data set based on the selected dimensions
        calculatePositions(*_positionDataset);

        // Pass the 2D points to the scatter plot widget
        _scatterPlotWidget->setData(&_positions);
        for (int i = 0; i < _projectionViews.size(); i++)
            _projectionViews[i]->setData(&_positions);
        _selectedView->setData(&_positions);

        updateSelection();

        // Update proj matrix
        _projMatrix = _fullProjMatrix(Eigen::all, std::vector<int> { xDim, yDim });

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

void ScatterplotPlugin::setFilterType(filters::FilterType type)
{
    _filterType = type;
}

std::uint32_t ScatterplotPlugin::getNumberOfPoints() const
{
    if (!_positionDataset.isValid())
        return 0;

    return _positionDataset->getNumPoints();
}

void ScatterplotPlugin::setXDimension(const std::int32_t& dimensionIndex)
{
    updateData();
}

void ScatterplotPlugin::setYDimension(const std::int32_t& dimensionIndex)
{
    updateData();
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

        _mousePos = QPoint(mouseEvent->x(), mouseEvent->y());

        hdps::Bounds bounds = _scatterPlotWidget->getBounds();

        if (!_positionDataset.isValid())
            return QObject::eventFilter(target, event);

        const auto w = _scatterPlotWidget->width();
        const auto h = _scatterPlotWidget->height();
        const auto size = w < h ? w : h;

        // Loop over all points and establish whether they are selected or not
        int closestIndex = 0;
        float minDist = std::numeric_limits<float>::max();
        for (std::uint32_t i = 0; i < _positions.size(); i++) {
            const auto pointUV = Vector2f((_positions[i].x - bounds.getLeft()) / bounds.getWidth(), (bounds.getTop() - _positions[i].y) / bounds.getHeight());
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

        // Create vector for target selection indices
        std::vector<std::uint32_t> targetSelectionIndices;
        targetSelectionIndices.push_back(closestIndex);

        // Apply the selection indices
        _positionDataset->setSelectionIndices(targetSelectionIndices);

        // Notify others that the selection changed
        _core->notifyDatasetSelectionChanged(_positionDataset);

        //_explanationWidget->update();

        break;
    }

    default:
        break;
    }

    return QObject::eventFilter(target, event);
}

void ScatterplotPlugin::onLineClicked(int dim)
{
    qDebug() << "Dim: " << dim;
    _selectedDimension = dim;

    const auto dimValues = _dataMatrix(Eigen::all, _selectedDimension);

    float minV = *std::min_element(dimValues.begin(), dimValues.end());
    float maxV = *std::max_element(dimValues.begin(), dimValues.end());

    std::vector<Vector3f> colors(_positions.size());
    for (int i = 0; i < dimValues.size(); i++)
    {
        float dimValue = dimValues[i] / (maxV - minV);

        colors[i] = (i == _selectedPoint) ? Vector3f(1, 0, 0) : Vector3f(1 - dimValue, 1 - dimValue, 1 - dimValue);
    }
    _selectedView->setColors(colors);

    _selectedView->setProjectionName(_enabledDimNames[_selectedDimension]);
}

QIcon ScatterplotPluginFactory::getIcon(const QColor& color /*= Qt::black*/) const
{
    return Application::getIconFont("FontAwesome").getIcon("braille", color);
}

ViewPlugin* ScatterplotPluginFactory::produce()
{
    return new ScatterplotPlugin(this);
}

hdps::DataTypes ScatterplotPluginFactory::supportedDataTypes() const
{
    DataTypes supportedTypes;
    supportedTypes.append(PointType);
    return supportedTypes;
}
