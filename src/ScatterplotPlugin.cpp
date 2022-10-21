#include "ScatterplotPlugin.h"
#include "ScatterplotWidget.h"
#include "ProjectionView.h"
#include "DataHierarchyItem.h"
#include "Application.h"

#include "util/PixelSelectionTool.h"

#include "PointData.h"
#include "ClusterData.h"
#include "ColorData.h"

#include "graphics/Vector2f.h"
//#include "graphics/Vector3f.h"
#include "widgets/DropWidget.h"

#include <Eigen/Dense>
#include "LocalDimensionality.h"
#include "RandomWalks.h"

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

Q_PLUGIN_METADATA(IID "nl.biovault.GradientExplorerPlugin")

using namespace hdps;
using namespace hdps::util;

namespace
{
    void convertToEigenMatrix(hdps::Dataset<Points> dataset, DataMatrix& dataMatrix)
    {
        int numPoints = dataset->getNumPoints();
        int numDimensions = dataset->getNumDimensions();

        dataMatrix.resize(numPoints, numDimensions);

        for (int i = 0; i < numPoints; i++)
        {
            for (int j = 0; j < numDimensions; j++)
            {
                int index = dataset->isFull() ? i * numDimensions + j : dataset->indices[i] * numDimensions + j;
                dataMatrix(i, j) = dataset->getValueAt(index);
            }
        }
    }

    void indicesToVectors(const std::vector<int>& indices, std::vector<Vector2f>& vec)
    {

    }

    void findPointsInRadius(Vector2f center, float radius, const std::vector<Vector2f>& points, std::vector<int>& indices)
    {
        float radiusSqr = powf(radius, 2);
        for (int i = 0; i < points.size(); i++)
        {
            const Vector2f& pos = points[i];

            Vector2f diff = center - pos;
            float len = diff.x * diff.x + diff.y * diff.y;

            if (len < radiusSqr)
            {
                indices.push_back(i);
            }
        }
    }
    
    void computeDimensionAverage(const DataMatrix& data, const std::vector<int>& indices, std::vector<float>& averages)
    {
        int numDimensions = data.cols();
        averages.resize(numDimensions);
        for (int d = 0; d < numDimensions; d++)
        {
            for (const int& index : indices)
            {
                float v = data(index, d);
                averages[d] += v;
            }
            averages[d] /= indices.size();
        }
    }

    void computeDirection(DataMatrix& dataMatrix, DataMatrix& projMatrix, KnnGraph& knnGraph, std::vector<Vector2f>& directions)
    {
        for (int p = 0; p < dataMatrix.rows(); p++)
        {
            std::vector<std::vector<int>> floodFill;
            compute::doFloodFill(dataMatrix, projMatrix, knnGraph, p, floodFill);

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
    _projectionViews(3),
    _dropWidget(nullptr),
    _settingsAction(this),
    _gradientGraph(new GradientGraph()),
    _selectedDimension(-1),
    _dimPicker(this)
{
    setObjectName("GradientExplorer");

    _dropWidget = new DropWidget(_scatterPlotWidget);

    for (int i = 0; i < _projectionViews.size(); i++)
    {
        _projectionViews[i] = new ProjectionView();
    }

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

    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    layout->addWidget(_settingsAction.createWidget(&getWidget()));
    gradientViewLayout->addWidget(_projectionViews[0], 33);
    gradientViewLayout->addWidget(_projectionViews[1], 33);
    gradientViewLayout->addWidget(_gradientGraph, 33);
    QPushButton* showRandomWalk = new QPushButton("Show random walks");
    connect(showRandomWalk, &QPushButton::pressed, &getScatterplotWidget(), &ScatterplotWidget::showRandomWalk);
    gradientViewLayout->addWidget(showRandomWalk);
    QPushButton* showDirections = new QPushButton("Show directions");
    connect(showDirections, &QPushButton::pressed, &getScatterplotWidget(), &ScatterplotWidget::showDirections);
    gradientViewLayout->addWidget(showDirections);
    QPushButton* showLocalDimensionality = new QPushButton("Show local dimensionality");
    connect(showLocalDimensionality, &QPushButton::pressed, this, &ScatterplotPlugin ::showLocalDimensionality);
    gradientViewLayout->addWidget(showLocalDimensionality);
    gradientViewLayout->addWidget(_dimPicker.createWidget(&getWidget()));
    QPushButton* updateFeatureSet = new QPushButton("Update feature set");
    connect(updateFeatureSet, &QPushButton::pressed, this, [this]() {
        computeStaticData();
    });
    gradientViewLayout->addWidget(updateFeatureSet);

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

    // Update the selection when the pixel selection tool selected area changed
    connect(&_scatterPlotWidget->getPixelSelectionTool(), &PixelSelectionTool::areaChanged, [this]() {
        if (_scatterPlotWidget->getPixelSelectionTool().isNotifyDuringSelection())
            selectPoints();
    });

    // Update the selection when the pixel selection process ended
    connect(&_scatterPlotWidget->getPixelSelectionTool(), &PixelSelectionTool::ended, [this]() {
        if (_scatterPlotWidget->getPixelSelectionTool().isNotifyDuringSelection())
            return;

        selectPoints();
    });

    _eventListener.setEventCore(Application::core());
    _eventListener.addSupportedEventType(static_cast<std::uint32_t>(EventType::DataSelectionChanged));
    _eventListener.registerDataEventByType(PointType, std::bind(&ScatterplotPlugin::onDataEvent, this, std::placeholders::_1));

    // Load points when the pointer to the position dataset changes
    connect(&_positionDataset, &Dataset<Points>::changed, this, &ScatterplotPlugin::positionDatasetChanged);

    // Update points when the position dataset data changes
    connect(&_positionDataset, &Dataset<Points>::dataChanged, this, &ScatterplotPlugin::updateData);

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
                // Extract the enabled dimensions from the data
                std::vector<bool> enabledDimensions = _dimPicker.getPickerAction().getEnabledDimensions();

                hdps::Dataset<Points> selection = _positionSourceDataset->getSelection();

                int numDimensions = _dataMatrix.cols();

                Bounds bounds = _scatterPlotWidget->getBounds();

                float size = bounds.getWidth() > bounds.getHeight() ? bounds.getWidth() : bounds.getHeight();
                
                if (selection->indices.size() > 0)
                {
                    int selectionIndex = selection->indices[0];
                    Vector2f center = _positions[selectionIndex];

                    ///////////////////////
                    // Do random walking //
                    ///////////////////////
                    //{
                    //    std::vector<std::vector<int>> randomWalks;
                    //    compute::doRandomWalksKNN(_dataMatrix, _projMatrix, _knnGraph, selectionIndex, randomWalks);

                    //    std::vector<std::vector<Vector2f>> randomWalkPoints;
                    //    randomWalkPoints.resize(randomWalks.size());
                    //    for (int i = 0; i < randomWalks.size(); i++)
                    //    {
                    //        randomWalkPoints[i].resize(randomWalks[i].size());
                    //        for (int j = 0; j < randomWalks[i].size(); j++)
                    //        {
                    //            int index = randomWalks[i][j];
                    //            randomWalkPoints[i][j] = Vector2f(_projMatrix(index, 0), _projMatrix(index, 1));
                    //        }
                    //    }

                    //    _scatterPlotWidget->setRandomWalks(randomWalkPoints);

                    //    //std::vector<Vector3f> colors(_dataMatrix.rows(), Vector3f(0, 0, 0));
                    //    //for (int i = 0; i < randomWalks.size(); i++)
                    //    //{
                    //    //    for (int j = 0; j < randomWalks[i].size(); j++)
                    //    //    {
                    //    //        int index = randomWalks[i][j];
                    //    //        colors[index] += Vector3f(0.1f, 0, 0);
                    //    //    }
                    //    //}

                    //    //_scatterPlotWidget->setColors(colors);
                    //}

                    //////////////////
                    // Do floodfill //
                    //////////////////
                    std::vector<std::vector<int>> floodFill;

                    {
                        compute::doFloodFill(_dataMatrix, _projMatrix, _knnGraph, selectionIndex, floodFill);

                        std::vector<float> ccolors(_dataMatrix.rows(), 0);
                        for (int i = 0; i < floodFill.size(); i++)
                        {
                            for (int j = 0; j < floodFill[i].size(); j++)
                            {
                                int index = floodFill[i][j];
                                ccolors[index] = 1 - 0.1 * i;
                            }
                        }

                        int numNodes = 0;
                        for (int i = 0; i < 3; i++)
                        {
                            numNodes += floodFill[i].size();
                        }

                        Eigen::MatrixXf nodeLocations(numNodes, 2);
                        int p = 0;
                        for (int i = 0; i < 3; i++)
                        {
                            for (int j = 0; j < floodFill[i].size(); j++)
                            {
                                int index = floodFill[i][j];

                                nodeLocations(p, 0) = _projMatrix(index, 0);
                                nodeLocations(p, 1) = _projMatrix(index, 1);
                                p++;
                            }
                        }
                        
                        getScatterplotWidget().setColorMap(_settingsAction.getColoringAction().getColorMapAction().getColorMapImage());
                        _scatterPlotWidget->setColoringMode(ScatterplotWidget::ColoringMode::Data);
                        getScatterplotWidget().setScalarEffect(PointEffect::Color);
                        _scatterPlotWidget->setScalars(ccolors);
                        _scatterPlotWidget->setColorMapRange(0, 1);

                        //getScatterplotWidget().setColorMap(_settingsAction.getColoringAction().getColorMapAction().getColorMapImage());
                        //_scatterPlotWidget->setColoringMode(ScatterplotWidget::ColoringMode::Data);
                        //getScatterplotWidget().setScalarEffect(PointEffect::Color);
                        //_scatterPlotWidget->setScalars(_colors);
                        //_scatterPlotWidget->setColorMapRange(0, 1);
                    }

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
                    {
                        for (int j = 0; j < floodFill[i].size(); j++)
                        {
                            floodNodes[n++] = floodFill[i][j];
                        }
                    }

                    //int currentNode = selectionIndex;
                    //Vector2f currentNodePos = center;

                    //// Find neighbours
                    //std::vector<int> neighbours;
                    //for (int i = 0; i < floodNodes.size(); i++)
                    //{
                    //    int floodIndex = floodNodes[i];
                    //    Vector2f nodePos = _positions[floodIndex];
                    //    if (abs(nodePos.x - currentNodePos.x) < 1.1f && abs(nodePos.y - currentNodePos.y) < 1.1f)
                    //    {
                    //        neighbours.push_back(floodIndex);
                    //    }
                    //}
                    //
                    //auto nodeValues = _dataMatrix.row(currentNode);

                    //for (int i = 0; i < neighbours.size(); i++)
                    //{
                    //    auto neighbourValues = _dataMatrix.row(neighbours[i]);


                    //}

                    //// Store dimension values for every flood node
                    //std::vector<std::vector<float>> dimValues(_dataMatrix.cols(), std::vector<float>(floodNodes.size()));
                    //for (int i = 0; i < floodNodes.size(); i++)
                    //{
                    //    int floodNode = floodNodes[i];
                    //    auto floodValues = _dataMatrix.row(floodNode);
                    //    for (int d = 0; d < floodValues.size(); d++)
                    //    {
                    //        dimValues[d][i] = floodValues(d);
                    //    }
                    //}

                    //for (int d = 0; d < dimValues.size(); d++)
                    //{
                    //    sort(dimValues[d].begin(), dimValues[d].end());
                    //}

                    //_gradientGraph->setValues(dimValues);

                    //std::vector<std::vector<Vector2f>> lineagePoints;

                    //auto centerValues = _dataMatrix.row(selectionIndex);
                    //for (int i = 1; i < 2; i++)
                    //{
                    //    for (int j = 0; j < floodFill[i].size(); j++)
                    //    {
                    //        int node = floodFill[i][j];
                    //        auto nodeValues = _dataMatrix.row(node);
                    //        std::vector<Vector2f> linPoints;
                    //        linPoints.push_back(Vector2f(_projMatrix(selectionIndex, 0), _projMatrix(selectionIndex, 1)));
                    //        linPoints.push_back(Vector2f(_projMatrix(node, 0), _projMatrix(node, 1)));
                    //        lineagePoints.push_back(linPoints);
                    //    }
                    //}
                    //_scatterPlotWidget->setRandomWalks(lineagePoints);

                    //////////////////
                    std::vector<std::vector<Vector2f>> linPoints(10, std::vector<Vector2f>());
                    for (int w = 0; w < 10; w++)
                    {
                        std::vector<int> lineage;
                        lineage.push_back(selectionIndex);
                        compute::traceLineage(_dataMatrix, floodFill, _positions, selectionIndex, lineage);

                        for (int i = 0; i < lineage.size(); i++)
                        {
                            linPoints[w].push_back(_positions[lineage[i]]);
                            //qDebug() << linPoints[w][i].x << linPoints[w][i].y;
                        }
                    }
                    _scatterPlotWidget->setRandomWalks(linPoints);

                    /////////////////////
                    // Gradient picker //
                    /////////////////////
                    // Small and large circle averages
                    std::vector<std::vector<float>> averages(2);
                    std::vector<std::vector<int>> circleIndices(2);

                    findPointsInRadius(center, 0.05f * size, _positions, circleIndices[0]);
                    computeDimensionAverage(_dataMatrix, circleIndices[0], averages[0]);
                    findPointsInRadius(center, 0.1f * size, _positions, circleIndices[1]);
                    computeDimensionAverage(_dataMatrix, circleIndices[1], averages[1]);

                    std::vector<float> diffAverages(numDimensions);
                    for (int d = 0; d < numDimensions; d++)
                    {
                        diffAverages[d] = fabs(averages[0][d] - averages[1][d]);
                    }
                    
                    // Sort averages from high to low
                    std::vector<size_t> idx(diffAverages.size());
                    std::iota(idx.begin(), idx.end(), 0);

                    std::stable_sort(idx.begin(), idx.end(), [&diffAverages](size_t i1, size_t i2) {return diffAverages[i1] > diffAverages[i2]; });

                    // Set appropriate coloring of gradient view, FIXME use colormap later
                    for (int pi = 0; pi < _projectionViews.size(); pi++)
                    {
                        auto dimValues = _dataMatrix(Eigen::all, idx[pi]);

                        float minV = *std::min_element(dimValues.begin(), dimValues.end());
                        float maxV = *std::max_element(dimValues.begin(), dimValues.end());
                        for (int i = 0; i < dimValues.size(); i++)
                        {
                            dimValues[i] /= maxV - minV;
                        }

                        std::vector<Vector3f> colors(_positions.size());
                        for (int i = 0; i < dimValues.size(); i++)
                        {
                            colors[i] = Vector3f(1-dimValues[i], 1 - dimValues[i], 1 - dimValues[i]);
                            if (i == selectionIndex)
                            {
                                colors[i] = Vector3f(1, 0, 0);
                            }
                        }
                        _projectionViews[pi]->setColors(colors);
                        const auto& dimNames = _positionSourceDataset->getDimensionNames();
                        auto enabledDimensions = _dimPicker.getPickerAction().getEnabledDimensions();
                        std::vector<QString> enabledDimNames;
                        for (int i = 0; i < enabledDimensions.size(); i++)
                        {
                            if (enabledDimensions[i])
                                enabledDimNames.push_back(dimNames[i]);
                        }

                        _projectionViews[pi]->setProjectionName(enabledDimNames[idx[pi]]);
                    }

                    //_gradientGraph->setDimension(_dataMatrix, idx[0]);

                    // Show dim values on floodnodes
                    std::vector<float> dimScalars(_dataMatrix.rows(), 0);

                    for (int i = 0; i < floodNodes.size(); i++)
                    {
                        dimScalars[floodNodes[i]] = _dataMatrix(floodNodes[i], idx[0]);
                    }
                    getScatterplotWidget().setScalars(dimScalars);
                }
            }
        }
    }
}

void ScatterplotPlugin::computeStaticData()
{
    DataMatrix dataMatrix;
    convertToEigenMatrix(_positionSourceDataset, dataMatrix);
    convertToEigenMatrix(_positionDataset, _projMatrix);

    // Subsample based on subset or full data
    if (!_positionDataset->isFull())
        _dataMatrix = dataMatrix(_positionDataset->indices, Eigen::all);
    else
        _dataMatrix = dataMatrix;

    // Subsample dimensions based on dimension picker
    std::vector<bool> enabledDimensions = _dimPicker.getPickerAction().getEnabledDimensions();
    std::vector<int> dimIndices;
    for (int i = 0; i < enabledDimensions.size(); i++)
        if (enabledDimensions[i]) dimIndices.push_back(i);
    DataMatrix featureMatrix = _dataMatrix(Eigen::all, dimIndices);
    _dataMatrix = featureMatrix;

    createKnnGraph(_dataMatrix);
    _knnGraph.build(_dataMatrix, _kdtree, 6);
    _largeKnnGraph.build(_dataMatrix, _kdtree, 30);

    //compute::computeSpatialLocalDimensionality(_dataMatrix, _projMatrix, _colors);
    compute::computeHDLocalDimensionality(_dataMatrix, _largeKnnGraph, _colors);
    getScatterplotWidget().setScalars(_colors);

    std::vector<Vector2f> directions;
    computeDirection(_dataMatrix, _projMatrix, _knnGraph, directions);
    getScatterplotWidget().setDirections(directions);
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

void ScatterplotPlugin::selectPoints()
{
    // Only proceed with a valid points position dataset and when the pixel selection tool is active
    if (!_positionDataset.isValid() || !_scatterPlotWidget->getPixelSelectionTool().isActive())
        return;

    // Get binary selection area image from the pixel selection tool
    auto selectionAreaImage = _scatterPlotWidget->getPixelSelectionTool().getAreaPixmap().toImage();

    // Get smart pointer to the position selection dataset
    auto selectionSet = _positionDataset->getSelection<Points>();

    // Create vector for target selection indices
    std::vector<std::uint32_t> targetSelectionIndices;

    // Reserve space for the indices
    targetSelectionIndices.reserve(_positionDataset->getNumPoints());

    // Mapping from local to global indices
    std::vector<std::uint32_t> localGlobalIndices;

    // Get global indices from the position dataset
    _positionDataset->getGlobalIndices(localGlobalIndices);

    const auto dataBounds   = _scatterPlotWidget->getBounds();
    const auto width        = selectionAreaImage.width();
    const auto height       = selectionAreaImage.height();
    const auto size         = width < height ? width : height;

    // Loop over all points and establish whether they are selected or not
    for (std::uint32_t i = 0; i < _positions.size(); i++) {
        const auto uvNormalized     = QPointF((_positions[i].x - dataBounds.getLeft()) / dataBounds.getWidth(), (dataBounds.getTop() - _positions[i].y) / dataBounds.getHeight());
        const auto uvOffset         = QPoint((selectionAreaImage.width() - size) / 2.0f, (selectionAreaImage.height() - size) / 2.0f);
        const auto uv               = uvOffset + QPoint(uvNormalized.x() * size, uvNormalized.y() * size);

        // Add point if the corresponding pixel selection is on
        if (selectionAreaImage.pixelColor(uv).alpha() > 0)
            targetSelectionIndices.push_back(localGlobalIndices[i]);
    }

    // Selection should be subtracted when the selection process was aborted by the user (e.g. by pressing the escape key)
    const auto selectionModifier = _scatterPlotWidget->getPixelSelectionTool().isAborted() ? PixelSelectionModifierType::Remove : _scatterPlotWidget->getPixelSelectionTool().getModifier();

    switch (selectionModifier)
    {
        case PixelSelectionModifierType::Replace:
            break;

        case PixelSelectionModifierType::Add:
        case PixelSelectionModifierType::Remove:
        {
            // Get reference to the indices of the selection set
            auto& selectionSetIndices = selectionSet->indices;

            // Create a set from the selection set indices
            QSet<std::uint32_t> set(selectionSetIndices.begin(), selectionSetIndices.end());

            switch (selectionModifier)
            {
                // Add points to the current selection
                case PixelSelectionModifierType::Add:
                {
                    // Add indices to the set 
                    for (const auto& targetIndex : targetSelectionIndices)
                        set.insert(targetIndex);

                    break;
                }

                // Remove points from the current selection
                case PixelSelectionModifierType::Remove:
                {
                    // Remove indices from the set 
                    for (const auto& targetIndex : targetSelectionIndices)
                        set.remove(targetIndex);

                    break;
                }

                default:
                    break;
            }

            // Convert set back to vector
            targetSelectionIndices = std::vector<std::uint32_t>(set.begin(), set.end());

            break;
        }

        default:
            break;
    }

    // Apply the selection indices
    _positionDataset->setSelectionIndices(targetSelectionIndices);

    // Notify others that the selection changed
    _core->notifyDatasetSelectionChanged(_positionDataset);
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

    // Enable pixel selection if the point positions dataset is valid
    _scatterPlotWidget->getPixelSelectionTool().setEnabled(_positionDataset.isValid());
    
    // Do not show the drop indicator if there is a valid point positions dataset
    _dropWidget->setShowDropIndicator(!_positionDataset.isValid());

    // Update positions data
    updateData();

    // Update the window title to reflect the position dataset change
    updateWindowTitle();

    _dimPicker.getPickerAction().setPointsDataset(_positionSourceDataset);
    _positionSourceDataset->addAction(_dimPicker);

    computeStaticData();
}

void ScatterplotPlugin::createKnnGraph(const DataMatrix& highDim)
{
    Eigen::MatrixXf tarray = highDim.transpose();
    _kdtree = new knncpp::KDTreeMinkowskiX<float, knncpp::ManhattenDistance<float>>(tarray, true);
    _kdtree->build();

    qDebug() << "Rows:" << tarray.rows() << tarray.cols();
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

        updateSelection();

        //getScatterplotWidget().setColors(colors);
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

void ScatterplotPlugin::showLocalDimensionality()
{
    getScatterplotWidget().setScalars(_colors);
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

        //_mousePressed = true;

        break;
    }

    case QEvent::MouseButtonRelease:
    {
        auto mouseEvent = static_cast<QMouseEvent*>(event);

        //_mousePressed = false;

        break;
    }

    case QEvent::MouseMove:
    {
        auto mouseEvent = static_cast<QMouseEvent*>(event);

        QPoint mousePos = QPoint(mouseEvent->x(), mouseEvent->y());

        _scatterPlotWidget->setCurrentPosition(mousePos);
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

            Vector2f diff = uv - Vector2f(mousePos.x(), mousePos.y());
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
