#pragma once

#include <ViewPlugin.h>

#include "util/PixelSelectionTool.h"

#include "Common.h"

#include "SettingsAction.h"
#include "graphics/Vector3f.h"
#include "GradientGraph.h"
#include "Filters.h"

#include <Eigen/Eigen>
#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>

#include "KnnIndex.h"
#include "KnnGraph.h"
#include "Filters.h"

#include <QPoint>

using namespace hdps::plugin;
using namespace hdps::util;

using DataMatrix = Eigen::MatrixXf;

class Points;

class ScatterplotWidget;
class ProjectionView;

namespace hdps
{
    class CoreInterface;
    class Vector2f;

    namespace gui {
        class DropWidget;
    }
}

enum class OverlayType
{
    NONE,
    DIM_VALUES,
    LOCAL_DIMENSIONALITY,
    DIRECTIONS
};

class ScatterplotPlugin : public ViewPlugin
{
    Q_OBJECT
    
public:
    ScatterplotPlugin(const PluginFactory* factory);
    ~ScatterplotPlugin() override;

    void init() override;

    void onDataEvent(hdps::DataEvent* dataEvent);
    void onPointSelection();

    void computeStaticData();

    /**
     * Load one (or more datasets in the view)
     * @param datasets Dataset(s) to load
     */
    void loadData(const Datasets& datasets) override;

    /** Get number of points in the position dataset */
    std::uint32_t getNumberOfPoints() const;

public:
    void createSubset(const bool& fromSourceData = false, const QString& name = "");

public: // Dimension picking
    void setXDimension(const std::int32_t& dimensionIndex);
    void setYDimension(const std::int32_t& dimensionIndex);

protected: // Data loading

    /** Invoked when the position points dataset changes */
    void positionDatasetChanged();

    void createBruteGraph(const DataMatrix& highDim);
    void createKnnGraph(const DataMatrix& highDim);
    //void createAnnGraph(const DataMatrix& highDim);
    void createFaissGraph(const DataMatrix& highDim);
    void createFaissIVFGraph(const DataMatrix& highDim);
    void createFaissGpuGraph(const DataMatrix& highDim);
    void createAnnoyIndex(const DataMatrix& highDim);

public: // Point colors

    /**
     * Load color from points dataset
     * @param points Smart pointer to points dataset
     * @param dimensionIndex Index of the dimension to load
     */
    void loadColors(const Dataset<Points>& points, const std::uint32_t& dimensionIndex);

    /**
     * Load color from clusters dataset
     * @param clusters Smart pointer to clusters dataset
     */
    void loadColors(const Dataset<Clusters>& clusters);

public: // Miscellaneous

    /** Get smart pointer to points dataset for point position */
    Dataset<Points>& getPositionDataset();

    /** Get smart pointer to source of the points dataset for point position (if any) */
    Dataset<Points>& getPositionSourceDataset();

protected:

    /** Updates the window title (displays the name of the view and the GUI name of the loaded points dataset) */
    void updateWindowTitle();

public:

    /** Get reference to the scatter plot widget */
    ScatterplotWidget& getScatterplotWidget();
    std::vector<ProjectionView*>& getProjectionViews() { return _projectionViews; }

    SettingsAction& getSettingsAction() { return _settingsAction; }

private:
    void updateData();
    void calculatePositions(const Points& points);
    void updateSelection();
    void setFilterType(filters::FilterType type);

    bool eventFilter(QObject* target, QEvent* event);

private slots:
    void onLineClicked(int dim);

private:
    Dataset<Points>                 _positionDataset;           /** Smart pointer to points dataset for point position */
    Dataset<Points>                 _positionSourceDataset;     /** Smart pointer to source of the points dataset for point position (if any) */
    std::vector<hdps::Vector2f>     _positions;                 /** Point positions */
    unsigned int                    _numPoints;                 /** Number of point positions */

    DataMatrix                      _dataMatrix;
    std::vector<std::vector<float>> _dimValues;
    std::vector<std::vector<float>> _normalizedData;
    std::vector<std::vector<int>>   _bins;

    knn::Index                      _knnIndex;

    DataMatrix                      _fullProjMatrix;
    DataMatrix                      _projMatrix;
    KnnGraph                        _knnGraph;
    KnnGraph                        _largeKnnGraph;
    KnnGraph                        _sourceKnnGraph;
    std::vector<QString>            _enabledDimNames;
    std::vector<float>              _variances;

    int _numFloodSteps;

protected:
    ScatterplotWidget*              _scatterPlotWidget;
    std::vector<ProjectionView*>    _projectionViews;
    ProjectionView*                 _selectedView;
    std::vector<Vector3f> colors;
    std::vector<float>              _localSpatialDimensionality;
    std::vector<float>              _localHighDimensionality;
    std::vector<Vector2f>           _directions;
    GradientGraph*                  _gradientGraph;
    int                             _selectedPoint;
    int                             _selectedDimension;
    float                           _projectionSize;
    filters::FilterType             _filterType;
    filters::SpatialPeakFilter      _spatialPeakFilter;
    filters::HDFloodPeakFilter      _hdFloodPeakFilter;
    bool                            _useSharedDistances = false;
    OverlayType                     _overlayType;
    QPoint                          _mousePos;
    bool                            _mousePressed = false;

    hdps::gui::DropWidget*      _dropWidget;
    SettingsAction              _settingsAction;
};

// =============================================================================
// Factory
// =============================================================================

class ScatterplotPluginFactory : public ViewPluginFactory
{
    Q_INTERFACES(hdps::plugin::ViewPluginFactory hdps::plugin::PluginFactory)
    Q_OBJECT
    Q_PLUGIN_METADATA(IID   "nl.biovault.GradientExplorerPlugin"
                      FILE  "GradientExplorerPlugin.json")
    
public:
    ScatterplotPluginFactory(void) {}
    ~ScatterplotPluginFactory(void) override {}

    /**
     * Get plugin icon
     * @param color Icon color for flat (font) icons
     * @return Icon
     */
    QIcon getIcon(const QColor& color = Qt::black) const override;

    ViewPlugin* produce() override;

    hdps::DataTypes supportedDataTypes() const override;
};
