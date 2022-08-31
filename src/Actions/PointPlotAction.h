#pragma once

#include "PluginAction.h"

#include "ScalarAction.h"

#include <QLabel>

using namespace hdps::gui;

class PointPlotAction : public PluginAction
{
protected: // Widget

    class Widget : public WidgetActionWidget {
    public:
        Widget(QWidget* parent, PointPlotAction* pointPlotAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this, widgetFlags);
    };

public:
    PointPlotAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

    /**
     * Add point size dataset
     * @param pointSizeDataset Smart pointer to point size dataset
     */
    void addPointSizeDataset(const Dataset<DatasetImpl>& pointSizeDataset);

    /**
     * Add point opacity dataset
     * @param pointOpacityDataset Smart pointer to point opacity dataset
     */
    void addPointOpacityDataset(const Dataset<DatasetImpl>& pointOpacityDataset);

protected:

    /** Update default datasets (candidates are children of points type and with matching number of points) */
    void updateDefaultDatasets();

    /** Update the scatter plot widget point size scalars */
    void updateScatterPlotWidgetPointSizeScalars();

    /** Update the scatter plot widget point opacity scalars */
    void updateScatterPlotWidgetPointOpacityScalars();

public: // Action getters

    ScalarAction& getSizeAction() { return _sizeAction; }
    ScalarAction& getOpacityAction() { return _opacityAction; }
    ToggleAction& getFocusSelection() { return _focusSelection; }

protected:
    ScalarAction            _sizeAction;                /** Point size action */
    ScalarAction            _opacityAction;             /** Point opacity action */
    std::vector<float>      _pointSizeScalars;          /** Cached point size scalars */
    std::vector<float>      _pointOpacityScalars;       /** Cached point opacity scalars */
    ToggleAction            _focusSelection;            /** Focus selection action */
    std::int32_t            _lastOpacitySourceIndex;    /** Last opacity source index that was selected */

    static constexpr double DEFAULT_POINT_SIZE      = 10.0;     /** Default point size */
    static constexpr double DEFAULT_POINT_OPACITY   = 50.0;     /** Default point opacity */

    friend class PlotAction;
};