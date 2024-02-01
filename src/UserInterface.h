#pragma once

#include "ScatterplotWidget.h"
#include "ProjectionView.h"

#include <Widgets/DropWidget.h>
#include <Actions/SettingsAction.h>
#include <Actions/HorizontalToolbarAction.h>
#include <Actions/ColorMap1DAction.h>

class UserInterface
{
public:
    UserInterface(QObject* parent);

    void init();
    void initializeLabels();

public:
    ScatterplotWidget&              getMainView()                   { return *_scatterPlotWidget; }
    std::vector<ProjectionView*>&   getProjectionViews()            { return _projectionViews; }
    ProjectionView&                 getSelectionView()              { return *_selectedView; }

    SettingsAction&                 getSettingsAction()             { return _settingsAction; }
    const SettingsAction&           getSettingsAction() const       { return _settingsAction; }
    ColorMap1DAction&               getColorMapAction()             { return _colorMapAction; }
    DropWidget&                     getDropWidget()                 { return *_dropWidget; }
    HorizontalToolbarAction&        getPrimaryToolbar()             { return _primaryToolbarAction; }
    HorizontalToolbarAction&        getSecondaryToolbar()           { return _secondaryToolbarAction; }

    QLabel&                         getFilterLabel()                { return *_filterLabel; }
    QLabel&                         getDimensionSelectionLabel()    { return *_dimensionSelectionLabel; }
    QLabel&                         getSortedExpressionGraphLabel() { return *_sortedExpressionGraphLabel; }

private:
    ScatterplotWidget*              _scatterPlotWidget;
    std::vector<ProjectionView*>    _projectionViews;
    ProjectionView*                 _selectedView;

    mv::gui::DropWidget*            _dropWidget;
    SettingsAction                  _settingsAction;
    ColorMap1DAction                _colorMapAction;            /** Color map action */
    HorizontalToolbarAction         _primaryToolbarAction;      /** Horizontal toolbar for primary content */
    HorizontalToolbarAction         _secondaryToolbarAction;    /** Secondary toolbar for secondary content */

    // Labels
    QLabel*                         _filterLabel;
    QLabel*                         _dimensionSelectionLabel;
    QLabel*                         _sortedExpressionGraphLabel;
};
