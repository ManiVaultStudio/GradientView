#pragma once

#include <actions/GroupAction.h>
#include <actions/IntegralAction.h>
#include <actions/ToggleAction.h>
#include <actions/TriggerAction.h>
#include <actions/TriggersAction.h>
#include <actions/ColorAction.h>
#include <actions/StringAction.h>
#include <actions/StatusAction.h>

#include <PointData/DimensionsPickerAction.h>

using namespace hdps::gui;

class ScatterplotPlugin;

/**
 * Export image action class
 *
 * Action for export to image
 *
 * @author Thomas Kroes
 */
class ExportImageAction : public GroupAction
{
public:

    /** Scale options */
    enum Scale {
        Eighth,         /** Scale by one-eight */
        Quarter,        /** Scale by a quarter */
        Half,           /** Scale by a half */
        One,            /** Do not scale */
        Twice,          /** Scale by a factor of two */
        Thrice,         /** Scale by a factor of three */
        Four,           /** Scale by a factor of four */
        Eight           /** Scale by a factor of eight */
    };

    static const QMap<Scale, TriggersAction::Trigger> triggers;     /** Maps scale enum to trigger */
    static const QMap<Scale, float> scaleFactors;                   /** Maps scale enum to scale factor */

public:

    /**
     * Constructor
     * @param parent Pointer to parent object
     * @param scatterplotPlugin Reference to the scatter plot plugin
     */
    ExportImageAction(QObject* parent, ScatterplotPlugin& scatterplotPlugin);

    /** Grab target size from scatter plot widget */
    void initializeTargetSize();

    /** Export images to disk */
    void exportImages();

protected:

    /** Update the input points dataset of the dimensions picker action */
    void updateDimensionsPickerAction();

    /**
     * Establish whether export can take place
     * @return Whether export can take place
     */
    bool mayExport() const;

    /**
     * Get the number of selected dimensions
     * @return Number of selected dimensions
     */
    std::int32_t getNumberOfSelectedDimensions() const;

    /** Updates the export trigger text, tooltip and read-only */
    void updateExportTrigger();

public: // Action getters

    DimensionsPickerAction& getDimensionsPickerAction() { return _dimensionSelectionAction; }
    IntegralAction& getTargetWidthAction() { return _targetWidthAction; }
    IntegralAction& getTargetHeightAction() { return _targetHeightAction; }
    ToggleAction& getLockAspectRatioAction() { return _lockAspectRatioAction; }
    TriggersAction& getScaleAction() { return _scaleAction; }
    ColorAction& getBackgroundColorAction() { return _backgroundColorAction; }
    ToggleAction& getOverrideRangesAction() { return _overrideRangesAction; }
    DecimalRangeAction& getFixedRangeAction() { return _fixedRangeAction; }
    DirectoryPickerAction& getDirectoryPickerAction() { return _outputDirectoryAction; }
    TriggersAction& getExportCancelAction() { return _exportCancelAction; }

protected:
    ScatterplotPlugin&          _scatterplotPlugin;             /** Reference to scatterplot plugin */
    DimensionsPickerAction      _dimensionSelectionAction;      /** Dimension selection picker action */
    IntegralAction              _targetWidthAction;             /** Screenshot target width action */
    IntegralAction              _targetHeightAction;            /** Screenshot target height action */
    ToggleAction                _lockAspectRatioAction;         /** Lock aspect ratio action */
    TriggersAction              _scaleAction;                   /** Scale action */
    ColorAction                 _backgroundColorAction;         /** Background color action */
    ToggleAction                _overrideRangesAction;          /** Override ranges action */
    DecimalRangeAction          _fixedRangeAction;              /** Fixed range action */
    DirectoryPickerAction       _outputDirectoryAction;         /** Output directory picker action */
    StringAction                _fileNamePrefixAction;          /** File name prefix action */
    StatusAction                _statusAction;                  /** Status action */
    TriggersAction              _exportCancelAction;            /** Create and cancel triggers action */
    float                       _aspectRatio;                   /** Export image aspect ratio */
};