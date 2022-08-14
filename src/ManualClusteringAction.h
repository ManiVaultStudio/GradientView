#pragma once

#include "PluginAction.h"

#include "actions/DatasetPickerAction.h"

using namespace hdps;
using namespace hdps::gui;
using namespace hdps::util;

class Clusters;

/**
 * Manual clustering action class
 *
 * Action class for manual clustering
 *
 * @author Thomas Kroes
 */
class ManualClusteringAction : public PluginAction
{
protected:

    class Widget : public WidgetActionWidget
    {
    public:

        /**
         * Constructor
         * @param parent Pointer to parent widget
         * @param manualClusteringAction Pointer to manual clustering action
         * @param widgetFlags Widget flags for the configuration of the widget (type)
         */
        Widget(QWidget* parent, ManualClusteringAction* manualClusteringAction, const std::int32_t& widgetFlags);
    };

    /**
     * Get widget representation of the cluster action
     * @param parent Pointer to parent widget
     * @param widgetFlags Widget flags for the configuration of the widget (type)
     */
    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this, widgetFlags);
    };

public:

    /**
     * Constructor
     * @param scatterplotPlugin Pointer to scatter plot plugin
     */
    ManualClusteringAction(ScatterplotPlugin* scatterplotPlugin);

    /** Adds a clusters dataset to the position dataset as a child */
    void createDefaultClusterDataset();

    /** Update the target cluster datasets action (creates default set if no cluster sets are available) */
    void updateTargetClusterDatasets();

protected:

    /** Update the state of the actions */
    void updateActions();

public: // Action getters

    StringAction& getNameAction() { return _nameAction; }
    ColorAction& getColorAction() { return _colorAction; }
    TriggerAction& getCreateCluster() { return _addClusterAction; }
    TriggerAction& getAddClusterAction() { return _addClusterAction; }
    DatasetPickerAction& getTargetClusterDataset() { return _targetClusterDataset; }

protected:
    StringAction            _nameAction;                /** Cluster name action */
    ColorAction             _colorAction;               /** Cluster color action */
    TriggerAction           _addClusterAction;          /** Add manual cluster action */
    DatasetPickerAction     _targetClusterDataset;      /** Target cluster dataset action */
};
