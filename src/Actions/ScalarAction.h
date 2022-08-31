#pragma once

#include "PluginAction.h"

#include "ScalarSourceAction.h"

using namespace hdps::gui;

/**
 * Scalar action class
 *
 * Action class picking scalar options (scalar by constant or dataset dimension)
 *
 * @author Thomas Kroes
 */
class ScalarAction : public PluginAction
{

    Q_OBJECT

protected: // Widget

    /** Widget class for scalar action */
    class Widget : public WidgetActionWidget {
    public:

        /**
         * Constructor
         * @param parent Pointer to parent widget
         * @param scalarAction Pointer to scalar action
         */
        Widget(QWidget* parent, ScalarAction* scalarAction);
    };

protected:

    /**
     * Get widget representation of the scalar action
     * @param parent Pointer to parent widget
     * @param widgetFlags Widget flags for the configuration of the widget (type)
     */
    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this);
    };

public:

    /**
     * Constructor
     * @param scatterplotPlugin Pointer to scatter plot plugin
     * @param minimum Scalar minimum value
     * @param maximum Scalar maximum value
     * @param value Scalar value
     * @param defaultValue Scalar default value
     */
    ScalarAction(ScatterplotPlugin* scatterplotPlugin, const QString& title, const float& minimum, const float& maximum, const float& value, const float& defaultValue);

    /**
     * Add dataset to the model
     * @param dataset Dataset to add
     */
    void addDataset(const Dataset<DatasetImpl>& dataset);

    /** Removes all datasets */
    void removeAllDatasets();

    /** Get smart pointer to current scalar dataset (if any) */
    Dataset<DatasetImpl> getCurrentDataset();

    /**
     * Set the current dataset
     * @param dataset Smart pointer to dataset
     */
    void setCurrentDataset(const Dataset<DatasetImpl>& dataset);

    /**
     * Set the current source index
     * @param sourceIndex Source index
     */
    void setCurrentSourceIndex(bool sourceIndex);

    /** Determines whether the scalar source is a constant */
    bool isSourceConstant() const;

    /** Determines whether the scalar source is a selection */
    bool isSourceSelection() const;

    /** Determines whether the scalar source is a dataset */
    bool isSourceDataset() const;

public: // Action getters

    DecimalAction& getMagnitudeAction() { return _magnitudeAction; }
    ScalarSourceAction& getSourceAction() { return _sourceAction; }

signals:

    /**
     * Signals that the source selection changed
     * @param sourceIndex Index of the selected source (0 is constant, 1 si selection, above is a dataset)
     */
    void sourceSelectionChanged(const std::uint32_t& sourceSelectionIndex);

    /**
     * Signals that the source data changed (only emitted when a source dataset is selected)
     * @param dataset Source dataset that changed
     */
    void sourceDataChanged(const Dataset<DatasetImpl>& dataset);

    /**
     * Signals that the scalar range changed
     * @param minimum Scalar range minimum
     * @param maximum Scalar range maximum
     */
    void scalarRangeChanged(const float& minimum, const float& maximum);

    /**
     * Signals that the scalar magnitude changed
     * @param magnitude Scalar magnitude
     */
    void magnitudeChanged(const float& magnitude);

    /**
     * Signals that the scalar offset changed
     * @param magnitude Scalar magnitude
     */
    void offsetChanged(const float& offset);

protected:
    DecimalAction           _magnitudeAction;   /** Scalar magnitude action */
    ScalarSourceAction      _sourceAction;      /** Scalar source action */
};