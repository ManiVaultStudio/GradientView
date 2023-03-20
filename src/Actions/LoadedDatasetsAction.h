#pragma once

#include "PluginAction.h"

#include "actions/DatasetPickerAction.h"

using namespace hdps::gui;

class LoadedDatasetsAction : public PluginAction
{
protected:

    class Widget : public WidgetActionWidget {
    public:
        Widget(QWidget* parent, LoadedDatasetsAction* currentDatasetAction, const std::int32_t& widgetFlags);
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this, widgetFlags);
    };

public:
    LoadedDatasetsAction(ScatterplotPlugin* scatterplotPlugin);

public: // Serialization

    /**
     * Load widget action from variant map
     * @param Variant map representation of the widget action
     */
    void fromVariantMap(const QVariantMap& variantMap) override;

    /**
     * Save widget action to variant map
     * @return Variant map representation of the widget action
     */
    QVariantMap toVariantMap() const override;

protected:
    DatasetPickerAction	    _positionDatasetPickerAction;
    DatasetPickerAction     _colorDatasetPickerAction;

    friend class Widget;
};
