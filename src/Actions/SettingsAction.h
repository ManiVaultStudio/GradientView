#pragma once

#include "PluginAction.h"

#include "RenderModeAction.h"
#include "PlotAction.h"
#include "PositionAction.h"
#include "SubsetAction.h"
#include "MiscellaneousAction.h"
#include "FilterAction.h"
#include "OverlayAction.h"

#include "actions/WidgetActionStateWidget.h"

using namespace hdps::gui;

class ScatterplotPlugin;

class SettingsAction : public PluginAction
{
public:
    class SpacerWidget : public QWidget {
    public:
        enum class Type {
            Divider,
            Spacer
        };

    public:
        SpacerWidget(const Type& type = Type::Divider);

        static Type getType(const WidgetActionWidget::State& widgetTypeLeft, const WidgetActionWidget::State& widgetTypeRight);
        static Type getType(const hdps::gui::WidgetActionStateWidget* stateWidgetLeft, const hdps::gui::WidgetActionStateWidget* stateWidgetRight);

        void setType(const Type& type);
        static std::int32_t getWidth(const Type& type);

    protected:
        Type            _type;
        QHBoxLayout*    _layout;
        QFrame*         _verticalLine;
    };

protected: // Widget

    class Widget : public WidgetActionWidget {
    public:
        Widget(QWidget* parent, SettingsAction* settingsAction);

        bool eventFilter(QObject* object, QEvent* event);

    protected:
        void addStateWidget(WidgetAction* widgetAction, const std::int32_t& priority = 0);

    private:
        void updateLayout();

    protected:
        QHBoxLayout                         _layout;
        QWidget                             _toolBarWidget;
        QHBoxLayout                         _toolBarLayout;
        QVector<WidgetActionStateWidget*>   _stateWidgets;
        QVector<SpacerWidget*>              _spacerWidgets;

        friend class SettingsAction;
    };

    QWidget* getWidget(QWidget* parent, const std::int32_t& widgetFlags) override {
        return new Widget(parent, this);
    };

public:
    SettingsAction(ScatterplotPlugin* scatterplotPlugin);

    QMenu* getContextMenu();

    RenderModeAction& getRenderModeAction() { return _renderModeAction; }
    PositionAction& getPositionAction() { return _positionAction; }
    SubsetAction& getSubsetAction() { return _subsetAction; }
    PlotAction& getPlotAction() { return _plotAction; }
    TriggerAction& getExportAction() { return _exportAction; }
    MiscellaneousAction& getMiscellaneousAction() { return _miscellaneousAction; }
    FilterAction& getFilterAction() { return _filterAction; }

protected:
    RenderModeAction            _renderModeAction;
    PositionAction              _positionAction;
    SubsetAction                _subsetAction;
    PlotAction                  _plotAction;
    TriggerAction               _exportAction;
    MiscellaneousAction         _miscellaneousAction;
    FilterAction                _filterAction;
    OverlayAction               _overlayAction;
};
