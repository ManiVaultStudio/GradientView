#pragma once

#include <QWidget>
#include <QListWidget>

#include <ClusterData/ClusterData.h>

class MetadataView : public QWidget
{
public:
    MetadataView();

    void addListItem(mv::Dataset<Clusters> dataset);
private:
    QListWidget* _listWidget;
};
