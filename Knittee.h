#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_Knittee.h"
#include "ObjLoader.h"
#include "Visualizer.h"

class Knittee : public QMainWindow
{
    Q_OBJECT

public:
    Knittee(QWidget *parent = nullptr);
    ~Knittee();
    void openFile();

private:
    Ui::KnitteeClass ui;
    ObjLoader object_loader;
    Visualizer* vis;
};
