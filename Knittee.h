#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_Knittee.h"

class Knittee : public QMainWindow
{
    Q_OBJECT

public:
    Knittee(QWidget *parent = nullptr);
    ~Knittee();

private:
    Ui::KnitteeClass ui;
};
