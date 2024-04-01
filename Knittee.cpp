#include "Knittee.h"

Knittee::Knittee(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);
    
    QWidget* topMenuWidget = new QWidget;
    setCentralWiget(topMenuWidget);


}

QMenuBar topMenu(QWidget* parent) 
{


}

Knittee::~Knittee()
{}
