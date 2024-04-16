#include "Knittee.h"
#include "Visualizer.h"

#include <qmenubar.h>
#include <qmenu.h>
#include <qboxlayout.h>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QFileDialog>

Knittee::Knittee(QWidget *parent)
    : QMainWindow(parent)
{
  
    ui.setupUi(this);

    auto topMenu = this->menuBar();//new QMenuBar(this);
    
    auto fileMenu = new QMenu(topMenu);

    QAction* loadAction = new QAction("&Load", fileMenu);
    QAction* saveAction = new QAction("&Save", fileMenu);
    QAction* exitAction = new QAction("E&xit", fileMenu);
    fileMenu->addAction(loadAction);
    fileMenu->addAction(saveAction);
    fileMenu->addAction(exitAction);
    fileMenu->setTitle("File");

    topMenu->addMenu(fileMenu);

    setMenuBar(topMenu);

    QObject::connect(exitAction, &QAction::triggered, this, &QApplication::quit);
    QObject::connect(loadAction, &QAction::triggered, this, &Knittee::openFile);
   

    auto mainLayout = new QVBoxLayout(this);
    

    vis = new Visualizer(this); // Set the parent to be the main window
    //vis->setMinimumSize(400, 400); // Set a minimum size for the widget
    mainLayout->addWidget(vis);

    QWidget* centralWidget = new QWidget(this);
    centralWidget->setLayout(mainLayout);
    setCentralWidget(centralWidget);
}

void Knittee::openFile()
{
    QString file_path = QFileDialog::getOpenFileName(this, "Open OBJ File", "", "OBJ Files (*.obj)");

    if (!file_path.isEmpty()) {       
        ObjectMesh mesh = object_loader.loadFile(file_path);
        vis->loadMesh(mesh);
    }
}

Knittee::~Knittee()
{}
