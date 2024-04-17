#pragma once


#include "ui_Knittee.h"
#include "ObjHandler.h"
#include "Visualizer.h"

#include <QtWidgets/QMainWindow>
#include <QDialog>
#include <QWidget>
#include <QRadioButton>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QHBoxLayout>
#include <QSlider>
#include <qmenubar.h>
#include <qmenu.h>
#include <qboxlayout.h>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QFileDialog>
#include <QPushButton>


/*
* This is the main class that will be used to create the GUI for the program
* It also includes helper classes that will be shown at some stages of the program
* such as options related only to one type of project (3D or 2D)
* or a dialog that will be shown when the user wants to create a new project
*/


// This class will be called to be shown when the user selects a 3D mesh project
class MeshToolBar : public QWidget {
	Q_OBJECT

    public:
        MeshToolBar(QWidget* parent = nullptr);
private:
    // Will have a slew of functions that pass the user specified options to the visualizer/algorithm program

};


struct ProjectInfo {
    int type;    // same as modellingType, will pass the user selected option
    QString objectFilePath; //path of the 3D object/ 2D sheet that is used by the project
    QString projectName; //name of the project
    //TODO: more options as the program expands
};

class NewProjectDialog : public QDialog
{
    Q_OBJECT

public:
    NewProjectDialog(QWidget* parent = nullptr);

signals:
    void projectConfigurationsSelected(ProjectInfo options);
private slots:
    void onConfirm();
    void onCancel();
    void onMeshRadioClicked();
    void onSheetRadioClicked();
    void openUserFile();

private:
    ProjectInfo projectInfo;
    QRadioButton* sheetRadio;
    QRadioButton* meshRadio;    //cannot write 3D :C


    //3D project option components
    QLabel* meshSelectionLabel;
    QLineEdit* meshSelectionLineEdit;
    QPushButton* meshSelectionButton;
    QHBoxLayout* meshSelectionLayout;
    QLineEdit* projectNameLineEdit;

    //2D project option components
    QHBoxLayout* sheetSelectionLayout;
    QLabel* sheetSizeLabel;
    QLabel* sheetSizeSeparator;
    QLineEdit* widthLineEdit;
    QLineEdit* heightLineEdit;
};



class Knittee : public QMainWindow
{
    Q_OBJECT




public:
    Knittee(QWidget *parent = nullptr);
    ~Knittee();
    //void openFile(QString filePath);

private:
    int modellingType = 0; //is the user operating on a 3D model (0) or a 2D sheet? (1), perhaps could be an enum?


    Ui::KnitteeClass ui;
    ObjHandler object_loader;
    Visualizer* vis;
    QWidget* toolsWidget;
    //QVBoxLayout* mainLayout;
    QHBoxLayout* visualizerLayout;

    void loadProject(QString filePath);
    void openOptionsWindow();
    void openHelpWindow();
    void openAboutWindow();
    void openNewProjectWindow();
    void handleNewProject(ProjectInfo options);
    void openFile(QString filePath);
    void setUpNew3DProject(ProjectInfo options);
    void start3DProject(ProjectInfo context);
    void selectProject();
};
