#pragma once
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
#include <QTextEdit>
#include <QComboBox>
#include <qbuttongroup.h>

#include "ui_Knittee.h"
#include "ObjHandler.h"
#include "Visualizer.h"
#include "MeshToolBar.h"
#include "KnitGrapher.h"
#include "Stitch.h"
#include "Link.h"
#include "LaceKnitter.h"
#include "FlatPoint.h"
#include "KnitOutStitch.h"
#include "TracedStitch.h"
#include "KnitoutSheduler.h"
#include "SheetToolBar.h"

/*
* This is the main class that will be used to create the GUI for the program
* It also includes helper classes that will be shown at some stages of the program
* such as options related only to one type of project (3D or 2D)
* or a dialog that will be shown when the user wants to create a new project
*/


/*
* This struct is used to define the project types, it is saved into .knittee files and can also be loaded into the program from them
*/
struct ProjectInfo {
    int type;    // same as modellingType, will pass the user selected option
    QString objectFilePath; //path of the 3D object/ 2D sheet that is used by the project
    QString projectName; //name of the project
    int width; //width of the sheet
    int height; //height of the sheet
    int racking; //racking value of the sheet
};


/*
* This class is used to create the new project popup window when the user selects 'new' in the project menu
*/
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
    void onMeshSelectionRadioClicked();
    void onPreformRadioClicked();

private:
    ProjectInfo projectInfo;
    QRadioButton* sheetRadio;
    QRadioButton* meshRadio;

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

    QRadioButton* meshSelectionRadio;
    QRadioButton* preformSelectionRadio;
    QComboBox* preformSelectionComboBox;
    QLabel* preformSelectionLabel;
    QLabel* meshTypeLabel;
    //QHBoxLayout* meshTypeLayout;
    //QSeparator preformSelectionSeparator;
};


/*
* Knittee class is the main 'application' class of the whole project, it includes an instance of most other classes of this project
*   it is also in charge of helping the different components interface with each other
*/
class Knittee : public QMainWindow
{
    Q_OBJECT

public:
    Knittee(QWidget* parent = nullptr);
    ~Knittee();

private slots:
    void setConstraintsMode();
    void handleToolbarDone();
    void startRemeshing();
    void meshInterpolated(ObjectMesh, std::vector<float>);
    void firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
        std::vector< std::vector< Stitch > >* active_stitches,
        RowColGraph* graph);
    void peelSliceDone(ObjectMesh*, std::vector< std::vector< uint32_t > > *, std::vector< std::vector< uint32_t > >*);
    void linkChainsDone(std::vector< std::vector< Stitch > >* , std::vector< Link >* );
    void nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >*);
    void knitGraphCreated();
    void knitGraphTraced(std::vector< TracedStitch >*);
    void helpBoxCommunication(QString);
    void saveConstraints();
    void resetButtonClicked();
    void instructionsCreated(std::vector<std::string> instructions);
    void generateKnitoutButtonClicked();
    void knitoutGenerated(std::vector<QString>);
private:
    void saveTraced(std::vector< TracedStitch >*);


    int modellingType = -1; //is the user operating on a 3D model (0) or a 2D sheet? (1), perhaps could be an enum?
    void closeEvent(QCloseEvent* event) override;
    Ui::KnitteeClass ui;
    ObjHandler object_loader;
    Visualizer* vis;
    KnitGrapher knitGrapher;
    KnitoutScheduler knitoutScheduler;
    LaceKnitter laceKnitter;
    MeshToolBar* meshToolsWidget;
    SheetToolBar* sheetToolsWidget;
    QHBoxLayout* visualizerLayout;
    QString projectPath;
    //bottom textedit for messages to user
    QTextEdit* messageTextEdit;

    ProjectInfo context;

    bool canSlice = false;

    void loadProject(QString);
    void openOptionsWindow();
    void openHelpWindow();
    void openAboutWindow();
    void openNewProjectWindow();
    void handleNewProject(ProjectInfo info);
    void selectProject();

    void saveProject();

    //3D project functions
    void setUpNew3DProject();
    void start3DProject();
    void save3DProject();
    void loadConstraints();

    //2D project functions
    void setUpNew2DProject();
    void start2DProject();
    void save2DProject();
    void save2DContext();

protected:
    void KeyPressEvent(QKeyEvent* event);
};
