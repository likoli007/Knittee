#include "Knittee.h"


Knittee::Knittee(QWidget *parent)
    : QMainWindow(parent)
{
  
    ui.setupUi(this);

    auto topMenu = this->menuBar();//new QMenuBar(this);
    
    auto projectMenu = new QMenu(topMenu); //project menu: user can open a new project (2D or 3D), load an existing one etc.
    auto optionsMenu = new QMenu(topMenu); //options menu: will be able to select some basic options regarding the whole app
    auto helpAction = new QAction("&Help");    //help menu: get the current working context (2D or 3D modelling), and display related help
    auto aboutAction = new QAction("&About");   //about menu: display some basic information about the project
    



    //populate the projectMenu
    QAction* loadAction = new QAction("&Load", projectMenu);
    QAction* saveAction = new QAction("&Save", projectMenu);
    QAction* exitAction = new QAction("&Exit", projectMenu);


    QAction* newProjectAction = new QAction("&New", projectMenu);


    projectMenu->addAction(newProjectAction);
    projectMenu->addAction(loadAction);
    projectMenu->addAction(saveAction);
    projectMenu->addAction(exitAction);
    projectMenu->setTitle("Project");

    //populate the options menu
    //QAction* optionAction = new QAction("&Options", optionsMenu);
    //QAction* optionAction = optionsMenu->addAction("&options");
    //optionsMenu->addAction(optionAction);
    optionsMenu->setTitle("Options");

    //populate the help menu
    //helpMenu->setTitle("Help");

    //aboutAction->setTitle("About");

    topMenu->addMenu(projectMenu);
    topMenu->addMenu(optionsMenu);
    topMenu->addAction(helpAction);
    topMenu->addAction(aboutAction);

    setMenuBar(topMenu);

    QObject::connect(newProjectAction, &QAction::triggered, this, &Knittee::openNewProjectWindow);
    QObject::connect(loadAction, &QAction::triggered, this, &Knittee::selectProject);
    QObject::connect(exitAction, &QAction::triggered, this, &QApplication::quit);
    QObject::connect(optionsMenu, &QMenu::triggered, this, &Knittee::openOptionsWindow);
    QObject::connect(helpAction, &QAction::triggered, this, &Knittee::openHelpWindow);
    QObject::connect(aboutAction, &QAction::triggered, this, &Knittee::openAboutWindow);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    
    visualizerLayout = new QHBoxLayout(this);
    QVBoxLayout* optionsLayout = new QVBoxLayout(this);

    //QLabel* optionsLabel = new QLabel("Options about the visualization will be displayed here");
    //QWidget* spacer = new QWidget;
    //spacer->setMinimumSize(50, 0); 
    //toolsWidget = new QWidget(this);
    // Add the label and spacer to the layout
    //optionsLayout->addWidget(toolsWidget);
    //optionsLayout->addWidget(spacer);
    
    

    vis = new Visualizer(this); // Set the parent to be the main window
    vis->setMinimumSize(800, 600); // Set a minimum size for the widget
    visualizerLayout->addLayout(optionsLayout);
    visualizerLayout->addWidget(vis);

    mainLayout->addLayout(visualizerLayout);

    QWidget* centralWidget = new QWidget(this);
    centralWidget->setLayout(mainLayout);
    setCentralWidget(centralWidget);
}

void Knittee::selectProject() {
    QString filePath = QFileDialog::getOpenFileName(this, "Open Knittee Project File", "", "Custom Files (*.knittee)");
    loadProject(filePath);
}

void Knittee::loadProject(QString filePath) {
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		qDebug() << "Error opening project file!";
		return;
	}
    QTextStream in(&file);
	QString projectName = in.readLine();
	int projectType = in.readLine().toInt();
	file.close();
    QFileInfo projFileInfo(filePath); // Get the file info for the .proj file
    QString projFolder = projFileInfo.absolutePath(); // Get the folder containing the .proj file
    QString meshFilePath;
    
    if (projectType == 0)
    {
        meshFilePath = projFolder + "/mesh.obj"; // Construct the file path for mesh.obj
        ProjectInfo projectInfo;
        projectInfo.projectName = projectName;
        projectInfo.type = projectType;
        projectInfo.objectFilePath = meshFilePath;
        start3DProject(projectInfo);
    }

	
}

void Knittee::start3DProject(ProjectInfo context)
{
	qDebug() << "starting 3D project: " << context.objectFilePath;
    object_loader.setFilePath(context.objectFilePath);
    ObjectMesh mesh = object_loader.loadFile();
    vis->loadMesh(mesh);
    toolsWidget = new MeshToolBar(this);
    visualizerLayout->insertWidget(0, toolsWidget);
	// Handle selected options
	// context.type, context.file, etc. contain the selected values
}

void Knittee::openOptionsWindow()
{
    QMainWindow* optionsWindow = new QMainWindow(this);
    optionsWindow->setAttribute(Qt::WA_DeleteOnClose);
    optionsWindow->show();
    qDebug() << "Shown Options Menu...";
}

void Knittee::openAboutWindow()
{
    QDialog* aboutDialog = new QDialog(this); //needs to have 'this'?
    QString aboutText;

    aboutText = QString("Knittee is a 2D sheet and 3D object to Knit instruction CAD software \n") + 
        QString("designed by Alojz Holubek as part of his bachelor's thesis at Zhejiang University\n");
    

    aboutDialog->setWindowTitle("About This Software");
    QLabel* label = new QLabel(aboutText, aboutDialog);
    QVBoxLayout* layout = new QVBoxLayout(aboutDialog);
    layout->addWidget(label);
    aboutDialog->setLayout(layout);

    aboutDialog->exec();
    qDebug() << "Shown About Menu...";
}

void Knittee::openHelpWindow()
{
    QDialog* helpDialog = new QDialog(this); //needs to have 'this'?
    QString helpText;

    if (modellingType == 0) 
    {
        helpText = "Welcome to Knittee 2D/3D knitted objects modelling software!\n",
            "Instructions for 3D modelling: to be written here!";
    }
    else
    {
        helpText = "Welcome to Knittee 2D/3D knitted objects modelling software!\n",
            "Instructions for 2D modelling: to be written here!";
    }

    //show the help dialog with the relevant help information
    helpDialog->setWindowTitle("Knittee Help");


    QLabel* label = new QLabel(helpText, helpDialog);
    QVBoxLayout* layout = new QVBoxLayout(helpDialog);
    layout->addWidget(label);
    helpDialog->setLayout(layout);

    helpDialog->exec();
    qDebug() << "Shown Help Menu...";
}


void Knittee::openNewProjectWindow() 
{
    NewProjectDialog dialog;
    QObject::connect(&dialog, &NewProjectDialog::projectConfigurationsSelected, this, &Knittee::handleNewProject);
    dialog.exec();
}


void Knittee::setUpNew3DProject(ProjectInfo context)
{
    qDebug() << "Setting up 3D project with file: " << context.objectFilePath;
	//store the project information in the new folder
    QString projectFolderPath = QCoreApplication::applicationDirPath() + "/Projects/" + context.projectName;
    QDir().mkpath(projectFolderPath);
    QFile file(projectFolderPath+"/projectinfo.knittee");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
		qDebug() << "Error opening project file!";
		return;
	}
    file.write(context.projectName.toUtf8());
    file.write("\n");
    file.write(QString::number(context.type).toUtf8());
    file.write("\n");
    file.close();
    


    //copy the obj file to the project folder
    object_loader.setFilePath(context.objectFilePath);
    object_loader.copyObjFileToProject(context.projectName);

    context.objectFilePath = projectFolderPath + "/mesh.obj";
    start3DProject(context);
}

void Knittee::handleNewProject(ProjectInfo options)
{
    if (options.type == 0) {
        qDebug() << "Setting up 3D editing project...";
        setUpNew3DProject(options);
        //openFile(options.filePath);
    }
    else {
        qDebug() << "2D Sheet editing not yet implemented!";
    }

    // Handle selected options
    // options.type, options.file, etc. contain the selected values
}


Knittee::~Knittee()
{}



NewProjectDialog::NewProjectDialog(QWidget* parent) : QDialog(parent)
{
    setWindowTitle("Open New Project");
    setModal(true);         //dont allow the user to click on other things during project creation

    QVBoxLayout* mainLayout = new QVBoxLayout(this);


    QString projectNameLabelText = "Project Name:";
    QLabel* projectNameLabel = new QLabel(projectNameLabelText, this);
    projectNameLineEdit = new QLineEdit(this);
    
    QString projectTypeLabelText = "Select Your Project Type:";
    QLabel* projectTypeLabel = new QLabel(projectTypeLabelText, this);


    QHBoxLayout* projectTypeRadioLayout = new QHBoxLayout(this);
    meshRadio = new QRadioButton("3D Mesh editing", this);
    sheetRadio = new QRadioButton("2D Sheet editing", this);
    
    projectTypeRadioLayout->addWidget(meshRadio);
    projectTypeRadioLayout->addWidget(sheetRadio);

    QObject::connect(meshRadio, &QRadioButton::clicked, this, &NewProjectDialog::onMeshRadioClicked);
    QObject::connect(sheetRadio, &QRadioButton::clicked, this, &NewProjectDialog::onSheetRadioClicked);
    


    
    QHBoxLayout* projectConfirmDenyLayout = new QHBoxLayout(this);

    QPushButton* okButton = new QPushButton("OK", this);
    QPushButton* cancelButton = new QPushButton("Cancel", this);


    QObject::connect(okButton, &QPushButton::clicked, this, &NewProjectDialog::onConfirm);
    QObject::connect(cancelButton, &QPushButton::clicked, this, &NewProjectDialog::onCancel);


    projectConfirmDenyLayout->addWidget(okButton);
    projectConfirmDenyLayout->addWidget(cancelButton);


    /*Start of 3D project options
    * User must select an .obj mesh to be used in the project
    * TODO: add other options once the 3D portion grows more robust
    *       add a 'preform' option, package some basic .obj files with the program
    */
    meshSelectionLabel = new QLabel("Enter the path to the 3D mesh file (.obj):", this);
    meshSelectionLineEdit = new QLineEdit(this);
    meshSelectionButton = new QPushButton("Browse", this);
    meshSelectionLabel->setVisible(false);
    meshSelectionLineEdit->setVisible(false);
    
    meshSelectionLayout = new QHBoxLayout(this);
    meshSelectionLayout->addWidget(meshSelectionLineEdit);
    meshSelectionLayout->addWidget(meshSelectionButton);
    QObject::connect(meshSelectionButton, &QPushButton::clicked, this, &NewProjectDialog::openUserFile);



    /*Start of 2D project options
    * As the 2D portion gets developed, more options will be shown here, for now it is a simple sheet width/height set...
    * TODO: add other options once the 2D portion grows more robust
    */
    sheetSizeLabel = new QLabel("Enter the desired size of your sheet:");
    sheetSizeSeparator = new QLabel("/");
    widthLineEdit = new QLineEdit(this);
    heightLineEdit = new QLineEdit(this);
    
    

    QIntValidator* validator = new QIntValidator(this);
    widthLineEdit->setValidator(validator);
    heightLineEdit->setValidator(validator);

    sheetSizeLabel->setVisible(false);
    widthLineEdit->setVisible(false);
    heightLineEdit->setVisible(false);
    sheetSizeSeparator->setVisible(false);

    sheetSelectionLayout = new QHBoxLayout(this);
    sheetSelectionLayout->addWidget(widthLineEdit);
    sheetSelectionLayout->addWidget(sheetSizeSeparator);
    sheetSelectionLayout->addWidget(heightLineEdit);

    //mainLayout->addWidget(okButton);
    mainLayout->addWidget(projectNameLabel);
    mainLayout->addWidget(projectNameLineEdit);
    mainLayout->addWidget(projectTypeLabel);
    mainLayout->addLayout(projectTypeRadioLayout);
    mainLayout->addWidget(meshSelectionLabel);
    mainLayout->addLayout(meshSelectionLayout);
    mainLayout->addWidget(sheetSizeLabel);
    mainLayout->addLayout(sheetSelectionLayout);
    mainLayout->addLayout(projectConfirmDenyLayout);
    setLayout(mainLayout);
    meshRadio->setChecked(true);
    onMeshRadioClicked();
    
    qDebug() << "Shown Help Menu...";
}

void NewProjectDialog::openUserFile()
{
	QString file_path = QFileDialog::getOpenFileName(this, "Open OBJ File", "", "OBJ Files (*.obj)");
	meshSelectionLineEdit->setText(file_path);
}

void NewProjectDialog::onConfirm() 
{

    ProjectInfo projectInfo;
    projectInfo.projectName = projectNameLineEdit->text();

    if (meshRadio->isChecked()) {
        projectInfo.type = 0;
        projectInfo.objectFilePath = meshSelectionLineEdit->text();
    }
    else if (sheetRadio->isChecked()) {
        projectInfo.type = 1;
        qDebug() << "2D not yet implemented, selected size: " << widthLineEdit->text() << "/" << heightLineEdit->text();
    }
    close();
    emit projectConfigurationsSelected(projectInfo);
    
}

void NewProjectDialog::onCancel()
{
    close();
}
void NewProjectDialog::onMeshRadioClicked()
{
    meshSelectionLabel->setVisible(true);
    meshSelectionLineEdit->setVisible(true);
    meshSelectionButton->setVisible(true);
    sheetSizeLabel->setVisible(false);
    widthLineEdit->setVisible(false);
    heightLineEdit->setVisible(false);
    sheetSizeSeparator->setVisible(false);

}

void NewProjectDialog::onSheetRadioClicked()
{
    sheetSizeLabel->setVisible(true);
    widthLineEdit->setVisible(true);
    heightLineEdit->setVisible(true);
    sheetSizeSeparator->setVisible(true);
    meshSelectionLabel->setVisible(false);
    meshSelectionLineEdit->setVisible(false);
    meshSelectionButton->setVisible(false);
}


MeshToolBar::MeshToolBar(QWidget *parent) {
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QHBoxLayout* widthLayout = new QHBoxLayout(this);
    QHBoxLayout* heightLayout = new QHBoxLayout(this);

    QLabel* label = new QLabel("Stitch Width:", this);
    widthLayout->addWidget(label);

    QSlider* widthSlider = new QSlider(Qt::Horizontal);
    widthSlider->setMinimum(1);  // Set the minimum value
    widthSlider->setMaximum(10);  // Set the maximum value
    widthSlider->setValue(5);  // Set the initial value
    widthLayout->addWidget(widthSlider);

    QLabel* label2 = new QLabel("Stitch Height:", this);
    heightLayout->addWidget(label2);

    QSlider* heightSlider = new QSlider(Qt::Horizontal);
    heightSlider->setMinimum(1);  // Set the minimum value
    heightSlider->setMaximum(10);  // Set the maximum value
    heightSlider->setValue(5);  // Set the initial value
    heightLayout->addWidget(heightSlider);

    QLabel* label3 = new QLabel("Further Options...", this);

    mainLayout->addLayout(widthLayout);
    mainLayout->addLayout(heightLayout);
    mainLayout->addWidget(label3);
    setLayout(mainLayout);
    setMinimumSize(300,100);
}