#include "Knittee.h"


Knittee::Knittee(QWidget* parent)
    : QMainWindow(parent)
{

    ui.setupUi(this);
    setFocusPolicy(Qt::StrongFocus);

    //set up the menu bar
    auto topMenu = this->menuBar();        
    auto projectMenu = new QMenu(topMenu); //project menu: user can open a new project (2D or 3D), load an existing one etc.
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

    topMenu->addMenu(projectMenu);
    topMenu->addAction(helpAction);
    topMenu->addAction(aboutAction);
    
    setMenuBar(topMenu);


    //multithread the mesh->knitgraph process, as otherwise the GUI will freeze during operations, making the mesh uninteractable
    knitGrapherThread = new QThread();
    knitGrapher.moveToThread(knitGrapherThread);
    knitGrapherThread->start();


    //set up the toolbars (at the start both are not visible but have them prepared)
    meshToolsWidget = new MeshToolBar(this);
    sheetToolsWidget = new SheetToolBar(this);

    //set up the visualizer widget, which is comprised of a left toolbar and the visualizer itself
    visualizerLayout = new QHBoxLayout(this);
    vis = new Visualizer(this);         
    vis->setMinimumSize(800, 600);       
    visualizerLayout->addWidget(meshToolsWidget);
    visualizerLayout->addWidget(sheetToolsWidget);
    meshToolsWidget->hide();
    sheetToolsWidget->hide();
    visualizerLayout->addWidget(vis);

   
    //set up the layout used for the main window
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    

    

    //connect the different signals and slots of the member objects between themselves, allowing for communication between them

    //menu actions
    QObject::connect(newProjectAction, &QAction::triggered, this, &Knittee::openNewProjectWindow);
    QObject::connect(loadAction, &QAction::triggered, this, &Knittee::selectProject);
    QObject::connect(exitAction, &QAction::triggered, this, &QApplication::quit);
    QObject::connect(saveAction, &QAction::triggered, this, &Knittee::saveProject);
    QObject::connect(helpAction, &QAction::triggered, this, &Knittee::openHelpWindow);
    QObject::connect(aboutAction, &QAction::triggered, this, &Knittee::openAboutWindow);

    //3D toolbar actions
    QObject::connect(meshToolsWidget, SIGNAL(constraintsButtonClicked(bool)), vis, SLOT(setConstraintsMode(bool)));
    QObject::connect(meshToolsWidget, SIGNAL(remeshButtonClicked()), this, SLOT(startRemeshing()));
    QObject::connect(meshToolsWidget, SIGNAL(widthChanged(float)), &knitGrapher, SLOT(setStitchWidth(float)));
    QObject::connect(meshToolsWidget, SIGNAL(heightChanged(float)), &knitGrapher, SLOT(setStitchHeight(float)));
    QObject::connect(meshToolsWidget, SIGNAL(unitChanged(float)), &knitGrapher, SLOT(setModelUnitLength(float)));
    QObject::connect(meshToolsWidget, SIGNAL(stepButtonClicked(int)), &knitGrapher, SLOT(stepButtonClicked(int)));
    QObject::connect(meshToolsWidget, SIGNAL(traceButtonClicked()), &knitGrapher, SLOT(traceButtonClicked()));
    QObject::connect(meshToolsWidget, SIGNAL(showInterpolatedChanged(int)), vis, SLOT(showInterpolatedChanged(int)));
    QObject::connect(meshToolsWidget, SIGNAL(showGraphChanged(int)), vis, SLOT(showGraphChanged(int)));
    QObject::connect(meshToolsWidget, SIGNAL(showTracedChanged(int)), vis, SLOT(showTracedChanged(int)));
    QObject::connect(meshToolsWidget, SIGNAL(showYarnChanged(int)), vis, SLOT(showYarnChanged(int)));
    QObject::connect(meshToolsWidget, SIGNAL(helpBoxCommunication(QString)), this, SLOT(helpBoxCommunication(QString)));
    QObject::connect(meshToolsWidget, SIGNAL(resetButtonClicked()), this, SLOT(resetButtonClicked()));
    QObject::connect(meshToolsWidget, SIGNAL(generateKnitoutButtonClicked()), this, SLOT(generateKnitoutButtonClicked()));

    //2D toolbar actions
    QObject::connect(sheetToolsWidget, SIGNAL(rackingChanged(int)), &laceKnitter, SLOT(rackingChanged(int)));
    QObject::connect(sheetToolsWidget, SIGNAL(widthChanged(int, int)), &laceKnitter, SLOT(widthChanged(int, int)));
    QObject::connect(sheetToolsWidget, SIGNAL(heightChanged(int, int)), &laceKnitter, SLOT(heightChanged(int, int)));
    QObject::connect(sheetToolsWidget, SIGNAL(generateKnitoutSheet(int)), this, SLOT(generateKnitoutSheet(int)));

    //visualizer actions
    QObject::connect(vis, SIGNAL(requestConstraintsSave()), this, SLOT(saveConstraints()));
    QObject::connect(vis, SIGNAL(moveLoop(QPair<int, int>, QPair<int, int>)), &laceKnitter, SLOT(moveLoop(QPair<int, int>, QPair<int, int>)));
    
    //3D algorithmic actions
    QObject::connect(&knitGrapher, SIGNAL(knitGraphInterpolated(ObjectMesh, std::vector<float>)), this, SLOT(meshInterpolated(ObjectMesh, std::vector<float>)));
    QObject::connect(&knitGrapher, SIGNAL(firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >*, std::vector< std::vector< Stitch > >*, RowColGraph*)), 
        this, SLOT(firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >*, std::vector< std::vector< Stitch > >*, RowColGraph*)));
    QObject::connect(&knitGrapher, SIGNAL(peelSliceDone(ObjectMesh, std::vector< std::vector< uint32_t > >, std::vector< std::vector< uint32_t > >)), vis, SLOT(peelSliceDone(ObjectMesh, std::vector< std::vector< uint32_t > > , std::vector< std::vector< uint32_t > >)));
    QObject::connect(&knitGrapher, SIGNAL(linkChainsDone(std::vector< std::vector< Stitch > >*, std::vector< Link >*)), this, SLOT(linkChainsDone(std::vector< std::vector< Stitch > >*, std::vector< Link >*)));
    QObject::connect(&knitGrapher, SIGNAL(nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >*)), this, SLOT(nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >*)));
    QObject::connect(&knitGrapher, SIGNAL(knitGraphCreated()), this, SLOT(knitGraphCreated()));
    QObject::connect(&knitGrapher, SIGNAL(knitGraphTraced(std::vector< TracedStitch >*)), this, SLOT(knitGraphTraced(std::vector< TracedStitch >*)));
    QObject::connect(&knitGrapher, SIGNAL(knitGraphTraced(std::vector< TracedStitch >*)), meshToolsWidget, SLOT(knitGraphTraced()));
    QObject::connect(&knitGrapher, SIGNAL(knitGraphTraced(std::vector< TracedStitch >*)), vis, SLOT(knitGraphTraced(std::vector< TracedStitch >*)));
    QObject::connect(&knitGrapher, SIGNAL(unlockMeshButtons()), meshToolsWidget, SLOT(unlockMeshButtons()));
    QObject::connect(&knitGrapher, SIGNAL(helpBoxCommunication(QString)), this, SLOT(helpBoxCommunication(QString)));
    QObject::connect(&knitoutScheduler, SIGNAL(helpBoxCommunication(QString)), this, SLOT(helpBoxAppend(QString)));
    QObject::connect(&knitoutScheduler, SIGNAL(instructionsCreated(std::vector<std::string>)), this, SLOT(instructionsCreated(std::vector<std::string>)));
    QObject::connect(&knitoutScheduler, SIGNAL(knitoutGenerated(std::vector<QString>)), this, SLOT(knitoutGenerated(std::vector<QString>)));

    //2D algorithmic actions
    QObject::connect(&laceKnitter, SIGNAL(sheetChanged(std::vector<std::vector<FlatPoint>>*)), vis, SLOT(sheetChanged(std::vector<std::vector<FlatPoint>>*)));
    QObject::connect(&laceKnitter, SIGNAL(sheetDimensionsChanged(int, int, int)), this, SLOT(sheetDimensionsChanged(int, int, int)));
    QObject::connect(&laceKnitter, SIGNAL(knitoutGenerated(std::vector<QString>)), this, SLOT(knitoutGenerated(std::vector<QString>)));
    QObject::connect(&laceKnitter, SIGNAL(helpBoxCommunication(QString)), this, SLOT(helpBoxCommunication(QString)));

    //set up the message text box at the bottom of the window
    messageTextEdit = new QTextEdit(this);
    messageTextEdit->setReadOnly(true);
    messageTextEdit->setMinimumSize(800, 100);
    messageTextEdit->setText(HelperText::welcomeText);

    //set up the main layout of the window
    mainLayout->addLayout(visualizerLayout);
    mainLayout->addWidget(messageTextEdit);
    mainLayout->setStretchFactor(visualizerLayout, 4);
    mainLayout->setStretchFactor(messageTextEdit, 1);
    QWidget* centralWidget = new QWidget(this);
    centralWidget->setLayout(mainLayout);
    setCentralWidget(centralWidget);
}



/*
*   Function: start the 2D knitout generation process
*   Arguments: int algorithm - the algorithm to use for knitout generation
*/
void Knittee::generateKnitoutSheet(int algorithm) 
{
    
    laceKnitter.generateKnitout(algorithm);
}

/*
*  Function: save the generated knitout instructions into a file inside the project folder
*  Arguments: std::vector<QString> knitout - the knitout instructions to save
*/
void Knittee::knitoutGenerated(std::vector<QString> knitout) 
{


    //save to project folder under 'knitout.k'
    QString filePath = projectPath + "/knitout.k";
    QFile file(filePath);

    if (!file.open(QIODevice::WriteOnly))
    {
        qDebug() << "could not open file";
        return;
    }

    QTextStream out(&file);
    for (const QString& instruction : knitout)
    {
        out << instruction << '\n';
    }
    file.close();
}

/*
*   Function: start the 3D knitout generation process, is a slot function
*/
void Knittee::generateKnitoutButtonClicked() 
{
    helpBoxCommunication(HelperText::meshKnitoutStartText);


    //multithread the operation, as otherwise the whole GUI freezes while the (long) process is running
    knitoutGeneratorThread = new QThread();
    knitoutScheduler.setFilePath(projectPath + "/traced");
    knitoutScheduler.moveToThread(knitoutGeneratorThread);
    connect(knitoutGeneratorThread, &QThread::started, &knitoutScheduler, &KnitoutScheduler::schedule);
    knitoutGeneratorThread->start();
}

/*
*  Functions: show or append a message to the help box at the bottom of the window, are slot functions
*  Arguments: QString message - the message to append
*/
void Knittee::helpBoxAppend(QString message) 
{
    messageTextEdit->append(message);
}
void Knittee::helpBoxCommunication(QString message) 
{
    messageTextEdit->setText(message);
}

/*
*   Function: save the midpoint between the Traced graph and knitout, is a slot function
*   Arguments: std::vector<std::string> instructions - the scheduling instructions to save
*/
void Knittee::instructionsCreated(std::vector<std::string> instructions) 
{
    //save to project folder under 'kgenerator'
    QString filePath = projectPath + "/kgenerator";
    QFile file(filePath);

    if (!file.open(QIODevice::WriteOnly))
	{
		qDebug() << "could not open file";
		return;
	}

	QTextStream out(&file);
    for (const std::string& instruction : instructions)
	{
		out << QString::fromStdString(instruction) << '\n';
	}
	file.close();
}

/*
*   Function: make the user select a project file to open
*/
void Knittee::selectProject()
{
    QString filePath = QFileDialog::getOpenFileName(this, "Open Knittee Project File", "", "Custom Files (*.knittee)");
    QFileInfo fileInfo(filePath);
    projectPath = fileInfo.absolutePath();
    loadProject(filePath);
}

/*
*   Function: load a project file and set up the project based on the information in the file
*   Arguments: QString filePath - the path to the project file
*/
void Knittee::loadProject(QString filePath) 
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "Error opening project file!";
        return;
    }

    QTextStream in(&file);
    QString projectName = in.readLine();
    int projectType = in.readLine().toInt();
    QFileInfo projFileInfo(filePath); // Get the file info for the .proj file
    QString projFolder = projFileInfo.absolutePath(); // Get the folder containing the .proj file
    QString meshFilePath;

    if (projectType == 0)
    {
        meshFilePath = projFolder + "/mesh.obj"; // Construct the file path for mesh.obj
        context.projectName = projectName;
        context.type = projectType;
        context.objectFilePath = meshFilePath;
        file.close();
        start3DProject();
    }
    else
    {
        int width = in.readLine().toInt();
        int height = in.readLine().toInt();
        int rack = in.readLine().toInt();

        context.projectName = projectName;
        context.type = projectType;
        context.width = width;
        context.height = height;
        context.racking = rack;

        file.close();
        laceKnitter.setDimensions(width, height, rack);
        laceKnitter.loadFromFile(projFolder);
        start2DProject();
    }
}

/*
*   Function: load previously created constraints for the current 3D project
*/
void Knittee::loadConstraints() 
{
    std::vector<Constraint*> constraints;
    QString filePath = projectPath + "/constraints";
    QFile file(filePath);
    if (file.exists()) {
        if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QTextStream in(&file);
            while (!in.atEnd()) {
                Constraint* constraintPointer = new Constraint();
                QString line = in.readLine();
                constraintPointer->timeValue = line.toFloat();
                line = in.readLine();
                QStringList indices = line.split(';');
                if (!indices.isEmpty())
                {
                    for (const QString& index : indices) {
                        constraintPointer->vertices.push_back(index.toInt());
                    }
                    constraints.push_back(constraintPointer);
                }
            }
            file.close();
        }
        else {
            qDebug() << "Constraints file open error!";
        }
    }
    else {
        qDebug() << "Constraints file does not exist error!";
    }
    vis->setConstraints(constraints);
}

/*
*  Function: start up the relevant objects for a 3D project
*/
void Knittee::start3DProject()
{
    object_loader.setFilePath(context.objectFilePath);
    meshToolsWidget->resetButtonClicked();

    ObjectMesh mesh = object_loader.loadFile();

    loadConstraints();

    modellingType = 0;
    vis->projectType = 0;
    vis->loadMesh(mesh);
    knitGrapher.setOriginalMesh(mesh);

    helpBoxCommunication(HelperText::project3DText);

    sheetToolsWidget->hide();
    meshToolsWidget->show();
    visualizerLayout->setStretchFactor(vis, 3);
    visualizerLayout->setStretchFactor(sheetToolsWidget, 0);
    visualizerLayout->setStretchFactor(meshToolsWidget, 1);
}

/*
*  Function: start up the relevant objects for a 2D project
*/
void Knittee::start2DProject()
{
    vis->projectType = 1;
    modellingType = 1;

    laceKnitter.setDimensions(context.width, context.height, context.racking);
    laceKnitter.loadFromFile(projectPath);

    sheetToolsWidget->show();
    meshToolsWidget->hide();
    visualizerLayout->setStretchFactor(vis, 3);
    visualizerLayout->setStretchFactor(sheetToolsWidget, 1);
    visualizerLayout->setStretchFactor(meshToolsWidget, 0);
}

/*
*  Function: a cleaner function that saves the project when the user closes the application
*/
void Knittee::closeEvent(QCloseEvent* event) 
{   
    saveProject();
}


/*
*   Function: the following 2 functions save the sheet changes made by the user when the user either clicks 'save' or exits the program
*       Technically the user may not want to save changes when exiting...
*/
void Knittee::saveProject()
{
    //3D projects dont need to be saved as all their changes are saved in the mesh files, so only save 2D projects
    if (modellingType == 1){
		save2DProject();
	}
}
void Knittee::save2DProject() 
{
    QString projectFolderPath = QCoreApplication::applicationDirPath() + "/Projects/" + context.projectName;

    laceKnitter.saveToFile(projectFolderPath);
    
    context.height = laceKnitter.getHeight();
    context.width = laceKnitter.getWidth();
    context.racking = laceKnitter.getRacking();

    QFile file(projectFolderPath + "/" + context.projectName + ".knittee");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "Error opening project file!";
        return;
    }
    file.write(context.projectName.toUtf8());
    file.write("\n");
    file.write(QString::number(context.type).toUtf8());
    file.write("\n");
    file.write(QString::number(context.width).toUtf8());
    file.write("\n");
    file.write(QString::number(context.height).toUtf8());
    file.write("\n");
    file.write(QString::number(context.racking).toUtf8());
    file.close();
}

/*
*   Function: reset the 3D project such that the user can start over
*/
void Knittee::resetButtonClicked() 
{
	vis->reset();
	knitGrapher.reset();
}

/*
*   Function: handle the creation of a traced graph - save it to project, is a slot function
*   Arguments: std::vector< TracedStitch >* traced_stitches - the traced stitches to save
*/
void Knittee::knitGraphTraced(std::vector< TracedStitch >* traced_stitches) 
{
    saveTraced(traced_stitches);
}
void Knittee::saveTraced(std::vector< TracedStitch >* traced_stitches)  //actually save the traced stitches to a file
{
	QString filePath = projectPath + "/traced";
	QFile file(filePath);

    if (!file.open(QIODevice::WriteOnly))
    {
        qDebug() << "could not open file";
        return;
    }
   
    QTextStream out(&file);

    for (auto const& ts : *traced_stitches)
    {
        out << ts.yarn
            << ' ' << ts.type
            << ' ' << ts.dir
            << ' ' << (int32_t)ts.ins[0]
            << ' ' << (int32_t)ts.ins[1]
            << ' ' << (int32_t)ts.outs[0]
            << ' ' << (int32_t)ts.outs[1]
            << ' ' << ts.at.x << ' ' << ts.at.y << ' ' << ts.at.z << '\n';
    }

    file.close();
}

/*
*  Function: handle the creation of a knit graph - warn other member objects, is a slot function
*/
void Knittee::knitGraphCreated()
{
    meshToolsWidget->knitGraphCreated();
    vis->knitGraphCreated();
}

/*
*   Functions: handle the inputs outputs of different steps of the knitGrapher algorithm, are slot functions
*       four outputs: for now only the visualizer is warned, but later on the user may be warned as well, 
            maybe save a history of changes to be replayed?
*   Arguments: relevant data structures for the different steps
*/
void Knittee::nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >* active_chains) 
{
	vis->nextActiveChainsDone(active_chains);
}
void Knittee::linkChainsDone(std::vector< std::vector< Stitch > >* next_stitches, std::vector< Link >* links)
{
	vis->linkChainsDone(next_stitches, links);
}
void Knittee::peelSliceDone(ObjectMesh slice_, std::vector< std::vector< uint32_t > > slice_active_chains_, std::vector< std::vector< uint32_t > > slice_next_chains_) 
{
	vis->peelSliceDone(slice_, slice_active_chains_, slice_next_chains_);
}
void Knittee::firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
    std::vector< std::vector< Stitch > >* active_stitches,
    RowColGraph* graph) 
{  
    vis->firstActiveChainsCreated(active_chains, active_stitches, graph);
}
void Knittee::meshInterpolated(ObjectMesh mesh, std::vector<float> values) 
{
	vis->meshInterpolated(mesh, values);
}
void Knittee::startRemeshing() 
{
    std::vector<Constraint*> constraints = vis->getConstraints();
    knitGrapher.constructNewMesh(constraints);
}
void Knittee::saveConstraints()    //save new constraints to the project folder
{
    std::vector<Constraint*> constraints = vis->getConstraints();

    QString filePath = projectPath + "/constraints";
    QFile file(filePath);
    if (file.open(QIODevice::WriteOnly)) {
        QTextStream out(&file);
        for (Constraint* c : constraints) {
            int size = c->vertices.size();
            if (size > 0) {
                out << c->timeValue << "\n";
                for (int i = 0; i < size - 1; i++) {
                    out << c->vertices[i] << ";";
                }
                out << c->vertices[size - 1];
            }
            out << "\n";
        }
        file.close();
    }
    else {
        qDebug() << "Constraints file open error!";
        return;
    }
}


/*
*  Function: open the 'about' window, displaying some basic information about the project
*/
void Knittee::openAboutWindow()
{
    QDialog* aboutDialog = new QDialog(this); //needs to have 'this'?
    QString aboutText;

    aboutText = HelperText::aboutWindowText;
    aboutDialog->setWindowTitle("About This Software");
    QLabel* label = new QLabel(aboutText, aboutDialog);
    QVBoxLayout* layout = new QVBoxLayout(aboutDialog);
    layout->addWidget(label);
    aboutDialog->setLayout(layout);

    aboutDialog->exec();
}

/*
* Function: open the 'help' window, displaying steps on how to use the software
*/
void Knittee::openHelpWindow()
{
    QDialog* helpDialog = new QDialog(this); //needs to have 'this'?
    QString helpText;

    if (modellingType == 0)
    {
        helpText = HelperText::helpWindowText3D;
    }
    else
    {
        helpText = HelperText::helpWindowText2D;
    }

    //show the help dialog with the relevant help information
    helpDialog->setWindowTitle("Knittee Help");

    QLabel* label = new QLabel(helpText, helpDialog);
    QVBoxLayout* layout = new QVBoxLayout(helpDialog);
    layout->addWidget(label);
    helpDialog->setLayout(layout);

    helpDialog->exec();
}

/*
*  Function: open the new project window, allowing the user to select the project parameters
*/
void Knittee::openNewProjectWindow()
{
    NewProjectDialog dialog;
    QObject::connect(&dialog, &NewProjectDialog::projectConfigurationsSelected, this, &Knittee::handleNewProject);
    dialog.exec();
}
/*
*   Function: handle the new project parameters selected by the user
*/
void Knittee::handleNewProject(ProjectInfo info)
{
    context = info;
    if (context.type == 0) {
        setUpNew3DProject();
    }
    else {
        setUpNew2DProject();
    }
}

/*
*   Function: set up a new 3D project based on the user's input
*/
void Knittee::setUpNew3DProject()
{
    //store the project information in the new folder
    QString projectFolderPath = QCoreApplication::applicationDirPath() + "/Projects/" + context.projectName;
    QDir().mkpath(projectFolderPath);
    QFile file(projectFolderPath + "/" + context.projectName + ".knittee");
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
    projectPath = projectFolderPath;

    start3DProject();
}

/*
*  Function: set up a new 2D project based on the user's input
* 
*/
void Knittee::setUpNew2DProject() 
{
    QString projectFolderPath = QCoreApplication::applicationDirPath() + "/Projects/" + context.projectName;
    QDir().mkpath(projectFolderPath);
    QFile file(projectFolderPath + "/" + context.projectName + ".knittee");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "Error opening project file!";
        return;
    }
    file.write(context.projectName.toUtf8());
    file.write("\n");
    file.write(QString::number(context.type).toUtf8());
    file.write("\n");
    file.write(QString::number(context.width).toUtf8());
    file.write("\n");
	file.write(QString::number(context.height).toUtf8());
    file.write("\n");
    file.write(QString::number(3).toUtf8());
    file.close();

    laceKnitter.createSheet(context.width, context.height, 3);
    laceKnitter.saveToFile(projectFolderPath);

    projectPath = projectFolderPath;
    start2DProject();
}


//destructor, currently empty
Knittee::~Knittee()
{}


//////////////////////////////NewProjectDialog functions/////////////////////////////////////////

/*
*   The following functions are part of a simple 'form' that allow the user to set up a new project, can be intuitively understood so won't be commented
*/
void NewProjectDialog::onPreformRadioClicked() 
{
    meshSelectionButton->setVisible(false);
    meshSelectionLabel->setVisible(false);
    meshSelectionLineEdit->setVisible(false);

    preformSelectionLabel->setVisible(true);
    preformSelectionComboBox->setVisible(true);
}

void NewProjectDialog::onMeshSelectionRadioClicked() 
{
	preformSelectionLabel->setVisible(false);
	preformSelectionComboBox->setVisible(false);

	meshSelectionButton->setVisible(true);
	meshSelectionLabel->setVisible(true);
	meshSelectionLineEdit->setVisible(true);
}

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
    meshRadio = new QRadioButton("3D Mesh editing");
    sheetRadio = new QRadioButton("2D Sheet editing");

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
    */

    //QString projectTypeLabelText = "Select Your Project Type:";
    meshTypeLabel = new QLabel("Project Mesh Type:", this);

    QHBoxLayout* meshTypeRadioLayout = new QHBoxLayout(this);
    preformSelectionRadio = new QRadioButton("Knittee Preforms");
    meshSelectionRadio = new QRadioButton("Custom Mesh");

    meshTypeRadioLayout->addWidget(preformSelectionRadio);
    meshTypeRadioLayout->addWidget(meshSelectionRadio);

    QObject::connect(preformSelectionRadio, &QRadioButton::clicked, this, &NewProjectDialog::onPreformRadioClicked);
    QObject::connect(meshSelectionRadio, &QRadioButton::clicked, this, &NewProjectDialog::onMeshSelectionRadioClicked);


    meshSelectionLabel = new QLabel("Enter the path to the 3D mesh file (.obj):", this);
    meshSelectionLineEdit = new QLineEdit(this);
    meshSelectionButton = new QPushButton("Browse", this);
    meshSelectionLabel->setVisible(false);
    meshSelectionLineEdit->setVisible(false);
    meshSelectionButton->setVisible(false);

    meshSelectionLayout = new QHBoxLayout(this);
    meshSelectionLayout->addWidget(meshSelectionLineEdit);
    meshSelectionLayout->addWidget(meshSelectionButton);
    QObject::connect(meshSelectionButton, &QPushButton::clicked, this, &NewProjectDialog::openUserFile);

    preformSelectionLabel = new QLabel("Select a preform to start with:", this);
    preformSelectionComboBox = new QComboBox(this);

    QString meshFolderPath = QCoreApplication::applicationDirPath() + "/Meshes/";
    QDir meshDirectory(meshFolderPath);
    QStringList files = meshDirectory.entryList(QStringList() << "*.obj", QDir::Files);

    for (QString& file : files) {
        file = QFileInfo(file).baseName();
        preformSelectionComboBox->addItem(file);
    }
    
    QButtonGroup *projectTypeGroup = new QButtonGroup(this);
    QButtonGroup *meshTypeGroup = new QButtonGroup(this);

    //Add radio buttons to their respective groups
    projectTypeGroup->addButton(meshRadio);
    projectTypeGroup->addButton(sheetRadio);

    meshTypeGroup->addButton(preformSelectionRadio);
    meshTypeGroup->addButton(meshSelectionRadio);


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

    mainLayout->addWidget(projectNameLabel);
    mainLayout->addWidget(projectNameLineEdit);
    mainLayout->addWidget(projectTypeLabel);
    mainLayout->addLayout(projectTypeRadioLayout);
    mainLayout->addWidget(meshTypeLabel);
    mainLayout->addLayout(meshTypeRadioLayout);
    mainLayout->addWidget(preformSelectionLabel);
    mainLayout->addWidget(preformSelectionComboBox);
    mainLayout->addWidget(meshSelectionLabel);
    mainLayout->addLayout(meshSelectionLayout);
    mainLayout->addWidget(sheetSizeLabel);
    mainLayout->addLayout(sheetSelectionLayout);
    mainLayout->addLayout(projectConfirmDenyLayout);
    setLayout(mainLayout);
    meshRadio->setChecked(true);
    preformSelectionRadio->setChecked(true);
    onMeshRadioClicked();
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
        if (meshSelectionRadio->isChecked()) {
            projectInfo.objectFilePath = meshSelectionLineEdit->text();
        }
        else {
            projectInfo.objectFilePath = QCoreApplication::applicationDirPath() + "/Meshes/" + preformSelectionComboBox->currentText() + ".obj";
        }
    }
    else if (sheetRadio->isChecked()) {
        projectInfo.type = 1;
        projectInfo.height = heightLineEdit->text().toInt();
        projectInfo.width = widthLineEdit->text().toInt();
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
    sheetSizeLabel->setVisible(false);
    widthLineEdit->setVisible(false);
    heightLineEdit->setVisible(false);
    sheetSizeSeparator->setVisible(false);
    preformSelectionRadio->setVisible(true);
    meshSelectionRadio->setVisible(true);
    meshTypeLabel->setVisible(true);
    preformSelectionRadio->setChecked(true);
    onPreformRadioClicked();
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
    preformSelectionComboBox->setVisible(false);
    preformSelectionLabel->setVisible(false);
    preformSelectionRadio->setVisible(false);
    meshSelectionRadio->setVisible(false);
    meshTypeLabel->setVisible(false);
}


