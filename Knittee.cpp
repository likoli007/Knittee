#include "Knittee.h"
#include "HelperText.h"

Knittee::Knittee(QWidget* parent)
    : QMainWindow(parent)
{

    ui.setupUi(this);
    setFocusPolicy(Qt::StrongFocus);

    auto topMenu = this->menuBar();//new QMenuBar(this);

    auto projectMenu = new QMenu(topMenu); //project menu: user can open a new project (2D or 3D), load an existing one etc.
    //auto optionsMenu = new QMenu(topMenu); //options menu: will be able to select some basic options regarding the whole app
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


    //optionsMenu->setTitle("Options");

    meshToolsWidget = new MeshToolBar(this);
    sheetToolsWidget = new SheetToolBar(this);

    topMenu->addMenu(projectMenu);
    //topMenu->addMenu(optionsMenu);
    topMenu->addAction(helpAction);
    topMenu->addAction(aboutAction);

    setMenuBar(topMenu);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);

    visualizerLayout = new QHBoxLayout(this);
    QVBoxLayout* optionsLayout = new QVBoxLayout(this);


    vis = new Visualizer(this); // Set the parent to be the main window
    vis->setMinimumSize(800, 600); // Set a minimum size for the widget
    visualizerLayout->addLayout(optionsLayout);
    visualizerLayout->addWidget(vis);

    QObject::connect(newProjectAction, &QAction::triggered, this, &Knittee::openNewProjectWindow);
    QObject::connect(loadAction, &QAction::triggered, this, &Knittee::selectProject);
    QObject::connect(exitAction, &QAction::triggered, this, &QApplication::quit);
    QObject::connect(saveAction, &QAction::triggered, this, &Knittee::saveProject);
    //QObject::connect(optionsMenu, &QMenu::triggered, this, &Knittee::openOptionsWindow);
    QObject::connect(helpAction, &QAction::triggered, this, &Knittee::openHelpWindow);
    QObject::connect(aboutAction, &QAction::triggered, this, &Knittee::openAboutWindow);
    QObject::connect(meshToolsWidget, SIGNAL(constraintsButtonClicked(bool)), vis, SLOT(setConstraintsMode(bool)));
    //QObject::connect(toolsWidget, SIGNAL(doneButtonClicked()), this, SLOT(handleToolbarDone()));
    QObject::connect(meshToolsWidget, SIGNAL(remeshButtonClicked()), this, SLOT(startRemeshing()));
    QObject::connect(meshToolsWidget, SIGNAL(widthChanged(float)), &knitGrapher, SLOT(setStitchWidth(float)));
    QObject::connect(meshToolsWidget, SIGNAL(heightChanged(float)), &knitGrapher, SLOT(setStitchHeight(float)));
    QObject::connect(meshToolsWidget, SIGNAL(unitChanged(float)), &knitGrapher, SLOT(setModelUnitLength(float)));
    QObject::connect(meshToolsWidget, SIGNAL(stepButtonClicked()), &knitGrapher, SLOT(stepButtonClicked()));
    QObject::connect(meshToolsWidget, SIGNAL(traceButtonClicked()), &knitGrapher, SLOT(traceButtonClicked()));
    
    QObject::connect(vis, SIGNAL(requestConstraintsSave()), this, SLOT(saveConstraints()));
    



    QObject::connect(&knitGrapher, SIGNAL(knitGraphInterpolated(ObjectMesh, std::vector<float>)), this, SLOT(meshInterpolated(ObjectMesh, std::vector<float>)));
    QObject::connect(&knitGrapher, SIGNAL(firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >*, std::vector< std::vector< Stitch > >*, RowColGraph*)), 
        this, SLOT(firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >*, std::vector< std::vector< Stitch > >*, RowColGraph*)));
    QObject::connect(&knitGrapher, SIGNAL(peelSliceDone(ObjectMesh*, std::vector< std::vector< uint32_t > >*, std::vector< std::vector< uint32_t > >*)), this, SLOT(peelSliceDone(ObjectMesh*, std::vector< std::vector< uint32_t > > *, std::vector< std::vector< uint32_t > >*)));
    QObject::connect(&knitGrapher, SIGNAL(linkChainsDone(std::vector< std::vector< Stitch > >*, std::vector< Link >*)), this, SLOT(linkChainsDone(std::vector< std::vector< Stitch > >*, std::vector< Link >*)));
    QObject::connect(&knitGrapher, SIGNAL(nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >*)), this, SLOT(nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >*)));
    QObject::connect(&knitGrapher, SIGNAL(knitGraphCreated()), this, SLOT(knitGraphCreated()));
   
    QObject::connect(&knitGrapher, SIGNAL(knitGraphTraced(std::vector< TracedStitch >*)), this, SLOT(knitGraphTraced(std::vector< TracedStitch >*)));
    QObject::connect(&knitGrapher, SIGNAL(knitGraphTraced(std::vector< TracedStitch >*)), meshToolsWidget, SLOT(knitGraphTraced()));
    QObject::connect(&knitGrapher, SIGNAL(knitGraphTraced(std::vector< TracedStitch >*)), vis, SLOT(knitGraphTraced(std::vector< TracedStitch >*)));
    //QObject::connect(&knitGrapher, SIGNAL)
    QObject::connect(meshToolsWidget, SIGNAL(showInterpolatedChanged(int)), vis, SLOT(showInterpolatedChanged(int)));
    QObject::connect(meshToolsWidget, SIGNAL(showGraphChanged(int)), vis, SLOT(showGraphChanged(int)));
    QObject::connect(meshToolsWidget, SIGNAL(showTracedChanged(int)), vis, SLOT(showTracedChanged(int)));


    QObject::connect(meshToolsWidget, SIGNAL(helpBoxCommunication(QString)),this, SLOT(helpBoxCommunication(QString)));
    QObject::connect(&knitGrapher, SIGNAL(helpBoxCommunication(QString)), this, SLOT(helpBoxCommunication(QString)));

    QObject::connect(meshToolsWidget, SIGNAL(resetButtonClicked()), this, SLOT(resetButtonClicked()));

    //QObject::connect(sheetToolsWidget, SIGNAL(sheetSizeChanged(int, int)), &laceKnitter, SLOT(createSheet(int, int)));

    QObject::connect(sheetToolsWidget, SIGNAL(rackingChanged(int)), &laceKnitter, SLOT(rackingChanged(int)));
    QObject::connect(sheetToolsWidget, SIGNAL(widthChanged(int, int)), &laceKnitter, SLOT(widthChanged(int, int)));
    QObject::connect(sheetToolsWidget, SIGNAL(heightChanged(int, int)), &laceKnitter, SLOT(heightChanged(int, int)));

    QObject::connect(vis, SIGNAL(moveLoop(QPair<int, int>, QPair<int, int>)), &laceKnitter, SLOT(moveLoop(QPair<int, int>, QPair<int, int>)));
    QObject::connect(&laceKnitter, SIGNAL(sheetChanged(std::vector<std::vector<FlatPoint>>*)), vis, SLOT(sheetChanged(std::vector<std::vector<FlatPoint>>*)));
    QObject::connect(&laceKnitter, SIGNAL(sheetDimensionsChanged(int, int, int)), this, SLOT(sheetDimensionsChanged(int, int, int)));
    
    QObject::connect(&knitoutScheduler, SIGNAL(instructionsCreated(std::vector<std::string>)), this, SLOT(instructionsCreated(std::vector<std::string>)));
    QObject::connect(&knitoutScheduler, SIGNAL(knitoutGenerated(std::vector<QString>)), this, SLOT(knitoutGenerated(std::vector<QString>)));


    QObject::connect(meshToolsWidget, SIGNAL(generateKnitoutButtonClicked()), this, SLOT(generateKnitoutButtonClicked()));

    messageTextEdit = new QTextEdit(this);
    messageTextEdit->setReadOnly(true);
    messageTextEdit->setMinimumSize(800, 100);
    messageTextEdit->setText("Welcome to Knittee! Please open a project to get started.");

    mainLayout->addLayout(visualizerLayout);
    mainLayout->addWidget(messageTextEdit);

    mainLayout->setStretchFactor(visualizerLayout, 4);
    mainLayout->setStretchFactor(messageTextEdit, 1);

    QWidget* centralWidget = new QWidget(this);
    centralWidget->setLayout(mainLayout);

    setCentralWidget(centralWidget);
}

void Knittee::knitoutGenerated(std::vector<QString> knitout) {
    qDebug() << "instructions created!";

    //save to project folder under 'knitoutgenerator'
    QString filePath = projectPath + "/knitout";
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

void Knittee::generateKnitoutButtonClicked() {
	qDebug() << "generate knitout button clicked!";
	knitoutScheduler.schedule(projectPath + "/traced");
}

void Knittee::instructionsCreated(std::vector<std::string> instructions) {
	qDebug() << "instructions created!";
	
    //save to project folder under 'knitoutgenerator'
    QString filePath = projectPath + "/knitoutgenerator";
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

void Knittee::closeEvent(QCloseEvent* event) {
    
    saveProject();
}

void Knittee::save3DProject() {
    qDebug() << "saving later yo";
}

void Knittee::saveProject(){
    qDebug() << "savingyo";
	if (modellingType == 0) {
		save3DProject();
	}
	else if (modellingType == 1){
		save2DProject();
	}
}

void Knittee::save2DProject() {
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

void Knittee::resetButtonClicked() {
	vis->reset();
	knitGrapher.reset();
}

void Knittee::helpBoxCommunication(QString message) {
	messageTextEdit->setText(message);
}

void Knittee::knitGraphTraced(std::vector< TracedStitch >* traced_stitches) {
	//vis->knitGraphTraced(traced_stitches);

    //qDebug() << "saving knitgraph to file";

    saveTraced(traced_stitches);

    //qDebug() << "sending traced to visualizer";
}

void Knittee::saveTraced(std::vector< TracedStitch >* traced_stitches) {
    qDebug() << "Saving, traced size: " << traced_stitches->size();



	QString filePath = projectPath + "/traced";
	QFile file(filePath);

    if (!file.open(QIODevice::WriteOnly))
    {
        qDebug() << "could not open file";
        return;
    }
   
    QTextStream out(&file);
    //KnitOutStitch* stitch;

    for (auto const& ts : *traced_stitches)
    {
        //stitch = new KnitOutStitch();
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

void Knittee::knitGraphCreated()
{
    meshToolsWidget->knitGraphCreated();
}

void Knittee::nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >* active_chains) {
	qDebug() << "knittee caught the nextactivechains emit!";
	vis->nextActiveChainsDone(active_chains);
}

void Knittee::linkChainsDone(std::vector< std::vector< Stitch > >* next_stitches, std::vector< Link >* links) {
    qDebug() << "knittee caught the linkchainsdone emits!";
	vis->linkChainsDone(next_stitches, links);
}

void Knittee::peelSliceDone(ObjectMesh* slice_, std::vector< std::vector< uint32_t > >* slice_active_chains_, std::vector< std::vector< uint32_t > >* slice_next_chains_) {
	vis->peelSliceDone(slice_, slice_active_chains_, slice_next_chains_);

}

void Knittee::firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
    std::vector< std::vector< Stitch > >* active_stitches,
    RowColGraph* graph) {
    qDebug() << "knittee caught the firstachtivechains emit!";
    vis->firstActiveChainsCreated(active_chains, active_stitches, graph);

}

void Knittee::meshInterpolated(ObjectMesh mesh, std::vector<float> values) {
	vis->meshInterpolated(mesh, values);
    canSlice = true;
}

void Knittee::selectProject() {
    QString filePath = QFileDialog::getOpenFileName(this, "Open Knittee Project File", "", "Custom Files (*.knittee)");
    QFileInfo fileInfo(filePath);
    projectPath = fileInfo.absolutePath();
    loadProject(filePath);
}

void Knittee::startRemeshing() {
    messageTextEdit->setText("Remeshing caught...");
    std::vector<Constraint*> constraints = vis->getConstraints();

    knitGrapher.constructNewMesh(constraints);
}


void Knittee::saveConstraints()
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

void Knittee::handleToolbarDone() {
    //currently nothing but constraints to handle, but later on user clicking done may mean done constraints,
    //done interpolation, done editing etc..
    qDebug() << "handling toolbar done!";
    saveConstraints();

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

void Knittee::loadConstraints() {

    qDebug() << "loading constraints";
    std::vector<Constraint*> constraints;
    qDebug() << "project path: " << projectPath;
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

    visualizerLayout->removeItem(visualizerLayout->itemAt(0));
    visualizerLayout->insertWidget(0, meshToolsWidget);
    visualizerLayout->setStretchFactor(vis, 3); 
    visualizerLayout->setStretchFactor(meshToolsWidget, 1);

   

}


void Knittee::start2DProject()
{
    qDebug() << "starting 2D project: " << context.width << "x" << context.height;
    vis->projectType = 1;
    modellingType = 1;

    laceKnitter.setDimensions(context.width, context.height, context.racking);
    laceKnitter.loadFromFile(projectPath);

    visualizerLayout->removeItem(visualizerLayout->itemAt(0));
    visualizerLayout->insertWidget(0, sheetToolsWidget);
    visualizerLayout->setStretchFactor(vis, 3);
    visualizerLayout->setStretchFactor(sheetToolsWidget, 1);
}

void Knittee::openOptionsWindow()
{
    QMainWindow* optionsWindow = new QMainWindow(this);
    optionsWindow->setAttribute(Qt::WA_DeleteOnClose);
    optionsWindow->show();
}

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


void Knittee::openNewProjectWindow()
{
    NewProjectDialog dialog;
    QObject::connect(&dialog, &NewProjectDialog::projectConfigurationsSelected, this, &Knittee::handleNewProject);
    dialog.exec();
}


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



void Knittee::setUpNew2DProject() {
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

void Knittee::handleNewProject(ProjectInfo info)
{
    context = info;
    if (context.type == 0) {
        qDebug() << "Setting up 3D editing project...";
        setUpNew3DProject();
        //openFile(options.filePath);
    }
    else {
        qDebug() << "Setting up 2D editing project...";
        setUpNew2DProject();
    }

    // Handle selected options
    // options.type, options.file, etc. contain the selected values
}

void Knittee::setConstraintsMode()
{
    //vis->setConstraintsMode();
}


void Knittee::KeyPressEvent(QKeyEvent* event)
{
    qDebug() << "Key pressesssd: " << event->key();
    if (event->key() == Qt::Key_N) {
        qDebug() << "pressed N";
        knitGrapher.stepButtonClicked();
    }
}

Knittee::~Knittee()
{}


void NewProjectDialog::onPreformRadioClicked() {
    meshSelectionButton->setVisible(false);
    meshSelectionLabel->setVisible(false);
    meshSelectionLineEdit->setVisible(false);

    preformSelectionLabel->setVisible(true);
    preformSelectionComboBox->setVisible(true);
}

void NewProjectDialog::onMeshSelectionRadioClicked() {
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
    * TODO: add other options once the 3D portion grows more robust
    *       add a 'preform' option, package some basic .obj files with the program
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
    //QStringList files = meshDirectory.entryList(QStringList() << "*.obj", QDir::Files);
    for (QString& file : files) {
        file = QFileInfo(file).baseName();
        preformSelectionComboBox->addItem(file);
    }
    
    QButtonGroup *projectTypeGroup = new QButtonGroup(this);
    QButtonGroup *meshTypeGroup = new QButtonGroup(this);

    // Add radio buttons to their respective groups
    projectTypeGroup->addButton(meshRadio);
    projectTypeGroup->addButton(sheetRadio);

    meshTypeGroup->addButton(preformSelectionRadio);
    meshTypeGroup->addButton(meshSelectionRadio);

    //projectTypeGroup->setExclusive(false);
    //meshTypeGroup->setExclusive(false);

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
        qDebug() << "2D project starting, selected size: " << widthLineEdit->text() << "/" << heightLineEdit->text();
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
    //meshSelectionLabel->setVisible(true);
    //meshSelectionLineEdit->setVisible(true);
    //meshSelectionButton->setVisible(true);
    sheetSizeLabel->setVisible(false);
    widthLineEdit->setVisible(false);
    heightLineEdit->setVisible(false);
    sheetSizeSeparator->setVisible(false);
    preformSelectionRadio->setVisible(true);
    meshSelectionRadio->setVisible(true);

    //preformSelectionRadio->setVisible(false);
    //meshSelectionRadio->setVisible(false);
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


