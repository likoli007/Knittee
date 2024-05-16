#include "MeshToolBar.h"

void MeshToolBar::initializeUI() {
    // Set the buttons to be checkable
    showInterpolatedCheck->setDisabled(true);
    showGraphCheck->setDisabled(true);
    showTracedCheck->setDisabled(true);
    showYarnCheck->setDisabled(true);
    generateKnitoutButton->setDisabled(true);
    traceButton->setDisabled(true);
    stepButton->setDisabled(true);
    autoStepButton->setDisabled(true);

}

MeshToolBar::MeshToolBar(QWidget* parent) {

    QFrame* meshSeparator = new QFrame();
    meshSeparator->setFrameShape(QFrame::HLine);
    meshSeparator->setFrameShadow(QFrame::Raised);

    QFrame* sizeSeparator = new QFrame();
    sizeSeparator->setFrameShape(QFrame::HLine);
    sizeSeparator->setFrameShadow(QFrame::Raised);

    QFrame* exportSeparator = new QFrame();
    exportSeparator->setFrameShape(QFrame::HLine);
    exportSeparator->setFrameShadow(QFrame::Raised);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QHBoxLayout* widthLayout = new QHBoxLayout(this);
    QHBoxLayout* heightLayout = new QHBoxLayout(this);
    QHBoxLayout* unitLayout = new QHBoxLayout(this);
    QHBoxLayout* stepLayout = new QHBoxLayout(this);

    widthSpinBox = new QDoubleSpinBox();
    widthSpinBox->setDecimals(2);
    widthSpinBox->setSingleStep(1.0);
    widthSpinBox->setMinimum(0.05);
    widthSpinBox->setMaximum(2000.0);
    widthSpinBox->setValue(3.66f);
    widthSpinBox->setSuffix(" mm");

    heightSpinBox = new QDoubleSpinBox();
    heightSpinBox->setDecimals(2);
    heightSpinBox->setSingleStep(1.0);
    heightSpinBox->setMinimum(0.05);
    heightSpinBox->setMaximum(2000.0);
    heightSpinBox->setValue(1.73f);
    heightSpinBox->setSuffix(" mm");

    unitSpinBox = new QDoubleSpinBox();
    unitSpinBox->setDecimals(2);
    unitSpinBox->setSingleStep(1.0);
    unitSpinBox->setMinimum(0.05);
    unitSpinBox->setMaximum(2000.0);
    unitSpinBox->setValue(10.0f);
    unitSpinBox->setSuffix(" mm");

    //QIcon icon = QIcon::fromTheme("help-about");

    QLabel* label = new QLabel("Stitch Width:", this);
    QLabel* label2 = new QLabel("Stitch Height:", this);
    QLabel* label4 = new QLabel("Unit Length:", this);

    widthLayout->addWidget(label);
    widthLayout->addWidget(widthSpinBox);
    heightLayout->addWidget(label2);
    heightLayout->addWidget(heightSpinBox);
    unitLayout->addWidget(label4);
    unitLayout->addWidget(unitSpinBox);


    constraintsModeButton = new QPushButton("Constraints Mode", this);
    remeshButton = new QPushButton("Interpolate Mesh", this);
    stepButton = new QPushButton("Step", this);




    autoStepButton = new QPushButton("Auto Step", this);
    

    stepSpinBox = new QSpinBox();
    stepSpinBox->setMinimum(1);
    stepSpinBox->setMaximum(500);
    stepSpinBox->setValue(1);
    stepSpinBox->setSuffix(" steps");

    stepLayout->addWidget(stepSpinBox);
    stepLayout->addWidget(autoStepButton);

    traceButton = new QPushButton("Trace Graph", this);
    


    showInterpolatedCheck = new QCheckBox("Mesh");
    showGraphCheck = new QCheckBox("Graph");
    showTracedCheck = new QCheckBox("Traced");
    showYarnCheck = new QCheckBox("Yarn");

    

    // Create a layout for the radio buttons
    QHBoxLayout* showLayout = new QHBoxLayout;
    


    resetButton = new QPushButton("Reset", this);

    QObject::connect(autoStepButton, &QPushButton::clicked, this, &MeshToolBar::onAutoStepButtonClicked);
    QObject::connect(constraintsModeButton, &QPushButton::clicked, this, &MeshToolBar::onConstraintsButtonClicked);
    QObject::connect(remeshButton, &QPushButton::clicked, this, &MeshToolBar::onRemeshButtonClicked);
    QObject::connect(stepButton, &QPushButton::clicked, this, &MeshToolBar::onStepButtonClicked);
    QObject::connect(widthSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
        this, &MeshToolBar::onWidthValueChanged);
    QObject::connect(heightSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
        this, &MeshToolBar::onHeightValueChanged);
    QObject::connect(unitSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
        this, &MeshToolBar::onUnitValueChanged);
    QObject::connect(traceButton, &QPushButton::clicked, this, &MeshToolBar::onTraceButtonClicked);
    QObject::connect(showInterpolatedCheck, &QCheckBox::stateChanged, this, &MeshToolBar::onShowInterpolatedButtonClicked);
    QObject::connect(showGraphCheck, &QCheckBox::stateChanged, this, &MeshToolBar::onShowGraphButtonClicked);
    QObject::connect(showTracedCheck, &QCheckBox::stateChanged, this, &MeshToolBar::onShowTracedButtonClicked);
    QObject::connect(showYarnCheck, &QCheckBox::stateChanged, this, &MeshToolBar::onShowYarnButtonClicked);
    QObject::connect(resetButton, &QPushButton::clicked, this, &MeshToolBar::onResetButtonClicked);
    
  

   QLabel* label3 = new QLabel("Show:", this);
   label3->setStyleSheet("QLabel { padding: 0px; }");
   //QLabel* label = new QLabel("Label");
   // label3->setContentsMargins(0, 0, 0, 0);


   QLabel* parameterLabel = new QLabel("Parameters", this);
   parameterLabel->setAlignment(Qt::AlignCenter);
   //parameterLabel->setStyleSheet("QLabel { padding: 0px; }");
   //parameterLabel->setContentsMargins(0, 0, 0, 0);
   parameterLabel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);

   QLabel* remeshLabel = new QLabel("Remesh Model", this);
   remeshLabel->setAlignment(Qt::AlignCenter);
   remeshLabel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);

   QLabel* buildLabel = new QLabel("Build KnitGraph", this);
   buildLabel->setAlignment(Qt::AlignCenter);
   buildLabel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);

   QLabel* ExportLabel = new QLabel("Export Instructions", this);
   ExportLabel->setAlignment(Qt::AlignCenter);
   ExportLabel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);



    generateKnitoutButton = new QPushButton("Generate Knitout");
   
    QObject::connect(generateKnitoutButton, &QPushButton::clicked, this, &MeshToolBar::onGenerateKnitoutButtonClicked);


    QWidget* spacer = new QWidget();
    spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

    showLayout->addWidget(label3);
    showLayout->addWidget(showInterpolatedCheck);
    showLayout->addWidget(showGraphCheck);
    showLayout->addWidget(showTracedCheck);
    showLayout->addWidget(showYarnCheck);
    showLayout->setContentsMargins(0, 0, 0, 0);
    

    //mainLayout->addWidget(icon);
    mainLayout->addWidget(parameterLabel);
    mainLayout->addLayout(widthLayout);
    mainLayout->addLayout(heightLayout);
    mainLayout->addLayout(unitLayout);
    mainLayout->addWidget(sizeSeparator);
    mainLayout->addWidget(remeshLabel);
    mainLayout->addWidget(constraintsModeButton);
    mainLayout->addWidget(remeshButton);
    mainLayout->addWidget(meshSeparator);
    mainLayout->addWidget(buildLabel);
    mainLayout->addWidget(stepButton);
    mainLayout->addLayout(stepLayout);
    mainLayout->addWidget(traceButton);
    mainLayout->addWidget(exportSeparator);
    mainLayout->addWidget(ExportLabel);
    mainLayout->addWidget(generateKnitoutButton);
    mainLayout->addWidget(resetButton);

    mainLayout->addWidget(spacer);
    //mainLayout->addWidget(label3);
    mainLayout->addLayout(showLayout);
    setLayout(mainLayout);
    setMinimumSize(300, 100);

    initializeUI();
}

void MeshToolBar::knitGraphCreated() {
    traceButton->setDisabled(false);
    showGraphCheck->setDisabled(false);
    showGraphCheck->setChecked(true);
    showInterpolatedCheck->setChecked(false);
}

