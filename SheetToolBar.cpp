#include "SheetToolBar.h"


SheetToolBar::SheetToolBar(QWidget* parent) {
    QFrame* rackingSeparator = new QFrame();
    rackingSeparator->setFrameShape(QFrame::HLine);
    rackingSeparator->setFrameShadow(QFrame::Raised);

    QFrame* resizeSeparator = new QFrame();
    resizeSeparator->setFrameShape(QFrame::HLine);
    resizeSeparator->setFrameShadow(QFrame::Raised);

    QFrame* exportSeparator = new QFrame();
    exportSeparator->setFrameShape(QFrame::HLine);
    exportSeparator->setFrameShadow(QFrame::Raised);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QVBoxLayout* resizeLayout = new QVBoxLayout(this);
    QHBoxLayout* leftRightLayout = new QHBoxLayout(this);
    QHBoxLayout* resizeRadioLayout = new QHBoxLayout(this);

    resizeDownButton = new QPushButton("v");
    resizeLeftButton = new QPushButton("<");
    resizeRightButton = new QPushButton(">");
    resizeUpButton = new QPushButton("^");

    increaseRadio = new QRadioButton("Increase");
    decreaseRadio = new QRadioButton("Decrease");
    increaseRadio->setChecked(true);


    AlgorithmComboBox = new QComboBox();
    AlgorithmComboBox->addItem("Collapse-Shift-Expand");
    AlgorithmComboBox->addItem("Schoolbus");

    generateKnitoutButton = new QPushButton("Generate Knitout");
    exportButton = new QPushButton("Export");

    QLabel* resizeLabel = new QLabel("Resize Sheet");
    QLabel* algorithmLabel = new QLabel("Algorithm Selection");
    QLabel* exportLabel = new QLabel("Export Machine Instructions");
    QLabel* sheetLabel = new QLabel("Sheet");

    resizeLabel->setAlignment(Qt::AlignCenter);
    resizeLabel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    algorithmLabel->setAlignment(Qt::AlignCenter);
    algorithmLabel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    exportLabel->setAlignment(Qt::AlignCenter);
    exportLabel->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sheetLabel->setAlignment(Qt::AlignCenter);


    leftRightLayout->addWidget(resizeLeftButton);
    leftRightLayout->addWidget(sheetLabel);
    leftRightLayout->addWidget(resizeRightButton);
    
    resizeLayout->addWidget(resizeUpButton);
    resizeLayout->addLayout(leftRightLayout);
    resizeLayout->addWidget(resizeDownButton);

    resizeRadioLayout->addWidget(increaseRadio);
    resizeRadioLayout->addWidget(decreaseRadio);

    QWidget* spacer = new QWidget();
    spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);


    QLabel* rackingLabel = new QLabel("Maximum Racking:");
    QHBoxLayout* rackingLayout = new QHBoxLayout();
    rackingLayout->addWidget(rackingLabel);
    rackingSpinBox = new QSpinBox();
    rackingSpinBox->setSingleStep(1);
    rackingSpinBox->setMinimum(1);
    rackingSpinBox->setValue(3);
    rackingLayout->addWidget(rackingSpinBox);

    mainLayout->addLayout(rackingLayout);
    mainLayout->addWidget(rackingSeparator);
    mainLayout->addWidget(resizeLabel);
    mainLayout->addLayout(resizeRadioLayout);
    mainLayout->addLayout(resizeLayout);
    mainLayout->addWidget(resizeSeparator);
    mainLayout->addWidget(algorithmLabel);
    mainLayout->addWidget(AlgorithmComboBox);
    mainLayout->addWidget(generateKnitoutButton);
    mainLayout->addWidget(exportSeparator);
    mainLayout->addWidget(exportLabel);
    mainLayout->addWidget(exportButton);
    mainLayout->addWidget(spacer);

    setLayout(mainLayout);
    setMinimumSize(300, 100);

    QObject::connect(resizeUpButton, &QPushButton::clicked, this, &SheetToolBar::upButtonClicked);
    QObject::connect(resizeDownButton, &QPushButton::clicked, this, &SheetToolBar::downButtonClicked);
    QObject::connect(resizeLeftButton, &QPushButton::clicked, this, &SheetToolBar::leftButtonClicked);
    QObject::connect(resizeRightButton, &QPushButton::clicked, this, &SheetToolBar::rightButtonClicked);
    //connect radio buttons to slot
    QObject::connect(increaseRadio, &QRadioButton::clicked, this, &SheetToolBar::increaseRadioClicked);
    QObject::connect(decreaseRadio, &QRadioButton::clicked, this, &SheetToolBar::decreaseRadioClicked);

    QObject::connect(rackingSpinBox, QOverload<int>::of(&QSpinBox::valueChanged), this, &SheetToolBar::rackingSpinBoxChanged);

    QObject::connect(generateKnitoutButton, &QPushButton::clicked, this, &SheetToolBar::generateKnitoutClicked);

}