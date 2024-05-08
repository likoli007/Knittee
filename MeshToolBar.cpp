#include "MeshToolBar.h"

MeshToolBar::MeshToolBar(QWidget* parent) {
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QHBoxLayout* widthLayout = new QHBoxLayout(this);
    QHBoxLayout* heightLayout = new QHBoxLayout(this);
    QHBoxLayout* unitLayout = new QHBoxLayout(this);
    QHBoxLayout* stepLayout = new QHBoxLayout(this);

    widthSpinBox = new QDoubleSpinBox();
    widthSpinBox->setDecimals(2);
    widthSpinBox->setSingleStep(1.0);
    widthSpinBox->setMinimum(0.05);
    widthSpinBox->setMaximum(20.0);
    widthSpinBox->setValue(3.66f);
    widthSpinBox->setSuffix(" mm");

    heightSpinBox = new QDoubleSpinBox();
    heightSpinBox->setDecimals(2);
    heightSpinBox->setSingleStep(1.0);
    heightSpinBox->setMinimum(0.05);
    heightSpinBox->setMaximum(20.0);
    heightSpinBox->setValue(1.73f);
    heightSpinBox->setSuffix(" mm");

    unitSpinBox = new QDoubleSpinBox();
    unitSpinBox->setDecimals(2);
    unitSpinBox->setSingleStep(1.0);
    unitSpinBox->setMinimum(0.05);
    unitSpinBox->setMaximum(20.0);
    unitSpinBox->setValue(10.0f);
    unitSpinBox->setSuffix(" mm");

    QLabel* label = new QLabel("Stitch Width:", this);
    QLabel* label2 = new QLabel("Stitch Height:", this);
    QLabel* label4 = new QLabel("Unit Length:", this);

    widthLayout->addWidget(label);
    widthLayout->addWidget(widthSpinBox);
    heightLayout->addWidget(label2);
    heightLayout->addWidget(heightSpinBox);
    unitLayout->addWidget(label4);
    unitLayout->addWidget(unitSpinBox);


    constraintsModeButton = new QPushButton("Add Constraints", this);
    remeshButton = new QPushButton("Interpolate Mesh", this);
    stepButton = new QPushButton("Step", this);
    stepButton->setDisabled(true);




    autoStepButton = new QPushButton("Auto Step", this);
    autoStepButton->setDisabled(true);

    stepSpinBox = new QSpinBox();
    stepSpinBox->setMinimum(1);
    stepSpinBox->setMaximum(150);
    stepSpinBox->setValue(1);
    stepSpinBox->setSuffix(" steps");

    stepLayout->addWidget(stepSpinBox);
    stepLayout->addWidget(autoStepButton);

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

    QLabel* label3 = new QLabel("Further Options...", this);

    mainLayout->addLayout(widthLayout);
    mainLayout->addLayout(heightLayout);
    mainLayout->addLayout(unitLayout);
    mainLayout->addWidget(constraintsModeButton);
    mainLayout->addWidget(remeshButton);
    mainLayout->addWidget(stepButton);
    mainLayout->addLayout(stepLayout);
    mainLayout->addWidget(label3);
    setLayout(mainLayout);
    setMinimumSize(300, 100);
}