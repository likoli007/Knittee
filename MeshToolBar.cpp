#include "MeshToolBar.h"

MeshToolBar::MeshToolBar(QWidget* parent) {
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    QHBoxLayout* widthLayout = new QHBoxLayout(this);
    QHBoxLayout* heightLayout = new QHBoxLayout(this);

    QLabel* label = new QLabel("Stitch Width:", this);
    widthLayout->addWidget(label);

    widthSlider = new QSlider(Qt::Horizontal);

    widthSlider->setRange(0, 100);
    widthSlider->setSingleStep(1);


    //widthSlider->valueChanged->connect(valueHandler);
    //widthSlider->setMinimum(1);  // Set the minimum value
    //swidthSlider->setMaximum(10);  // Set the maximum value
    widthSlider->setValue(5);  // Set the initial value
    widthLayout->addWidget(widthSlider);

    QLabel* label2 = new QLabel("Stitch Height:", this);
    heightLayout->addWidget(label2);

    heightSlider = new QSlider(Qt::Horizontal);
    heightSlider->setRange(0, 100);
    heightSlider->setSingleStep(1);
    //heightSlider->setMinimum(1);  // Set the minimum value
    //heightSlider->setMaximum(10);  // Set the maximum value
    heightSlider->setValue(5);  // Set the initial value
    heightLayout->addWidget(heightSlider);


    QLabel* label4 = new QLabel("Unit Length:", this);
    heightLayout->addWidget(label4);
    unitSlider = new QSlider(Qt::Horizontal);
    unitSlider->setRange(0, 100);
    unitSlider->setSingleStep(1);
    //unitSlider->setMinimum(1);  // Set the minimum value
    //unitSlider->setMaximum(10);  // Set the maximum value
    unitSlider->setValue(5);  // Set the initial value
    heightLayout->addWidget(unitSlider);


    constraintsModeButton = new QPushButton("Add Constraints", this);
    QPushButton* interpolateButton = new QPushButton("Interpolate", this);
    //interpolateButton->setEnabled(false);
    doneButton = new QPushButton("Done", this);

    remeshButton = new QPushButton("Remesh", this);


    stepButton = new QPushButton("Step", this);


    //QObject::connect(constraintsModeButton, SIGNAL(clicked()), parent, SLOT(setConstraintsMode()));
    //connect(button, &QPushButton::clicked, this, &SubClass::onButtonClicked);
    QObject::connect(constraintsModeButton, &QPushButton::clicked, this, &MeshToolBar::onConstraintsButtonClicked);
    QObject::connect(doneButton, &QPushButton::clicked, this, &MeshToolBar::onDoneButtonClicked);
    QObject::connect(widthSlider, &QSlider::sliderReleased, this, &MeshToolBar::onWidthSliderReleased);
    QObject::connect(heightSlider, &QSlider::sliderReleased, this, &MeshToolBar::onHeightSliderReleased);
    QObject::connect(unitSlider, &QSlider::sliderReleased, this, &MeshToolBar::onUnitSliderReleased);
    QObject::connect(remeshButton, &QPushButton::clicked, this, &MeshToolBar::onRemeshButtonClicked);
    QObject::connect(stepButton, &QPushButton::clicked, this, &MeshToolBar::onStepButtonClicked);

    QLabel* label3 = new QLabel("Further Options...", this);

    mainLayout->addLayout(widthLayout);
    mainLayout->addLayout(heightLayout);
    mainLayout->addWidget(constraintsModeButton);
    mainLayout->addWidget(interpolateButton);
    mainLayout->addWidget(doneButton);
    mainLayout->addWidget(remeshButton);
    mainLayout->addWidget(stepButton);
    mainLayout->addWidget(label3);
    setLayout(mainLayout);
    setMinimumSize(300, 100);
}