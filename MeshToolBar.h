#pragma once
// This class will be called to be shown when the user selects a 3D mesh project
#include <QWidget>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QSlider>
#include <qobject.h>
#include <qspinbox.h>


//TODO: the following
/*TODO
autostep, edit knitgraph
perhaps a saving of the different steps such that it can be replayed?
generate knitout
output to a machine file
*/


class MeshToolBar : public QWidget {
    Q_OBJECT

public:
    MeshToolBar(QWidget* parent = nullptr);
signals:
    void constraintsButtonClicked();
    void doneButtonClicked();
    void widthChanged(float width);
    void heightChanged(float height);
    void unitChanged(float unit);
    void remeshButtonClicked();
    void stepButtonClicked();
private slots:
    void onStepButtonClicked() {
        qDebug() << "Step button clicked";
        emit stepButtonClicked();
    }
    
    void onWidthValueChanged(float v) {
        qDebug() << "width changed!";
        emit widthChanged(v);
    }
    void onHeightValueChanged(float v) {
        qDebug() << "height changed!";
        emit heightChanged(v);
    }
    void onUnitValueChanged(float v) {
        qDebug() << "unit changed!";
        emit unitChanged(v);
    }

    void onConstraintsButtonClicked() {
        qDebug() << "Constraints button clicked";

        if (constraintChoice != 1) {
            constraintsModeButton->setText("Cancel Constraint Adding");
            constraintChoice = 1;
            emit constraintsButtonClicked();
        }
        else {
            constraintsModeButton->setText("Add Constraints");
            constraintChoice = -1;
            emit doneButtonClicked();
        }

        
    }
    void onRemeshButtonClicked() {
        qDebug() << "Remesh button clicked";
        stepButton->setDisabled(false);
        emit remeshButtonClicked();
    }
private:
    // Will have a slew of functions that pass the user specified options to the visualizer/algorithm program
    QPushButton* constraintsModeButton;
    QPushButton* stepButton;
    QPushButton* remeshButton;


    QDoubleSpinBox* widthSpinBox;
    QDoubleSpinBox* heightSpinBox;
    QDoubleSpinBox* unitSpinBox;


    // Max values for the stitch width and height, will be used to scale the values
    // Might need to change these values later to something more appropriate
    const float widthMin = 1;
    const float widthMax = 10;
    const float heightMin = 1;
    const float heightMax = 10;


    int constraintChoice = 0;

};
