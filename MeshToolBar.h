#pragma once
// This class will be called to be shown when the user selects a 3D mesh project
#include <QWidget>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QSlider>
#include <qobject.h>

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
    
    void onConstraintsButtonClicked() {
        qDebug() << "Constraints button clicked";
        emit constraintsButtonClicked();
    }
    void onDoneButtonClicked() {
        qDebug() << "Done button clicked";
        emit doneButtonClicked();
    }
    void onRemeshButtonClicked() {
        qDebug() << "Remesh button clicked";
        emit remeshButtonClicked();
    }
    void onWidthSliderReleased() {
        qDebug() << "Width slider released";

        int value = widthSlider->value();

        // Map the value from 0 to 100 as a value of float from 1.0f to 6.0f
        float width = (value / 100.0f) * 5.0f + 1.0f;

        emit widthChanged(width);
    }
    void onHeightSliderReleased() {
        qDebug() << "Height slider released";
        int value = heightSlider->value();

        // Map the value from 0 to 100 as a value of float from 1.0f to 6.0f
        float height = (value / 100.0f) * 5.0f + 1.0f;

        emit heightChanged(height);

    }
    void onUnitSliderReleased() {
        qDebug() << "Unit slider released";
        int value = unitSlider->value();

        // Map the value from 0 to 100 as a value of float from 1.0f to 6.0f
        float unit = (value / 100.0f) * 5.0f + 1.0f;

        emit unitChanged(unit);

    }

private:
    // Will have a slew of functions that pass the user specified options to the visualizer/algorithm program
    QPushButton* constraintsModeButton;
    QPushButton* stepButton;
    QPushButton* doneButton;
    QSlider* widthSlider;
    QSlider* heightSlider;
    QSlider* unitSlider;
    QPushButton* remeshButton;

    // Max values for the stitch width and height, will be used to scale the values
    // Might need to change these values later to something more appropriate
    const float widthMin = 1;
    const float widthMax = 10;
    const float heightMin = 1;
    const float heightMax = 10;

};
