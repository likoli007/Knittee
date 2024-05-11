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
#include <qtimer.h>
#include <qeventloop.h>
#include <qradiobutton.h>
#include <qcheckbox.h>


#include "HelperText.h"

class MeshToolBar : public QWidget {
    Q_OBJECT

public:
    MeshToolBar(QWidget* parent = nullptr);
    void knitGraphCreated();
    void initializeUI();
signals:
    void constraintsButtonClicked(bool);
    void doneButtonClicked();
    void widthChanged(float width);
    void heightChanged(float height);
    void unitChanged(float unit);
    void remeshButtonClicked();
    void stepButtonClicked();
    void traceButtonClicked();
    void showInterpolatedChanged(int state);
    void showGraphChanged(int state);
    void showTracedChanged(int state);
    void helpBoxCommunication(QString message);
    void resetButtonClicked();
    void generateKnitoutButtonClicked();
private slots:
    void knitGraphTraced() {
        generateKnitoutButton->setDisabled(false);
        showTracedCheck->setDisabled(false);
        showGraphCheck->setChecked(false);
        showInterpolatedCheck->setChecked(false);
        showTracedCheck->setChecked(true);
    }

    void onAutoStepButtonClicked() {
        stepButton->setDisabled(true);
        autoStepButton->setDisabled(true);
        qDebug() << "Auto Step button clicked";
        for (int i = 0; i < stepSpinBox->value(); i++) {
			emit stepButtonClicked();
            QEventLoop loop;
            QTimer::singleShot(200, &loop, &QEventLoop::quit); // Wait for 1000 milliseconds (1 second)
            loop.exec();
		}
        stepButton->setDisabled(false);
        autoStepButton->setDisabled(false);
    }

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
        if (constraintChoice != 1) {
            constraintsModeButton->setText("Cancel Constraints Mode");
            constraintChoice = 1;
            emit helpBoxCommunication(HelperText::constraintText);
            emit constraintsButtonClicked(true);
        }
        else {
            constraintsModeButton->setText("Constraints Mode");
            constraintChoice = -1;
            emit helpBoxCommunication(HelperText::constraintEndText);
            emit constraintsButtonClicked(false);
        }

        
    }
    void onRemeshButtonClicked() {
        emit helpBoxCommunication(HelperText::interpolateText);
        stepButton->setDisabled(false);
        autoStepButton->setDisabled(false);
        showInterpolatedCheck->setDisabled(false);
        showInterpolatedCheck->setChecked(true);
        
        emit remeshButtonClicked();
    }
    void onTraceButtonClicked() {
        emit traceButtonClicked();
    }

    void onShowInterpolatedButtonClicked(int state) {
        emit showInterpolatedChanged(state);
	}
    void onShowGraphButtonClicked(int state) {
        emit showGraphChanged(state);
	}
    void onShowTracedButtonClicked(int state) {
        emit showTracedChanged(state);
    }

    void onResetButtonClicked() {
		emit helpBoxCommunication(HelperText::resetText);
        initializeUI();
        constraintsModeButton->setText("Constraints Mode");
        constraintChoice = 0;
        showInterpolatedCheck->setChecked(false);
        showTracedCheck->setChecked(false);
        showGraphCheck->setChecked(false);
		emit resetButtonClicked();
	}
    void onGenerateKnitoutButtonClicked() {
		//emit helpBoxCommunication(HelperText::knitoutText);
		emit generateKnitoutButtonClicked();
	}
private:
    // Will have a slew of functions that pass the user specified options to the visualizer/algorithm program
    QPushButton* constraintsModeButton;
    QPushButton* stepButton;
    QPushButton* remeshButton;
    QPushButton* traceButton;

    QDoubleSpinBox* widthSpinBox;
    QDoubleSpinBox* heightSpinBox;
    QDoubleSpinBox* unitSpinBox;


    QSpinBox* stepSpinBox;
    QPushButton* autoStepButton;

    QPushButton* generateKnitoutButton;
    
    QPushButton* resetButton;

    QCheckBox* showInterpolatedCheck;
    QCheckBox* showGraphCheck;
    QCheckBox* showTracedCheck;

    // Max values for the stitch width and height, will be used to scale the values
    // Might need to change these values later to something more appropriate
    const float widthMin = 1;
    const float widthMax = 10;
    const float heightMin = 1;
    const float heightMax = 10;


    int constraintChoice = 0;



};
