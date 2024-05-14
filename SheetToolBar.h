#pragma once
#include <QWidget>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QSlider>
#include <qradiobutton.h>
#include <qcombobox.h>
#include <qspinbox.h>

class SheetToolBar : public QWidget
{
	Q_OBJECT

private slots:
	void upButtonClicked() {
		if (adding == 1) {
			emit heightChanged(1, 0) ;
		}
		else {
			emit heightChanged(-1, 0);
		}
	}

	void downButtonClicked() {
		if (adding == 1) {
			emit heightChanged(1, 1);
		}
		else {
			emit heightChanged(-1, 1);
		}
	}

	void leftButtonClicked() {
		if (adding == 1) {
			emit widthChanged(1, 1);
		}
		else {
			emit widthChanged(-1, 1);
		}
	}

	void rightButtonClicked() {
		if (adding == 1) {
			emit widthChanged(1, 0);
		}
		else {
			emit widthChanged(-1, 0);
		}
	}

	void increaseRadioClicked() {
		adding = 1;
		resizeUpButton->setText("^");
		resizeDownButton->setText("v");
		resizeLeftButton->setText("<");
		resizeRightButton->setText(">");
	}
	void decreaseRadioClicked() {
		adding = 0;
		resizeUpButton->setText("v");
		resizeDownButton->setText("^");
		resizeLeftButton->setText(">");
		resizeRightButton->setText("<");
	}

	void rackingSpinBoxChanged() {
		int value = rackingSpinBox->value();
		emit rackingChanged(value);
	}

	void generateKnitoutClicked() {
		int algorithm = AlgorithmComboBox->currentIndex();
		emit generateKnitoutSheet(algorithm);
	}

signals:
	void generateKnitoutSheet(int algorithm);
	void rackingChanged(int value);
	void heightChanged(int height, int side);
	void widthChanged(int width, int side);


public:
	SheetToolBar(QWidget* parent = nullptr);


private:
	
	QPushButton* resizeUpButton;
	QPushButton* resizeDownButton;
	QPushButton* resizeLeftButton;
	QPushButton* resizeRightButton;

	QRadioButton* increaseRadio;
	QRadioButton* decreaseRadio;

	QComboBox* AlgorithmComboBox;

	QPushButton* generateKnitoutButton;

	QPushButton* exportButton;

	QSpinBox* rackingSpinBox;

	int adding = 1;
};

