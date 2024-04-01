#include "Knittee.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Knittee w;
    w.show();
    return a.exec();
}
