#include <QDebug>
#include <QApplication>

#include <unistd.h>

#include "plot2d.h"
#include "runtest.h"



int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    RunTest::fourierSeriesTest();

    return app.exec();
}
