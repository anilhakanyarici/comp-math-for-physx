#include <QDebug>
#include <QApplication>

#include <unistd.h>

#include "plot2d.h"
#include "runtest.h"
#include "polynomial.h"
#include "fourier.h"

const double planck          = 6.62607004e-30; //kg * cm2 / s
const double speed_of_light  = 29979245800; //cm / s


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    qDebug() << mp::Polynomial::hermite(6).toString();

    return app.exec();
}
