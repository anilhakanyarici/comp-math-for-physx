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

    qDebug() << mp::Polynomial::bessel(4).toString();

//    mp::Range x = mp::linspace(-4, 4, 1025);
//    mp::Range e = mp::exp(-x.pow(2));
//    for(int i = 1; i < 9; ++i){
//        mp::Range h = mp::Polynomial::hermite(i)(x);
//        mp::curve2d(x, mp::normalize(h * e), QString("n=%1").arg(QString::number(i)));
//    }
//    mp::show2d();

    return app.exec();
}
