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

//    double B1 = 1;
//    double B2 = 2;

//    mp::Range l = mp::linspace(0, 5, 10241);
//    mp::Range E1 = B1 * l * (l + 1);
//    mp::Range E2 = B2 * l * (l + 1);
//    mp::Range f1 = (2 * l) * mp::exp(-E1);
//    mp::Range f2 = (2 * l) * mp::exp(-E2);

//    mp::curve2d(l, f1, "B = 1");
//    mp::curve2d(l, f2, "B = 2");
//    mp::show2d();


//    RunTest::fourierSeriesTest();

    mp::Polynomial p({ -10, 5, 2, -8, 2 });
    qDebug() << p.toString();
    qDebug() << p.derivate().toString();
    qDebug() << p.integrate().toString();


//    mp::Range x = mp::linspace(-2, 4, 1025);
//    mp::Range y2 = mp::poly(x, { -10, 5, 2, -8, 2 });
//    mp::Range y = p(x);
//    mp::curve2d(x, y);
//    mp::curve2d(x, y2);
//    mp::show2d();

    return app.exec();
}
