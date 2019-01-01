#include "runtest.h"
#include "plot2d.h"

#include <iostream>
#include <unistd.h>

void RunTest::smoothingAnimation()
{
    mp::Range x = mp::linspace(1, 10, 1025);
    mp::Range y = mp::lb(x) + mp::random(x.size(), 0, 0.7);
    mp::Range s = mp::smooth(y, 2, 1);

    for(int i = 1; i < 1000; ++i) {
        mp::clear2d(1620, 980);
        mp::curve2d(x, y, "func");
        mp::curve2d(x, s, QString("smoothing(iteration:%1)").arg(QString::number(i)));
        mp::show2d();
        s = mp::smooth(s, 2, 1);
        ::usleep(50000);
    }
}

void RunTest::fourierTransformTest()
{
    mp::Range x = mp::linspace(-4, 4, 1001);
    mp::Range y = mp::cos(2 * mp::pi * 3 * x) * mp::exp(-mp::pi * x.pow(2));

    mp::Fourier fourier(x, y);
    mp::Range v = mp::linspace(-6, 6, 1001);
    mp::Range f = fourier.transformReal(v);

    mp::curve2d(x, y);
    mp::curve2d(v, f);

    mp::show2d();
}

void RunTest::fourierSeriesTest()
{
    mp::Range x = mp::linspace(-4, 4, 5001);
    mp::Range y = mp::cos(2 * mp::pi * 3 * x) * mp::exp(-mp::pi * x.pow(2));
    RunTest::fourierApproximationAnimation(x, y, 60, 100000);
}

void RunTest::fourierApproximationAnimation(const mp::Range &x, const mp::Range &y, int iteration, int frame_wait_time_us)
{
    mp::Fourier fourier(x, y);
    mp::Range f = fourier(x);
    double e = mp::mse(f, y);
    qDebug() << "iteration:" << fourier.iteration() << "Mean-Square Error:" << e;
    mp::clear2d(1620, 980);
    mp::curve2d(x, y, "func");
    mp::curve2d(x, f, QString("fourier(iteration:%1) error:%2").arg(QString::number(fourier.iteration()), QString::number(e)));
    mp::show2d();
    for(int i = fourier.iteration(); i < iteration; ++i) {
        ::usleep(frame_wait_time_us);
        fourier.iterate();
        f = fourier(x);
        e = mp::mse(f, y);
        qDebug() << "iteration:" << fourier.iteration() << "Mean-Square Error:" << e;
        mp::clear2d(1620, 980);
        mp::curve2d(x, y, "func");
        mp::curve2d(x, f, QString("fourier(iteration:%1) error:%2").arg(QString::number(fourier.iteration()), QString::number(e)));
        mp::show2d();
    }
}
