#ifndef RUNTEST_H
#define RUNTEST_H

#include "fourier.h"

class RunTest
{
public:
    static void smoothingAnimation();
    static void fourierTransformTest();
    static void fourierSeriesTest();
    static void fourierApproximationAnimation(const mp::Range &x, const mp::Range &y, int iteration = 50, int frame_wait_time_us = 1000000);
};

#endif // RUNTEST_H
