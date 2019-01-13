#include <QDebug>
#include <QApplication>

#include <unistd.h>

#include "plot2d.h"
#include "runtest.h"
#include "polynomial.h"
#include "fourier.h"
#include "videorecorder.h"

#include "Algebra/vector3.h"
#include "Algebra/vector4.h"
#include "Algebra/quaternion.h"
#include "Algebra/matrix3x3.h"
#include "Algebra/matrix4x4.h"
#include "Algebra/matrix.h"


const double planck          = 6.62607004e-30; //kg * cm2 / s
const double speed_of_light  = 29979245800; //cm / s


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    mp::Quaternion q(mp::pi / 4, mp::Vector3::k());
    mp::Matrix3x3 m = mp::Matrix3x3::random();
    std::cout << m.toString() << "\n";
    float *d = m.data()->data();
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j)
            std::cout << std::to_string(d[i + j * 3]) << " \t ";
        std::cout << "\n";
    }




    return 0;
}
