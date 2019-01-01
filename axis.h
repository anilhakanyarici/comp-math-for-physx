#ifndef AXIS_H
#define AXIS_H

#include <QSharedPointer>
#include <QVector>
#include <QString>

#include "func.h"

namespace mp {
class Axis
{
    struct pimpl;
    QSharedPointer<pimpl> _pimpl;

public:
    Axis();

    const Range ticks() const; //Marks on axis.
    int precision() const; //Maximum digit count (ticks) after comma.
    int rank() const; //All numbers (ticks) can write as k * 10 ^ r. When 1 <= k < 10, r is rank.
    double step() const; //Difference between 2 sequential tick.
    const QString &name() const; //Axis name.
    double exponent() const; //10 power of this->rank();

    void setInterval(double a, double b, bool embrace = true); //Compute and set new ticks embrace a and b.
    void setName(const QString &name);

    Axis clone() const;
};
}

#endif // AXIS_H
