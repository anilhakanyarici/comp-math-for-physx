#include "axis.h"
#include "func.h"

#include <math.h>
#include <stdexcept>

using namespace mp;

struct Axis::pimpl
{
    mp::Range ticks;
    int rank = 0;
    int precision = 0;
    double step = 0;
    double exponent = 1;
    QString name;
};

Axis::Axis()
{
    this->_pimpl = QSharedPointer<pimpl>(new pimpl());
    this->setInterval(0, 0);
}

const mp::Range Axis::ticks() const
{
    return this->_pimpl->ticks;
}
int Axis::precision() const
{
    return this->_pimpl->precision;
}
int Axis::rank() const
{
    return this->_pimpl->rank;
}
double Axis::step() const
{
    return this->_pimpl->step;
}
const QString &Axis::name() const
{
    return this->_pimpl->name;
}
double Axis::exponent() const
{
    return this->_pimpl->exponent;
}

void Axis::setInterval(double a, double b, bool embrace)
{
    double min = fmin(a, b);
    double max = fmax(a, b);
    if(min == -INFINITY || min == INFINITY || max == -INFINITY || max == INFINITY || std::isnan(min) || std::isnan(max))
        throw std::runtime_error("Any interval value cannot be infinity or NaN.");

    double epsilon = 0.000001;

    double abs_min = mp::abs(min);
    int min_rank = 0;
    if(min != 0){
        while(abs_min > 10) { min_rank++; abs_min /= 10; }
        while(abs_min < 1) { min_rank--; abs_min *= 10; }
    }

    double abs_max = mp::abs(max);
    int max_rank = 0;
    if(max != 0){
        while(abs_max > 10) { max_rank++; abs_max /= 10; }
        while(abs_max < 1) { max_rank--; abs_max *= 10; }
    }

    int rank = max_rank;
    if(min_rank > rank)
        rank = min_rank;
    double exp_rank = ::pow(10, rank);

    double min_mantissa = min / exp_rank;
    double max_mantissa = max / exp_rank;

    double delta_mantissa = mp::abs(max_mantissa - min_mantissa);
    double delta_rank = 0;
    if(delta_mantissa == 0) {
        min -= 1;
        max += 1;
        this->setInterval(min, max);
        return;
    } else {
        while(delta_mantissa > 10) { delta_rank++; delta_mantissa /= 10; }
        while(delta_mantissa < 1) { delta_rank--; delta_mantissa *= 10; }
    }
    double exp_delta = ::pow(10, delta_rank);

    double step_list[] { 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5 };
    double tick_step = step_list[0] * exp_delta;

    double tick_min = floor(min_mantissa / tick_step) * tick_step;
    double tick_max = ceil(max_mantissa / tick_step) * tick_step;
    double tick_interval = tick_max - tick_min;
    int tickCount = (1 + round(tick_interval / tick_step)) - ((int)!embrace * 2);

    int attempt = 1;
    while(tickCount > 9) {
        tick_step = step_list[attempt++] * exp_delta;
        tick_min = floor(min_mantissa / tick_step) * tick_step;
        tick_max = ceil(max_mantissa / tick_step) * tick_step;
        tick_interval = tick_max - tick_min;
        tickCount = 1 + round(tick_interval / tick_step) - ((int)!embrace * 2);
    }
    this->_pimpl->ticks = mp::Range(tickCount);

    if(embrace) {
        this->_pimpl->ticks[0] = tick_min;
        for(int i = 1; i < tickCount; i++){
            double tick = i * tick_step + tick_min;
            this->_pimpl->ticks[i] = tick;
        }
    } else {
        for(int i = 0; i < tickCount; i++){
            double tick = (i + 1) * tick_step + tick_min;
            this->_pimpl->ticks[i] = tick;
        }
    }

    this->_pimpl->step = tick_step;
    int prec = 0;
    while(mp::abs(tick_step - (int)tick_step) > epsilon){
        tick_step *= 10;
        prec++;
    }
    this->_pimpl->rank = rank;
    this->_pimpl->precision = prec;
    this->_pimpl->exponent = exp_rank;
}
void Axis::setName(const QString &name)
{
    this->_pimpl->name = name;
}

Axis Axis::clone() const
{
    Axis axis;
    for(int i = 0; i < this->_pimpl->ticks.size(); i++)
        axis._pimpl->ticks.append(this->_pimpl->ticks[i]);
    axis._pimpl->exponent = this->_pimpl->exponent;
    axis._pimpl->name = this->_pimpl->name;
    axis._pimpl->precision = this->_pimpl->precision;
    axis._pimpl->rank = this->_pimpl->rank;
    axis._pimpl->step = this->_pimpl->step;
    return axis;
}
