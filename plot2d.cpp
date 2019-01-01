#include "plot2d.h"

#include <QApplication>

namespace mp {
mp::Plot2D *plot2d_gInstance = nullptr;

void curve2d(const mp::Range &x, const mp::Range &y, const QString &name)
{
    if(mp::plot2d_gInstance == nullptr){
        mp::plot2d_gInstance = new mp::Plot2D();
        mp::plot2d_gInstance->resize(800, 600);
    }
    mp::plot2d_gInstance->addCurve(x, y, name);
}

void show2d()
{
    if(mp::plot2d_gInstance == nullptr){
        mp::plot2d_gInstance = new mp::Plot2D();
        mp::plot2d_gInstance->resize(800, 600);
    }
    mp::plot2d_gInstance->repaint();
    QApplication::processEvents();
    mp::plot2d_gInstance->show();
    QApplication::processEvents();
}

void dot2d(const Range &x, const Range &y, const double &rad, const QString &name)
{
    if(mp::plot2d_gInstance == nullptr){
        mp::plot2d_gInstance = new mp::Plot2D();
        mp::plot2d_gInstance->resize(800, 600);
    }
    mp::plot2d_gInstance->addDot(x, y, rad, name);
}

void clear2d(int w, int h)
{
    if(mp::plot2d_gInstance == nullptr){
        mp::plot2d_gInstance = new mp::Plot2D();
        mp::plot2d_gInstance->resize(800, 600);
    }
    mp::plot2d_gInstance->clear();
    mp::plot2d_gInstance->resize(w, h);
}

QPixmap drawnPixmap2d()
{
    if(mp::plot2d_gInstance == nullptr)
        return QPixmap();
    return mp::plot2d_gInstance->drawnPixmap();
}

}

mp::Plot2D::Plot2D(QWidget *parent) : QLabel(parent) {
    this->graph.setEmbrace(false);
    this->resize(800, 600);
}

void mp::Plot2D::addCurve(const mp::Range &x, const mp::Range &y, const QString &name) {
    this->graph.curve(x, y, name);
}

void mp::Plot2D::addDot(const mp::Range &x, const mp::Range &y, const double &rad, const QString &name)
{
    this->graph.dot(x, y, rad, name);
}

void mp::Plot2D::redraw() {
    this->graph.calculate(this->size());
    this->pixmap = QPixmap(this->size());
    this->graph.draw(&this->pixmap);
    this->setPixmap(this->pixmap);
}

void mp::Plot2D::clear()
{
    this->graph = GraphDrawer2D();
}

void mp::Plot2D::resizeEvent(QResizeEvent *e) {
    this->QLabel::resizeEvent(e);
    this->redraw();
}

void mp::Plot2D::paintEvent(QPaintEvent *e)
{
    this->QLabel::paintEvent(e);
    this->redraw();
}

void mp::Plot2D::curve2d(const mp::Range &x, const mp::Range &y, const QString &name)
{
    this->addCurve(x, y, name);
}

void mp::Plot2D::dot2d(const mp::Range &x, const mp::Range &y, const double &rad, const QString &name)
{
    this->addDot(x, y, rad, name);
}

void mp::Plot2D::clear2d(int w, int h)
{
    this->clear();
    this->resize(w, h);
}

void mp::Plot2D::show2d()
{
    this->repaint();
    QApplication::processEvents();
    this->show();
    QApplication::processEvents();
}
