#ifndef PLOT2D_H
#define PLOT2D_H

#include <QLabel>
#include <QDebug>

#include "graphdrawer2d.h"

namespace mp {
class Plot2D : public QLabel
{
    Q_OBJECT

    GraphDrawer2D graph;
    QPixmap pixmap;
public:
    explicit Plot2D(QWidget *parent = nullptr);

    const QPixmap &drawnPixmap() const { return this->pixmap; }

    QString apsisName() { return this->graph.xAxis().name(); }
    void setApsisName(const QString &name) { this->graph.xAxis().setName(name); }
    QString ordinateName() { return this->graph.yAxis().name(); }
    void setOrdinateName(const QString &name) { this->graph.yAxis().setName(name); }

    void addCurve(const Range &x, const Range &y, const QString &name = "");
    void addDot(const Range &x, const Range &y, const double &rad, const QString &name = "");
    void redraw();
    void clear();

    virtual void resizeEvent(QResizeEvent *e) override;
    virtual void paintEvent(QPaintEvent *e) override;

    void curve2d(const Range &x, const Range &y, const QString &name = "");
    void dot2d(const Range &x, const Range &y, const double &rad = 4, const QString &name = "");
    void clear2d(int w = 800, int h = 600);
    void show2d();
};

void curve2d(const Range &x, const Range &y, const QString &name = "");
void dot2d(const Range &x, const Range &y, const double &rad = 4, const QString &name = "");
void clear2d(int w = 800, int h = 600);
void show2d();
QPixmap drawnPixmap2d();
}

#endif // PLOT2D_H
