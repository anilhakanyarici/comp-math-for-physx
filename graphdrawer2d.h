#ifndef GRAPHDRAWER2D_H
#define GRAPHDRAWER2D_H

#include <QSharedPointer>
#include <QVector>
#include <QString>
#include <QColor>
#include <QPaintDevice>
#include <QPainter>

#include "axis.h"

namespace mp {

class GraphDrawer2D
{
    struct pimpl;
    QSharedPointer<pimpl> _pimpl;

public:
    GraphDrawer2D();

    Axis &xAxis() const;
    Axis &yAxis() const;
    double zoom() const;
    double curveWidth() const;
    bool embrace() const;
    QFont &font();
    double axisWidth() const;
    QColor &axisColor() const;
    double gridWidth() const;
    QColor &gridColor() const;

    void setZoom(double zoom);
    void setCurveWidth(double w);
    void setEmbrace(bool s);
    void setFont(const QFont &font);
    void setFontSize(int s);
    void setAxisWidth(double w);
    void setAxisColor(const QColor &color);
    void setGridWidth(double w);
    void setGridColor(const QColor &color);
    void setGridIntensities(int x, int y);
    void setEnvelopeColor(const QColor &color);
    void setBackgroundColor(const QColor &color);

    void curve(const Range &x, const Range &y, const QString &name = "");
    void curve(const Range &x, const Range &y, const QString &name, const QColor &color);
    void dot(const Range &x, const Range &y, const double &rad, const QString &name = "");
    void dot(const Range &x, const Range &y, const double &rad, const QString &name, const QColor &color);


    void calculate(const QSizeF &size_p);
    void calculate(const QRectF &rect_p); //rect_p = Tüm çizimlerin yapılacağı Rect alanını belirler.
    void draw(QPaintDevice *device, QPainter::RenderHint hint = QPainter::HighQualityAntialiasing);

    GraphDrawer2D clone() const;

private:
    void drawAxis(QPainter &painter);
    void drawLegend(QPainter &painter);
    void addPlot(const Range &x, const Range &y, const QString &name, const QColor &color, int draw_type, double width = -1);
};
}

#endif // GRAPHDRAWER2D_H
