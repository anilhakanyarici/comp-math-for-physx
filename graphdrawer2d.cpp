#include "graphdrawer2d.h"

#include <math.h>

using namespace mp;

namespace graph_drawer_cpp {
enum PlotType { PLOT_DOT = 1, PLOT_CURVE = 2 };

enum SyColor
{
    Sy_COLOR_WHITE = 0xffffff,
    Sy_COLOR_BLACK = 0x000000,
    Sy_COLOR_RED = 0xff0000,
    Sy_COLOR_LIME = 0x00ff00,
    Sy_COLOR_BLUE = 0x0000ff,
    Sy_COLOR_YELLOW = 0xb7b736,
    Sy_COLOR_CYAN = 0x00ffff,
    Sy_COLOR_MAGENTA = 0xff00ff,
    Sy_COLOR_ORANGE = 0xffa500,
    Sy_COLOR_SKY = 0xcceeff,
    Sy_COLOR_DARKBLUE = 0x0000a0,
    Sy_COLOR_LIGHTBLUE = 0xadd8e6,
    Sy_COLOR_PURPLE = 0x800080,
    Sy_COLOR_SILVER = 0xc0c0c0,
    Sy_COLOR_GRAY = 0x808080,
    Sy_COLOR_BROWN = 0xa52a2a,
    Sy_COLOR_MAROON = 0x810541,
    Sy_COLOR_OLIVE = 0x808000,
    Sy_COLOR_NAVYBLUE = 0x000080,
    Sy_COLOR_STEELBLUE = 0x4863a0,
    Sy_COLOR_GLACIALBLUE = 0x368bc1,
    Sy_COLOR_CRYSTALBLUE = 0x5cb3ff,
    Sy_COLOR_AQUAMARINE = 0x7fffd4,
    Sy_COLOR_TURQUOISE = 0x43c6db,
    Sy_COLOR_TEAL = 0x008080,
    Sy_COLOR_DARKFORESTGREEN = 0x254117,
    Sy_COLOR_FORESTGREEN = 0x4e9258,
    Sy_COLOR_SPRINGGREEN = 0x4aa02c,
    Sy_COLOR_JADE = 0x5efb6e,
    Sy_COLOR_DRAGONGREEN = 0x6afb92,
    Sy_COLOR_GOLDENROD = 0xedda74,
    Sy_COLOR_CREAM = 0xffffcc,
    Sy_COLOR_BEIGE = 0xf5f5dc,
    Sy_COLOR_MUSTARD = 0xffdb58,
    Sy_COLOR_CARAMEL = 0xc68e17,
    Sy_COLOR_CINNAMON = 0xc58917,
    Sy_COLOR_MOCHA = 0x493d26,
    Sy_COLOR_COFFEE = 0x6f4e37,
    Sy_COLOR_REDFOX = 0xc35817,
    Sy_COLOR_PUMPKIN = 0xf87217,
    Sy_COLOR_MANGO = 0xff8040,
    Sy_COLOR_CORAL = 0xff7f50,
    Sy_COLOR_VALENTINE = 0xe55451,
    Sy_COLOR_SCARLET = 0xff2400,
    Sy_COLOR_CRANBERRY = 0x9f000f,
    Sy_COLOR_CHESTNUT = 0x954535,
    Sy_COLOR_SIENNA = 0x8a4117,
    Sy_COLOR_PINK = 0xfaafbe,
    Sy_COLOR_VIOLET = 0xf6358a,
    Sy_COLOR_JASMINE = 0xa23bec,
    Sy_COLOR_CRIMSON = 0xe238ec,
    Sy_COLOR_INDIGO = 0x4b0082,
};

bool isValid(double d) { return d != INFINITY && d != -INFINITY && !std::isnan(d); }
double find_max(const double *array, int n) {
    if(n == 0)
        return 0;
    double max = array[0];
    int i = 0;
    for( ; i < n; i++) {
        if(isValid(array[i])){
            max = array[i];
            break;
        }
    }
    for( ; i < n; i++){
        if(array[i] > max && isValid(array[i]))
            max = array[i];
    }
    return max;
}
double find_min(const double *array, int n){
    if(n == 0)
        return 0;
    double min = array[0];
    int i = 0;
    for( ; i < n; i++) {
        if(isValid(array[i])){
            min = array[i];
            break;
        }
    }
    for( ; i < n; i++){
        if(array[i] < min && isValid(array[i]))
            min = array[i];
    }
    return min;
}
QVector<graph_drawer_cpp::SyColor> Sy_distinguishable_colors { graph_drawer_cpp::Sy_COLOR_ORANGE, graph_drawer_cpp::Sy_COLOR_LIME, graph_drawer_cpp::Sy_COLOR_RED, graph_drawer_cpp::Sy_COLOR_CRYSTALBLUE, graph_drawer_cpp::Sy_COLOR_YELLOW, graph_drawer_cpp::Sy_COLOR_MAGENTA, graph_drawer_cpp::Sy_COLOR_CYAN, graph_drawer_cpp::Sy_COLOR_WHITE, graph_drawer_cpp::Sy_COLOR_PURPLE, graph_drawer_cpp::Sy_COLOR_BLUE, graph_drawer_cpp::Sy_COLOR_PINK, graph_drawer_cpp::Sy_COLOR_BROWN, graph_drawer_cpp::Sy_COLOR_CORAL };
}


struct sy_plot {
    mp::Range y, x;
    QString name;
    QColor color;
    double width;
    int type;
};

struct GraphDrawer2D::pimpl {
    double y_min = 0;
    double y_max = 0;
    double x_min = 0;
    double x_max = 0;
    QMatrix R;

    double graph_ratio;
    double curve_width;

    bool embrace;

    QVector<sy_plot> pixel_space_plots;
    QVector<sy_plot> plots;
    Axis xAxis;
    Axis yAxis;

    double axisWidth;
    QColor axisColor;
    double gridWidth;
    QColor gridColor;
    QFont font;

    int grid_intensity_x;
    int grid_intensity_y;

    QColor envelope_color;
    QColor background_color;

    QRectF calculated_rect;
};

GraphDrawer2D GraphDrawer2D::clone() const
{
    GraphDrawer2D drawer;

    for(int i = 0; i < this->_pimpl->pixel_space_plots.size(); i++) drawer._pimpl->pixel_space_plots.append(this->_pimpl->pixel_space_plots[i]);
    for(int i = 0; i < this->_pimpl->plots.size(); i++) drawer._pimpl->plots.append(this->_pimpl->plots[i]);
    drawer._pimpl->xAxis = this->_pimpl->xAxis.clone();
    drawer._pimpl->yAxis = this->_pimpl->yAxis.clone();

    drawer._pimpl->y_min = this->_pimpl->y_min;
    drawer._pimpl->y_max = this->_pimpl->y_max;
    drawer._pimpl->x_min = this->_pimpl->x_min;
    drawer._pimpl->x_max = this->_pimpl->x_max;
    drawer._pimpl->R = this->_pimpl->R;
    drawer._pimpl->graph_ratio = this->_pimpl->graph_ratio;
    drawer._pimpl->curve_width = this->_pimpl->curve_width;
    drawer._pimpl->embrace = this->_pimpl->embrace;
    drawer._pimpl->axisWidth = this->_pimpl->axisWidth;
    drawer._pimpl->axisColor = this->_pimpl->axisColor;
    drawer._pimpl->gridWidth = this->_pimpl->gridWidth;
    drawer._pimpl->gridColor = this->_pimpl->gridColor;
    drawer._pimpl->font = this->_pimpl->font;
    drawer._pimpl->grid_intensity_x = this->_pimpl->grid_intensity_x;
    drawer._pimpl->grid_intensity_y = this->_pimpl->grid_intensity_y;
    drawer._pimpl->envelope_color = this->_pimpl->envelope_color;
    drawer._pimpl->background_color = this->_pimpl->background_color;
    drawer._pimpl->calculated_rect = this->_pimpl->calculated_rect;

    return drawer;
}

GraphDrawer2D::GraphDrawer2D()
{
    this->_pimpl = QSharedPointer<pimpl>(new pimpl());
    this->_pimpl->xAxis.setName("x");
    this->_pimpl->yAxis.setName("y");
    this->_pimpl->graph_ratio = 0.8;
    this->_pimpl->curve_width = 2.0;
    this->_pimpl->embrace = false;
    this->_pimpl->axisColor = QColor(255, 255, 255);
    //this->_pimpl->axisColor = QColor(0, 0, 0);
    this->_pimpl->axisWidth = 1;
    this->_pimpl->font.setPointSize(12);
    this->_pimpl->font.setPixelSize(12);
    this->_pimpl->grid_intensity_x = 10;
    this->_pimpl->grid_intensity_y = 10;
    this->_pimpl->gridWidth = 1;
    this->_pimpl->gridColor = QColor(150, 150, 150);
    //this->_pimpl->background_color = QColor(255, 255, 255);
    this->_pimpl->background_color = QColor(0, 0, 0);
    this->_pimpl->envelope_color = QColor(0, 128, 128);
}

Axis &GraphDrawer2D::xAxis() const
{
    return this->_pimpl->xAxis;
}

Axis &GraphDrawer2D::yAxis() const
{
    return this->_pimpl->yAxis;
}

double GraphDrawer2D::curveWidth() const
{
    return this->_pimpl->curve_width;
}

bool GraphDrawer2D::embrace() const
{
    return this->_pimpl->embrace;
}

QFont &GraphDrawer2D::font()
{
    return this->_pimpl->font;
}

double GraphDrawer2D::axisWidth() const
{
    return this->_pimpl->axisWidth;
}

QColor &GraphDrawer2D::axisColor() const
{
    return this->_pimpl->axisColor;
}

double GraphDrawer2D::gridWidth() const
{
    return this->_pimpl->gridWidth;
}

QColor &GraphDrawer2D::gridColor() const
{
    return this->_pimpl->gridColor;
}

void GraphDrawer2D::setCurveWidth(double w)
{
    this->_pimpl->curve_width = w;
}

void GraphDrawer2D::setEmbrace(bool s)
{
    if(this->_pimpl->plots.size() > 0)
        throw std::runtime_error("Embrace can able to set when there is no added plot.");
    this->_pimpl->embrace = s;
}

void GraphDrawer2D::setFont(const QFont &font)
{
    this->_pimpl->font = font;
}

void GraphDrawer2D::setFontSize(int s)
{
    this->font().setPixelSize(s);
}

void GraphDrawer2D::setAxisWidth(double w)
{
    this->_pimpl->axisWidth = w;
}

void GraphDrawer2D::setAxisColor(const QColor &color)
{
    this->_pimpl->axisColor = color;
}

void GraphDrawer2D::setGridWidth(double w)
{
    this->_pimpl->gridWidth = w;
}

void GraphDrawer2D::setGridColor(const QColor &color)
{
    this->_pimpl->gridColor = color;
}

void GraphDrawer2D::setGridIntensities(int x, int y)
{
    if(x < 0) x = 0;
    if(y < 0) y = 0;
    this->_pimpl->grid_intensity_x = x;
    this->_pimpl->grid_intensity_y = y;
}

void GraphDrawer2D::setEnvelopeColor(const QColor &color)
{
    this->_pimpl->envelope_color = color;
}

void GraphDrawer2D::setBackgroundColor(const QColor &color)
{
    this->_pimpl->background_color = color;
}

void GraphDrawer2D::curve(const mp::Range &x, const mp::Range &y, const QString &name)
{
    QColor color = QColor((QRgb)graph_drawer_cpp::Sy_distinguishable_colors[this->_pimpl->plots.size() % graph_drawer_cpp::Sy_distinguishable_colors.size()]);
    this->curve(x, y, name, color);
}
void GraphDrawer2D::curve(const mp::Range &x, const mp::Range &y, const QString &name, const QColor &color)
{
    this->addPlot(x, y, name, color, graph_drawer_cpp::PLOT_CURVE);
}

void GraphDrawer2D::dot(const Range &x, const Range &y, const double &rad, const QString &name)
{
    QColor color = QColor((QRgb)graph_drawer_cpp::Sy_distinguishable_colors[this->_pimpl->plots.size() % graph_drawer_cpp::Sy_distinguishable_colors.size()]);
    this->dot(x, y, rad, name, color);
}

void GraphDrawer2D::dot(const Range &x, const Range &y, const double &rad, const QString &name, const QColor &color)
{
    this->addPlot(x, y, name, color, graph_drawer_cpp::PLOT_DOT, rad);
}

void GraphDrawer2D::calculate(const QSizeF &size_p)
{
    this->calculate(QRectF(QPointF(0, 0), size_p));
}
void GraphDrawer2D::calculate(const QRectF &rect_p)
{
    Axis x_ax = this->xAxis();
    Axis y_ax = this->yAxis();

    this->_pimpl->calculated_rect = QRectF(rect_p.x(), rect_p.y(), rect_p.width(), rect_p.height());

    QFontMetrics font_metric = QFontMetrics(this->font());
    QString x_label_str = QString("%1 (x10^%2)").arg(x_ax.name(), QString::number(x_ax.rank()));
    int xl_w = font_metric.width(x_label_str);
    QString y_label_str = QString("%1 (x10^%2)").arg(y_ax.name(), QString::number(y_ax.rank()));
    font_metric.width(y_label_str);
    int yl_h = font_metric.height();

    double graph_ratio = this->_pimpl->graph_ratio;
    double xp = rect_p.x() - xl_w / 8; double wp = rect_p.width();
    double yp = rect_p.y() + yl_h / 4; double hp = rect_p.height() - 50;
    double wd = graph_ratio * wp - xl_w / 2;
    double hd = graph_ratio * hp - yl_h / 2;

    double x_min, x_max, y_min, y_max;

    if(this->_pimpl->embrace){
        x_min = x_ax.exponent() * x_ax.ticks()[0];
        x_max = x_ax.exponent() * x_ax.ticks()[x_ax.ticks().size() - 1];
        y_min = y_ax.exponent() * y_ax.ticks()[0];
        y_max = y_ax.exponent() * y_ax.ticks()[y_ax.ticks().size() - 1];
    } else {
        x_min = this->_pimpl->x_min;
        x_max = this->_pimpl->x_max;
        y_min = this->_pimpl->y_min;
        y_max = this->_pimpl->y_max;
    }

    double dx = x_max - x_min;
    double dy = y_max - y_min;

    double xr = wd / dx;
    double yr = -hd / dy;

    double xt = (wp - wd) / 2 + xp - x_min * xr;
    double yt = hp - (hp - hd) / 2 + yp - y_min * yr;

    this->_pimpl->R = QMatrix(xr, 0, 0, yr, xt, yt);

    this->_pimpl->pixel_space_plots = QVector<sy_plot>();

    QVector<sy_plot> plots = this->_pimpl->plots;
    for(int i = 0; i < plots.size(); i++){
        sy_plot plot = plots[i];
        int count = plot.x.size();
        double *xs = plot.x.data();
        double *ys = plot.y.data();
        sy_plot calculated_plot; //Converted plot from data space to pixel space.
        calculated_plot.width = plot.width;
        calculated_plot.color = plot.color;
        calculated_plot.name = plot.name;
        calculated_plot.type = plot.type;
        calculated_plot.x = mp::Range(count);
        calculated_plot.y = mp::Range(count);
        double *cx = calculated_plot.x.data();
        double *cy = calculated_plot.y.data();

        for(int j = 0; j < count; j++){
            cx[j] = xs[j] * xr + xt;
            cy[j] = ys[j] * yr + yt;
        }
        this->_pimpl->pixel_space_plots.append(calculated_plot);
    }
}

void GraphDrawer2D::draw(QPaintDevice *device, QPainter::RenderHint hint)
{
    QPainter painter(device);
    painter.setRenderHint(hint);

    painter.fillRect(QRectF(0, 0, painter.device()->width(), painter.device()->height()), QBrush(this->_pimpl->envelope_color));

    if(this->_pimpl->plots.size() > 0)
        this->drawAxis(painter);

    QVector<sy_plot> plots = this->_pimpl->pixel_space_plots;
    QPen pen = painter.pen();
    for(int i = 0; i < plots.size(); i++) {
        sy_plot plot = plots[i];
        pen.setColor(plot.color);
        pen.setWidthF(plot.width);
        painter.setPen(pen);
        int count = plot.x.size();
        double *x = plot.x.data();
        double *y = plot.y.data();

        if(plot.type == graph_drawer_cpp::PLOT_CURVE){
            for(int j = 1; j < count; j++) {
                if(!graph_drawer_cpp::isValid(y[j - 1]) || !graph_drawer_cpp::isValid(y[j]))
                    continue;
                QLineF line(x[j - 1], y[j - 1], x[j], y[j]);
                painter.drawLine(line);
            }
        }
        else if(plot.type == graph_drawer_cpp::PLOT_DOT){
            for(int j = 0; j < count; j++) {
                if(!graph_drawer_cpp::isValid(y[j]))
                    continue;
                painter.drawEllipse(QPointF(x[j], y[j]), plot.width / 2, plot.width / 2);
            }
        }
    }
    this->drawLegend(painter);
    painter.end();
}

void GraphDrawer2D::drawLegend(QPainter &painter)
{
    QRectF rect_p = this->_pimpl->calculated_rect;
    rect_p = QRect(0, rect_p.height() - 50, rect_p.width(), 25);

    QPen pen = painter.pen();
    painter.fillRect(rect_p, QBrush(QColor(64, 64, 64)));

    double L = rect_p.width();
    double H = rect_p.height();
    int n = this->_pimpl->plots.size();
    double legend_line_pixel = 40;
    QFontMetrics metrics = QFontMetrics(this->font());
    for(int i = 0; i < n; i++) {
        sy_plot plot = this->_pimpl->pixel_space_plots[i];
        double w = metrics.width(plot.name);
        w += legend_line_pixel + 5;
        double h = metrics.height();

        double p_x = rect_p.x() + ((i + 1) * L / (n + 1) - w / 2);
        double p_y = rect_p.y() + (H - h) / 2;

        pen.setWidthF(plot.width);
        pen.setColor(plot.color);
        painter.setPen(pen);
        painter.drawLine(QLineF(p_x, p_y + h / 2, p_x + legend_line_pixel, p_y + h / 2));
        pen.setColor(this->_pimpl->axisColor);
        painter.setPen(pen);
        painter.setFont(this->_pimpl->font);
        painter.drawText(p_x + legend_line_pixel + 5, p_y + h / 2 + metrics.height() / 4, plot.name);
    }
}

void GraphDrawer2D::addPlot(const Range &x, const Range &y, const QString &name, const QColor &color, int draw_type, double width)
{
    sy_plot plot;
    plot.x = x;
    plot.y = y;
    plot.name = name;
    plot.type = draw_type;
    if(name.isEmpty())
        plot.name = "(noname)";
    plot.color = color;
    plot.width = width <= 0 ? this->_pimpl->curve_width : width;

    double y_min = graph_drawer_cpp::find_min(y.data(), y.size());
    double y_max = graph_drawer_cpp::find_max(y.data(), y.size());
    double x_min = graph_drawer_cpp::find_min(x.data(), x.size());
    double x_max = graph_drawer_cpp::find_max(x.data(), x.size());

    if(this->_pimpl->plots.size() == 0){
        this->_pimpl->x_max = x_max;
        this->_pimpl->x_min = x_min;
        this->_pimpl->y_max = y_max;
        this->_pimpl->y_min = y_min;
    }{
        this->_pimpl->x_max = fmax(this->_pimpl->x_max, x_max);
        this->_pimpl->x_min = fmin(this->_pimpl->x_min, x_min);
        this->_pimpl->y_max = fmax(this->_pimpl->y_max, y_max);
        this->_pimpl->y_min = fmin(this->_pimpl->y_min, y_min);
    }

    this->_pimpl->plots.append(plot);
    this->_pimpl->xAxis.setInterval(this->_pimpl->x_min, this->_pimpl->x_max, this->_pimpl->embrace);
    this->_pimpl->yAxis.setInterval(this->_pimpl->y_min, this->_pimpl->y_max, this->_pimpl->embrace);
}

void GraphDrawer2D::drawAxis(QPainter &painter)
{
    QMatrix R = this->_pimpl->R;
    double xr = R.m11();
    double yr = R.m22();
    double xt = R.dx();
    double yt = R.dy();

    double x_min, x_max, y_min, y_max;

    Axis x_ax = this->xAxis();
    Axis y_ax = this->yAxis();

    if(this->_pimpl->embrace){
        x_min = x_ax.exponent() * x_ax.ticks()[0];
        x_max = x_ax.exponent() * x_ax.ticks()[x_ax.ticks().size() - 1];
        y_min = y_ax.exponent() * y_ax.ticks()[0];
        y_max = y_ax.exponent() * y_ax.ticks()[y_ax.ticks().size() - 1];
    } else {
        x_min = this->_pimpl->x_min;
        x_max = this->_pimpl->x_max;
        y_min = this->_pimpl->y_min;
        y_max = this->_pimpl->y_max;
    }

    //Pixel Space
    double x_Ax_start = x_min * xr + xt; //Pixel Space
    double y_Ax_start = y_min * yr + yt; //Pixel Space
    double x_Ax_end = x_max * xr + xt; //Pixel Space
    double y_Ax_end = y_max * yr + yt; //Pixel Space

    painter.fillRect(QRectF(x_Ax_start, y_Ax_end, x_Ax_end - x_Ax_start, y_Ax_start - y_Ax_end), QBrush(this->_pimpl->background_color));

    QPen pen = painter.pen();
    pen.setColor(this->_pimpl->axisColor);
    pen.setWidthF(this->_pimpl->axisWidth);
    painter.setPen(pen);

    painter.setFont(this->_pimpl->font);
    QFontMetrics font_metric = painter.fontMetrics();

    QLineF x_axis_line(x_Ax_start, y_Ax_start, x_Ax_end, y_Ax_start);
    QLineF y_axis_line(x_Ax_start, y_Ax_start, x_Ax_start, y_Ax_end);
    painter.drawLine(x_axis_line);
    painter.drawLine(y_axis_line);

    for(int i = 0 ; i < x_ax.ticks().size(); i++) { //x üzerindeki tickler
        double x_tick = x_ax.exponent() * x_ax.ticks()[i] * xr + xt; //Pixel Space
        QLineF tick = QLineF(x_tick, y_Ax_start - 5, x_tick, y_Ax_start + 5);
        painter.drawLine(tick);
        QString tick_label = QString::number(x_ax.ticks()[i], 'f', x_ax.precision());
        int tl_w = font_metric.width(tick_label);
        int tl_h = font_metric.height();
        painter.drawText(QPointF(x_tick - tl_w / 2, y_Ax_start + 5 + tl_h), tick_label);
    }

    for(int i = 0 ; i < y_ax.ticks().size(); i++) { //y üzerindeki tickler
        double y_tick = y_ax.exponent() * y_ax.ticks()[i] * yr + yt; //Pixel Space
        QLineF tick = QLineF(x_Ax_start - 5, y_tick, x_Ax_start + 5, y_tick);
        painter.drawLine(tick);
        QString tick_label = QString::number(y_ax.ticks()[i], 'f', y_ax.precision());
        int tl_w = font_metric.width(tick_label);
        int tl_h = font_metric.height();
        painter.drawText(QPointF(x_Ax_start - 5 - tl_w, y_tick + tl_h / 2), tick_label);
    }

    QString x_label_str = QString("%1 (x10^%2)").arg(x_ax.name(), QString::number(x_ax.rank()));
    font_metric.width(x_label_str);
    int xl_h = font_metric.height();
    painter.drawText(QPointF(x_Ax_end + 5, y_Ax_start + xl_h / 2), x_label_str);

    QString y_label_str = QString("%1 (x10^%2)").arg(y_ax.name(), QString::number(y_ax.rank()));
    int yl_w = font_metric.width(y_label_str);
    int yl_h = font_metric.height();
    painter.drawText(QPointF(x_Ax_start - yl_w / 2, y_Ax_end - yl_h / 2), y_label_str);


    pen.setWidthF(this->_pimpl->gridWidth);
    if(this->_pimpl->grid_intensity_x > 0){
        pen.setColor(QColor(this->_pimpl->gridColor.red(), this->_pimpl->gridColor.green(), this->_pimpl->gridColor.blue(), 255 / 3));
        painter.setPen(pen);
        double x_tick = x_ax.exponent() * x_ax.ticks()[0] * xr + xt; //Pixel Space
        double grid_interval_x = ((x_ax.exponent() * x_ax.ticks()[1] * xr + xt) - x_tick) / this->_pimpl->grid_intensity_x ; //Pixel Space
        double curr_grid_point = x_tick;
        while((float)curr_grid_point <= (float)x_Ax_end){ //x_tick, 1. tick değerini alır. Döngü, +x yönünde ilerleyerek y'ye paralel gridleri çizer.
            QLineF grid_line(curr_grid_point, y_Ax_start, curr_grid_point, y_Ax_end);
            painter.drawLine(grid_line);
            curr_grid_point += grid_interval_x;
        }
        curr_grid_point = x_tick;
        while((float)curr_grid_point > (float)x_Ax_start){ //x_tick, 1. tick değerini alır. Döngü, -x yönünde ilerleyerek y'ye paralel gridleri çizer.
            QLineF grid_line(curr_grid_point, y_Ax_start, curr_grid_point, y_Ax_end);
            painter.drawLine(grid_line);
            curr_grid_point -= grid_interval_x;
        }
        pen.setColor(this->_pimpl->gridColor);
        painter.setPen(pen);
        for(int i = 0 ; i < x_ax.ticks().size(); i++) { //ticklere karşılık gelen grid kenarları.
            x_tick = x_ax.exponent() * x_ax.ticks()[i] * xr + xt; //Pixel Space
            QLineF tick = QLineF(x_tick, y_Ax_start, x_tick, y_Ax_end);
            painter.drawLine(tick);
        }
    }
    if(this->_pimpl->grid_intensity_y > 0){
        pen.setColor(QColor(this->_pimpl->gridColor.red(), this->_pimpl->gridColor.green(), this->_pimpl->gridColor.blue(), 255 / 3));
        painter.setPen(pen);
        double y_tick = y_ax.exponent() * y_ax.ticks()[0] * yr + yt; //Pixel Space
        double grid_interval_y = ((y_ax.exponent() * y_ax.ticks()[1] * yr + yt) - y_tick) / this->_pimpl->grid_intensity_y; //Pixel Space
        double curr_grid_point = y_tick;
        while((float)curr_grid_point >= (float)y_Ax_end){ //y_tick, 1. tick değerini alır. Döngü, +y yönünde ilerleyerek x'e paralel gridleri çizer.
            QLineF grid_line(x_Ax_start, curr_grid_point, x_Ax_end, curr_grid_point);
            painter.drawLine(grid_line);
            curr_grid_point += grid_interval_y;
        }
        curr_grid_point = y_tick;
        while((float)curr_grid_point < (float)y_Ax_start){ //y_tick, 1. tick değerini alır. Döngü, -y yönünde ilerleyerek x'e paralel gridleri çizer.
            QLineF grid_line(x_Ax_start, curr_grid_point, x_Ax_end, curr_grid_point);
            painter.drawLine(grid_line);
            curr_grid_point -= grid_interval_y;
        }
        pen.setColor(this->_pimpl->gridColor);
        painter.setPen(pen);
        for(int i = 0 ; i < y_ax.ticks().size(); i++) { //ticklere karşılık gelen grid kenarları.
            y_tick = y_ax.exponent() * y_ax.ticks()[i] * yr + yt; //Pixel Space
            QLineF tick = QLineF(x_Ax_start, y_tick, x_Ax_end, y_tick);
            painter.drawLine(tick);
        }
    }
}
