#include "videorecorder.h"

using namespace mp;

VideoRecorder::VideoRecorder(const QString &filename, int fps, const QSize &size)
{
    this->writer = cv::VideoWriter(filename.toUtf8().data(), CV_FOURCC('M', 'J', 'P', 'G'), fps, cv::Size(size.width(), size.height()), false);
}

void VideoRecorder::write(const QImage &frame)
{
    unsigned int height = frame.height();
    unsigned int width = frame.width();

    cv::Mat3b dest(height, width);
    for (unsigned int y = 0; y < height; ++y) {
        cv::Vec3b *destrow = dest[y];
        for (unsigned int x = 0; x < width; ++x) {
            QRgb pxl = frame.pixel(x, y);
            destrow[x] = cv::Vec3b(::qBlue(pxl), ::qGreen(pxl), ::qRed(pxl));
        }
    }
    this->writer.write(dest);
}

void VideoRecorder::release()
{
    if(this->writer.isOpened())
        this->writer.release();
}
