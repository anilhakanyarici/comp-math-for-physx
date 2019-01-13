#ifndef VIDEORECORDER_H
#define VIDEORECORDER_H

#include <opencv2/opencv.hpp>

#include <QSize>
#include <QImage>
#include <QString>

namespace mp {
class VideoRecorder
{
    cv::VideoWriter writer;
public:
    VideoRecorder(const QString &filename, int fps, const QSize &size);

    void write(const QImage &frame);
    void release();
};
}

#endif // VIDEORECORDER_H
