QT += gui widgets

CONFIG += c++11 console
CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

LIBS += -L/home/anil/Desktop/opencv-3.4.5/install/lib -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_photo -lopencv_imgproc -lopencv_core
INCLUDEPATH += /home/anil/Desktop/opencv-3.4.5/install/include

SOURCES += \
        main.cpp \
    axis.cpp \
    func.cpp \
    graphdrawer2d.cpp \
    plot2d.cpp \
    range.cpp \
    runtest.cpp \
    complex.cpp \
    fourier.cpp \
    polynomial.cpp \
    videorecorder.cpp \
    Algebra/matrix3x3.cpp \
    Algebra/vector3.cpp \
    Algebra/quaternion.cpp \
    Algebra/matrix.cpp \
    Algebra/vector.cpp \
    Algebra/vector4.cpp \
    Algebra/matrix4x4.cpp \
    biginteger.cpp \
    drbg.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    axis.h \
    func.h \
    graphdrawer2d.h \
    plot2d.h \
    range.h \
    runtest.h \
    complex.h \
    fourier.h \
    polynomial.h \
    videorecorder.h \
    Algebra/matrix3x3.h \
    Algebra/vector3.h \
    Algebra/quaternion.h \
    Algebra/matrix.h \
    Algebra/vector.h \
    Algebra/vector4.h \
    Algebra/matrix4x4.h \
    biginteger.h \
    drbg.h
