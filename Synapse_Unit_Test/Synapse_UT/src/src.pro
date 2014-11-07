include(../defaults.pri)
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
TEMPLATE = lib
TARGET = myapp
SOURCES += tridiag.cpp \
    Closed_form.cpp \
    Cra-Nic.cpp \
    Explicit.cpp \
    Implicit.cpp
HEADERS += tridiag.h \
    Closed_form.h \
    Cra-Nic.h \
    Explicit.h \
    Implicit.h

