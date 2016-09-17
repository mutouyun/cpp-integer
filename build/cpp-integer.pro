TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++1y

INCLUDEPATH += ./../

HEADERS += ./../include/integer.hpp

SOURCES += ./../main.cpp
