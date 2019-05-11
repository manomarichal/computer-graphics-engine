//
// Created by mano on 29.03.19.
//

#include "ZBuffer.h"

ZBuffer::ZBuffer(const int width, const int height) {

    double posInf = std::numeric_limits<double>::infinity();

    xW = width;
    yH = height;

    zVals = std::vector<double>(xW*yH, posInf);
}

double ZBuffer::getZVal(int x, int y) {
    return zVals.at(y*xW + x);
}

bool ZBuffer::setVal(int x, int y, double val) {
    if (val < zVals.at(y*xW + x)) {
        zVals.at(y*xW + x) = val;
        return true;
    }
    else return false;
}

bool ZBuffer::compare(unsigned int x, unsigned int y, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, double zA, double zB) {

    double p;

    if (x0 + y0 - x1 - y1 != 0) {
        p = (x + y - x1 - y1) / (x0 + y0 - x1 - y1);
    }
    else p = 0;

    double temp = p/zA + (1-p)/zB;

    return this->setVal(x,y,temp);
}

