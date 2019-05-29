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

double ZBuffer::getZVal(int x, int y) const
{
    //std::cout << x << "   " << y << std::endl;
    return zVals.at(y*xW + x);
}

bool ZBuffer::setVal(int x, int y, double val) {
    if (val < zVals.at(y*xW + x))
    {
        //std::cout  << val << std::endl;
        zVals.at(y*xW + x) = val;
        return true;
    }
    else {
//        double temp = zVals.at(y*xW + x);
        return false;
    }
}

bool ZBuffer::compare(unsigned int x, unsigned int y, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, double zA, double zB) {

    double posInf = std::numeric_limits<double>::infinity();

    double p;

    if (x0 + y0 - x1 - y1 != 0) {
        p = (double)(x + y - x1 - y1) / (double)(x0 + y0 - x1 - y1);
    }
    else return this->setVal(x, y, 1/zA);

    double temp = p/zA + (1-p)/zB;

    return this->setVal(x,y,temp);
}

