//
// Created by mano on 29.03.19.
//

#include "ZBuffer.h"

ZBuffer::ZBuffer(const int width, const int height) {

    double posInf = std::numeric_limits<double>::infinity();

    for (int x=0;x<width;x++) {
        std::vector<double> temp;
        for (int y = 0; y < height; y++) {
            temp.emplace_back(posInf);
        }
        zVals.emplace_back(temp);
    }
}

double ZBuffer::getZVal(int x, int y) {
    return zVals[x][y];
}

void ZBuffer::setVal(int x, int y, double val) {
    zVals[x][y] = val;
}

bool ZBuffer::compare(unsigned int x, unsigned int y, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, double zA, double zB) {

    double p;

    if (x0 + y0 - x1 - y1 != 0) {
        p = (x + y - x1 - y1) / (x0 + y0 - x1 - y1);
    }
    else p = 0;

    double temp = p/zA + (1-p)/zB;

    if (temp < zVals[x][y]) {
        zVals[x][y] = temp;
        return true;
    }
    else return false;
}

