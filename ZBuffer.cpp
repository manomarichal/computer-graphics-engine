//
// Created by mano on 29.03.19.
//

#include "ZBuffer.h"

ZBuffer::ZBuffer(const int width, const int height) {
    double posInf = std::numeric_limits<double>::infinity();
    double negInf = -std::numeric_limits<double>::infinity();

    for (int x=0;x<width;x++) {
        std::vector<double> temp;
        for (int y = 0; y < height; y++) {
            temp.emplace_back(posInf);
        }
    }
}

double ZBuffer::getZVal(int x, int y) {
    return zVals[x][y];
}

void ZBuffer::setVal(int x, int y, double val) {
    zVals[x][y] = val;
}
