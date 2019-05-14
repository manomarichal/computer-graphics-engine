//
// Created by mano on 29.03.19.
//

#ifndef CGENGINE_ZBUFFER_H
#define CGENGINE_ZBUFFER_H

#include <vector>
#include <limits>
#include <iostream>

class ZBuffer
{
        std::vector<double> zVals;
    public:
        int xW;
        int yH;
        ZBuffer(const int width, const int height);

        ZBuffer()=default;

        double getZVal(int x, int y) const;

        bool setVal(int x, int y, double val);

        bool compare(unsigned int x, unsigned int y, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, double zA, double zB);

};


#endif //CGENGINE_ZBUFFER_H
