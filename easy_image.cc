/*
 * easy_image.cc
 * Copyright (C) 2011  Daniel van den Akker
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "easy_image.h"
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <iostream>

#define le32toh(x) (x)

namespace
{
	//structs borrowed from wikipedia's article on the BMP file format
	struct bmpfile_magic
	{
			uint8_t magic[2];
	};

	struct bmpfile_header
	{
			uint32_t file_size;
			uint16_t reserved_1;
			uint16_t reserved_2;
			uint32_t bmp_offset;
	};
	struct bmp_header
	{
			uint32_t header_size;
			int32_t width;
			int32_t height;
			uint16_t nplanes;
			uint16_t bits_per_pixel;
			uint32_t compress_type;
			uint32_t pixel_size;
			int32_t hres;
			int32_t vres;
			uint32_t ncolors;
			uint32_t nimpcolors;
	};
	//copy-pasted from lparser.cc to allow these classes to be used independently from each other
	class enable_exceptions
	{
		private:
			std::ios& ios;
			std::ios::iostate state;
		public:
			enable_exceptions(std::ios& an_ios, std::ios::iostate exceptions) :
				ios(an_ios)
			{
				state = ios.exceptions();
				ios.exceptions(exceptions);
			}
			~enable_exceptions()
			{
				ios.exceptions(state);
			}
	};
	//helper function to convert a number (char, int, ...) to little endian
	//regardless of the endiannes of the system
	//more efficient machine-dependent functions exist, but this one is more portable
	template<typename T> T to_little_endian(T value)
	{
		//yes, unions must be used with caution, but this is a case in which a union is needed
		union
		{
				T t;
				uint8_t bytes[sizeof(T)];
		} temp_storage;

		for (uint8_t i = 0; i < sizeof(T); i++)
		{
			temp_storage.bytes[i] = value & 0xFF;
			value >>= 8;
		}
		return temp_storage.t;
	}

	template<typename T> T from_little_endian(T value)
	{
		//yes, unions must be used with caution, but this is a case in which a union is needed
		union
		{
				T t;
				uint8_t bytes[sizeof(T)];
		} temp_storage;
		temp_storage.t = value;
		T retVal = 0;

		for (uint8_t i = 0; i < sizeof(T); i++)
		{
			retVal = (retVal << 8) | temp_storage.bytes[sizeof(T) - i - 1];
		}
		return retVal;
	}

}
img::Color::Color() :
	blue(0), green(0), red(0)
{
}
img::Color::Color(uint8_t r, uint8_t g, uint8_t b) :
	blue(b), green(g), red(r)
{
}
img::Color::~Color()
{
}

img::UnsupportedFileTypeException::UnsupportedFileTypeException(std::string const& msg) :
	message(msg)
{
}
img::UnsupportedFileTypeException::UnsupportedFileTypeException(const UnsupportedFileTypeException &original)
: std::exception(original)
, message(original.message)
{
}
img::UnsupportedFileTypeException::~UnsupportedFileTypeException() throw ()
{
}
img::UnsupportedFileTypeException& img::UnsupportedFileTypeException::operator=(UnsupportedFileTypeException const& original)
{
	this->message = original.message;
	return *this;
}
const char* img::UnsupportedFileTypeException::what() const throw ()
{
	return message.c_str();
}

img::EasyImage::EasyImage() :
	width(0), height(0), bitmap()
{
}

img::EasyImage::EasyImage(unsigned int _width, unsigned int _height, Color color) :
	width(_width), height(_height), bitmap(width * height, color)
{
}

img::EasyImage::EasyImage(EasyImage const& img) :
	width(img.width), height(img.height), bitmap(img.bitmap)
{
}

img::EasyImage::~EasyImage()
{
	bitmap.clear();
}

img::EasyImage& img::EasyImage::operator=(img::EasyImage const& img)
{
	width = img.width;
	height = img.height;
	bitmap.assign(img.bitmap.begin(),img.bitmap.end());
	return (*this);
}

unsigned int img::EasyImage::get_width() const
{
	return width;
}

unsigned int img::EasyImage::get_height() const
{
	return height;
}

void img::EasyImage::clear(Color color)
{
	for (std::vector<Color>::iterator i = bitmap.begin(); i != bitmap.end(); i++)
	{
		*i = color;
	}
}

img::Color& img::EasyImage::operator()(unsigned int x, unsigned int y)
{
	assert(x < this->width);
	assert(y < this->height);
	return bitmap.at(x * height + y);
}

img::Color const& img::EasyImage::operator()(unsigned int x, unsigned int y) const
{
	assert(x < this->width);
	assert(y < this->height);
	return bitmap.at(x * height + y);
}

void img::EasyImage::draw_line(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, Color color)
{
	assert(x0 < this->width && y0 < this->height);
	assert(x1 < this->width && y1 < this->height);
	if (x0 == x1)
	{
		//special case for x0 == x1
		for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
		{
			(*this)(x0, i) = color;
		}
	}
	else if (y0 == y1)
	{
		//special case for y0 == y1
		for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
		{
			(*this)(i, y0) = color;
		}
	}
	else
	{
		if (x0 > x1)
		{
			//flip points if x1>x0: we want x0 to have the lowest value
			std::swap(x0, x1);
			std::swap(y0, y1);
		}
		double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
		if (-1.0 <= m && m <= 1.0)
		{
			for (unsigned int i = 0; i <= (x1 - x0); i++)
			{
				(*this)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
			}
		}
		else if (m > 1.0)
		{
			for (unsigned int i = 0; i <= (y1 - y0); i++)
			{
				(*this)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
			}
		}
		else if (m < -1.0)
		{
			for (unsigned int i = 0; i <= (y0 - y1); i++)
			{
				(*this)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
			}
		}
	}
}

void img::EasyImage::draw_zbuf_line(ZBuffer &zBuf, unsigned int x0, unsigned int y0,
									double z0, unsigned int x1, unsigned int y1, double z1,
									Color color) {

	assert(x0 < this->width && y0 < this->height);
	assert(x1 < this->width && y1 < this->height);

	if (x0 == x1)
	{
		//special case for x0 == x1
		for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
		{
			if (zBuf.compare(x0,i,x0,y0,x1,x1,z0,z1)) {
				(*this)(x0, i) = color;
			}
		}
	}
	else if (y0 == y1)
	{
		//special case for y0 == y1
		for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
		{
			if (zBuf.compare(i,y0,x0,y0,x1,x1,z0,z1)) {
				(*this)(i, y0) = color;
			}
		}
	}
	else
	{
		if (x0 > x1)
		{
			//flip points if x1>x0: we want x0 to have the lowest value
			std::swap(x0, x1);
			std::swap(y0, y1);
		}
		double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
		if (-1.0 <= m && m <= 1.0)
		{
			for (unsigned int i = 0; i <= (x1 - x0); i++)
			{
				if (zBuf.compare(x0+i, (unsigned int) round(y0 + m * i) ,x0,y0,x1,x1,z0,z1)) {
					(*this)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
				}
			}
		}
		else if (m > 1.0)
		{
			for (unsigned int i = 0; i <= (y1 - y0); i++)
			{
				if (zBuf.compare((unsigned int) round(x0 + (i / m)), y0+i ,x0,y0,x1,x1,z0,z1)) {
					(*this)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
				}
			}
		}
		else if (m < -1.0)
		{
			for (unsigned int i = 0; i <= (y0 - y1); i++)
			{
				if (zBuf.compare((unsigned int) round(x0 - (i / m)),y0-i,x0,y0,x1,x1,z0,z1)) {
					(*this)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
				}
			}
		}
	}

}

void img::EasyImage::calculateXlXr(Point2D &P, Point2D &Q, double &xl, double &xr, double y) {

	if ( (y - P.y)*(y - Q.y) <= 0 and P.y != Q.y)
	{
		xl = Q.x + ( (P.x - Q.x) * ( (y - Q.y)/(P.y - Q.y) ) );
		xr = xl;
	}
}

void img::EasyImage::draw_zbuf_triangle(ZBuffer& zBuf,
										Vector3D const& A,
										Vector3D const& B,
										Vector3D const& C,

										double d,

										double dx,
										double dy,

										std::vector<double> ambientReflection,
										std::vector<double> diffuseReflection,
										std::vector<double> specularReflection,

										double reflectionCoeff,
										std::vector<Light>& lights,
										const Matrix &eyePointTrans, bool enableShadows) {




	double posInf = std::numeric_limits<double>::infinity();
	double negInf = -std::numeric_limits<double>::infinity();

	Point2D a((d * A.x/-A.z) + dx, (d * A.y/-A.z) + dy);
	Point2D b((d * B.x/-B.z) + dx, (d * B.y/-B.z) + dy);
	Point2D c((d * C.x/-C.z) + dx, (d * C.y/-C.z) + dy);

	int ymin = roundToInt(std::min(std::min(a.y, b.y), c.y)+0.5);
	int ymax = roundToInt(std::max(std::max(a.y, b.y), c.y)-0.5);

    //calculating of dzdx and dzdy
    Vector3D u = Vector3D::vector(B.x - A.x, B.y - A.y, B.z - A.z);
    Vector3D v = Vector3D::vector(C.x - A.x, C.y - A.y, C.z - A.z);

    Vector3D w = Vector3D::vector(u.cross_equals(v));

    // licht
	std::vector<double> color = {0, 0 ,0};
	Vector3D n = Vector3D::normalise(w);

    Point2D G;
    G.x = (a.x + b.x + c.x)/3;
    G.y = (a.y + b.y + c.y)/3;

    double zG = 1/(3*A.z) + 1/(3*B.z) + 1/(3*C.z);

	// rest
    double k = w.x * A.x + w.y * A.y + w.z * A.z;

    double dzdx = (w.x) / (-d*k);
    double dzdy = (w.y) / (-d*k);

	int xl; int xr;

	for (int y=ymin; y<=ymax; y++) {
        double xlab = posInf;
		double xlbc = posInf;
		double xlac = posInf;

		double xrab = negInf;
		double xrbc = negInf;
		double xrac = negInf;

		calculateXlXr(a,b,xlab, xrab,y);
		calculateXlXr(a,c,xlac, xrac,y);
		calculateXlXr(c,b,xlbc, xrbc,y);

		xl = roundToInt(std::min(std::min(xlab, xlac), xlbc) + 0.5);
		xr = roundToInt(std::max(std::max(xrab, xrac), xrbc) - 0.5);

		for (int x = xl;x<=xr;x++) {

			double zVal = 1.0001 * zG + (x - G.x)*dzdx + (y-G.y)*dzdy;

			// belichting en shaduw
			std::vector<double> colorTemp = color;

			for (const Light &light:lights)
			{
				if (light.amLight)
				{
					color[0] += light.ambientLight.red*ambientReflection[0];
					color[1] += light.ambientLight.green*ambientReflection[1];
					color[2] += light.ambientLight.blue*ambientReflection[2];
				}

				Vector3D point = Vector3D::point( (x - dx) / (d*(-zVal)), (y - dy) / (d*(-zVal)), 1/zVal);

                if (enableShadows)
                {
                    Vector3D L = point * eyePointTrans * light.eye;
                    Vector3D lA = Vector3D::point((light.d * L.x/-L.z)+light.dx, (light.d * L.y/-L.z)+light.dy, 0);

                    double ax = lA.x - std::floor(lA.x);
                    double ay = lA.y - std::floor(lA.y);

                    double az = light.shadowMask.getZVal(std::floor(lA.x), std::ceil(lA.y));
                    double bz = light.shadowMask.getZVal(std::ceil(lA.x), std::ceil(lA.y));
                    double cz = light.shadowMask.getZVal(std::floor(lA.x), std::floor(lA.y));
                    double dz = light.shadowMask.getZVal(std::ceil(lA.x), std::floor(lA.y));

                    double zeinvers = (1 - ax)*az + ax*bz;
                    double zfinvers = (1 - ax)*cz + ax*dz;

                    double zfinal = ay*zeinvers + (1-ay)*zfinvers;

                    if (std::abs(zfinal - (1/L.z)) > std::pow(10, -4)) continue;
                }

				if (light.difLight)
				{
					double cosAlpha;
					double cosBeita = -1;
					Vector3D l;

					if (light.infinity)
					{
						l = -light.ldVector;
						l.normalise();
						cosAlpha = Vector3D::dot(l ,n)*-1;
					}
					else {
                        l = Vector3D::vector(light.ldVector - point);
						l.normalise();
						cosAlpha = Vector3D::dot(l, n)*-1;
					}


					if (light.specLight)
					{
						Vector3D r = (2*cosAlpha*-1)*n - l;
						cosBeita = Vector3D::dot(r, Vector3D::normalise(Vector3D::point(0 , 0, 0) - point));
					}


                    if (cosBeita > 0)
                    {
                        double j = std::pow(cosBeita, reflectionCoeff)* cosBeita;
                        color[0] += light.specularLight.red * specularReflection[0]* j;
                        color[1] += light.specularLight.green * specularReflection[1] * j;
                        color[2] += light.specularLight.blue * specularReflection[2] * j;
                    }

                    if (cosAlpha > 0)
					{
						color[0] += light.diffuseLight.red * diffuseReflection[0] * cosAlpha;
						color[1] += light.diffuseLight.green * diffuseReflection[1] * cosAlpha;
						color[2] += light.diffuseLight.blue * diffuseReflection[2] * cosAlpha;
					}

				}
			}

			if (color[0] > 1) color[0] = 1;
			if (color[1] > 1) color[1] = 1;
			if (color[2] > 1) color[2] = 1;

            if (zBuf.setVal(x,y,zVal)) (*this)(x, y) = Color(color[0]*255, color[1]*255, color[2]*255);
            color = colorTemp;
		}
	}

}

void img::EasyImage::draw_zbuf_triangle_textures(ZBuffer& zBuf,
										Vector3D const& A,
										Vector3D const& B,
										Vector3D const& C,
										Vector3D const& O,

										double d,

										double dx,
										double dy,

										Matrix &eyepoint,
										const img::EasyImage &texture) {




	double posInf = std::numeric_limits<double>::infinity();
	double negInf = -std::numeric_limits<double>::infinity();

	Point2D a((d * A.x/-A.z) + dx, (d * A.y/-A.z) + dy);
	Point2D b((d * B.x/-B.z) + dx, (d * B.y/-B.z) + dy);
	Point2D c((d * C.x/-C.z) + dx, (d * C.y/-C.z) + dy);

	int ymin = roundToInt(std::min(std::min(a.y, b.y), c.y)+0.5);
	int ymax = roundToInt(std::max(std::max(a.y, b.y), c.y)-0.5);

	//calculating of dzdx and dzdy
	Vector3D u = Vector3D::vector(B.x - A.x, B.y - A.y, B.z - A.z);
	Vector3D v = Vector3D::vector(C.x - A.x, C.y - A.y, C.z - A.z);

	Vector3D w = Vector3D::vector(u.cross_equals(v));

	Point2D G;
	G.x = (a.x + b.x + c.x)/3;
	G.y = (a.y + b.y + c.y)/3;

	double zG = 1/(3*A.z) + 1/(3*B.z) + 1/(3*C.z);

	// rest
	double k = w.x * A.x + w.y * A.y + w.z * A.z;

	double dzdx = (w.x) / (-d*k);
	double dzdy = (w.y) / (-d*k);

	int xl; int xr;

	for (int y=ymin; y<=ymax; y++) {
		double xlab = posInf;
		double xlbc = posInf;
		double xlac = posInf;

		double xrab = negInf;
		double xrbc = negInf;
		double xrac = negInf;

		calculateXlXr(a,b,xlab, xrab,y);
		calculateXlXr(a,c,xlac, xrac,y);
		calculateXlXr(c,b,xlbc, xrbc,y);

		xl = roundToInt(std::min(std::min(xlab, xlac), xlbc) + 0.5);
		xr = roundToInt(std::max(std::max(xrab, xrac), xrbc) - 0.5);

		for (int x = xl;x<=xr;x++) {

            double zVal = 1.0001 * zG + (x - G.x)*dzdx + (y-G.y)*dzdy;

            Vector3D P = Vector3D::point( (x - dx) / (d*(-zVal)), (y - dy) / (d*(-zVal)), 1/zVal);

            Vector3D n = Vector3D::normalise(P - O);

            double u = std::asin(n.x) / M_PI + 0.5;
            double v = std::asin(n.y) / M_PI + 0.5;

            Color color = texture(roundToInt(1 + ( (texture.get_width() - 1) *u)), roundToInt(1 + ( (texture.get_height()-1) *v)));

            if (zBuf.setVal(x,y,zVal)) (*this)(x, y) = color;
		}
	}

}

void img::EasyImage::draw_zbuf_triangle_colorless(ZBuffer& zBuf,
                                        Vector3D const& A,
                                        Vector3D const& B,
                                        Vector3D const& C,

                                        double d,

                                        double dx,
                                        double dy) {


    double posInf = std::numeric_limits<double>::infinity();
    double negInf = -std::numeric_limits<double>::infinity();

    Point2D a((d * A.x/-A.z) + dx, (d * A.y/-A.z) + dy);
    Point2D b((d * B.x/-B.z) + dx, (d * B.y/-B.z) + dy);
    Point2D c((d * C.x/-C.z) + dx, (d * C.y/-C.z) + dy);

    int ymin = roundToInt(std::min(std::min(a.y, b.y), c.y)+0.5);
    int ymax = roundToInt(std::max(std::max(a.y, b.y), c.y)-0.5);

    //calculating of dzdx and dzdy
    Vector3D u = B - A;
    Vector3D v = C - A;
    Vector3D w = Vector3D::vector(u.cross_equals(v));

    double k = w.x * A.x + w.y * A.y + w.z * A.z;
    double dzdx = (w.x) / (-d*k);
    double dzdy = (w.y) / (-d*k);

    int xl; int xr;

    Point2D G;
    G.x = (a.x + b.x + c.x)/3;
    G.y = (a.y + b.y + c.y)/3;

    double zG = 1/(3*A.z) + 1/(3*B.z) + 1/(3*C.z);

    for (int y=ymin; y<=ymax; y++) {

        double xlab = posInf;
        double xlbc = posInf;
        double xlac = posInf;

        double xrab = negInf;
        double xrbc = negInf;
        double xrac = negInf;

        calculateXlXr(a,b,xlab, xrab,y);
        calculateXlXr(a,c,xlac, xrac,y);
        calculateXlXr(c,b,xlbc, xrbc,y);

        xl = roundToInt(std::min(std::min(xlab, xlac), xlbc) + 0.5);
        xr = roundToInt(std::max(std::max(xrab, xrac), xrbc) - 0.5);

        for (int x = xl;x<=xr;x++)
        {
            double zVal = zG + (x - G.x)*dzdx + (y-G.y)*dzdy;

            zBuf.setVal(x,y,zVal);
        }
    }
}
std::ostream& img::operator<<(std::ostream& out, EasyImage const& image)
{

	//temporaryily enable exceptions on output stream
	enable_exceptions(out, std::ios::badbit | std::ios::failbit);
	//declare some struct-vars we're going to need:
	bmpfile_magic magic;
	bmpfile_header file_header;
	bmp_header header;
	uint8_t padding[] =
	{ 0, 0, 0, 0 };
	//calculate the total size of the pixel data
	unsigned int line_width = image.get_width() * 3; //3 bytes per pixel
	unsigned int line_padding = 0;
	if (line_width % 4 != 0)
	{
		line_padding = 4 - (line_width % 4);
	}
	//lines must be aligned to a multiple of 4 bytes
	line_width += line_padding;
	unsigned int pixel_size = image.get_height() * line_width;

	//start filling the headers
	magic.magic[0] = 'B';
	magic.magic[1] = 'M';

	file_header.file_size = to_little_endian(pixel_size + sizeof(file_header) + sizeof(header) + sizeof(magic));
	file_header.bmp_offset = to_little_endian(sizeof(file_header) + sizeof(header) + sizeof(magic));
	file_header.reserved_1 = 0;
	file_header.reserved_2 = 0;
	header.header_size = to_little_endian(sizeof(header));
	header.width = to_little_endian(image.get_width());
	header.height = to_little_endian(image.get_height());
	header.nplanes = to_little_endian(1);
	header.bits_per_pixel = to_little_endian(24);//3bytes or 24 bits per pixel
	header.compress_type = 0; //no compression
	header.pixel_size = pixel_size;
	header.hres = to_little_endian(11811); //11811 pixels/meter or 300dpi
	header.vres = to_little_endian(11811); //11811 pixels/meter or 300dpi
	header.ncolors = 0; //no color palette
	header.nimpcolors = 0;//no important colors

	//okay that should be all the header stuff: let's write it to the stream
	out.write((char*) &magic, sizeof(magic));
	out.write((char*) &file_header, sizeof(file_header));
	out.write((char*) &header, sizeof(header));

	//okay let's write the pixels themselves:
	//they are arranged left->right, bottom->top, b,g,r
	for (unsigned int i = 0; i < image.get_height(); i++)
	{
		//loop over all lines
		for (unsigned int j = 0; j < image.get_width(); j++)
		{
			//loop over all pixels in a line
			//we cast &color to char*. since the color fields are ordered blue,green,red they should be written automatically
			//in the right order
			out.write((char*) &image(j, i), 3 * sizeof(uint8_t));
		}
		if (line_padding > 0)
			out.write((char*) padding, line_padding);
	}
	//okay we should be done
	return out;
}
std::istream& img::operator>>(std::istream& in, EasyImage & image)
{
	enable_exceptions(in, std::ios::badbit | std::ios::failbit);
	//declare some struct-vars we're going to need
	bmpfile_magic magic;
	bmpfile_header file_header;
	bmp_header header;
	//a temp buffer for reading the padding at the end of each line
	uint8_t padding[] =
	{ 0, 0, 0, 0 };

	//read the headers && do some sanity checks
	in.read((char*) &magic, sizeof(magic));
	if (magic.magic[0] != 'B' || magic.magic[1] != 'M')
		throw UnsupportedFileTypeException("Could not parse BMP File: invalid magic header");
	in.read((char*) &file_header, sizeof(file_header));
	in.read((char*) &header, sizeof(header));
	if (le32toh(header.pixel_size) + le32toh(file_header.bmp_offset) != le32toh(file_header.file_size))
		throw UnsupportedFileTypeException("Could not parse BMP File: file size mismatch");
	if (le32toh(header.header_size) != sizeof(header))
		throw UnsupportedFileTypeException("Could not parse BMP File: Unsupported BITMAPV5HEADER size");
	if (le32toh(header.compress_type) != 0)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only uncompressed BMP files can be parsed");
	if (le32toh(header.nplanes) != 1)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only one plane should exist in the BMP file");
	if (le32toh(header.bits_per_pixel) != 24)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only 24bit/pixel BMP's are supported");
	//if height<0 -> read top to bottom instead of bottom to top
	bool invertedLines = from_little_endian(header.height) < 0;
	image.height = std::abs(from_little_endian(header.height));
	image.width = std::abs(from_little_endian(header.width));
	unsigned int line_padding = from_little_endian(header.pixel_size) / image.height - (3 * image.width);
	//re-initialize the image bitmap
	image.bitmap.clear();
	image.bitmap.assign(image.height * image.width, Color());
	//okay let's read the pixels themselves:
	//they are arranged left->right., bottom->top if height>0, top->bottom if height<0, b,g,r
	for (unsigned int i = 0; i < image.get_height(); i++)
	{
		//loop over all lines
		for (unsigned int j = 0; j < image.get_width(); j++)
		{
			//loop over all pixels in a line
			//we cast &color to char*. since the color fields are ordered blue,green,red, the data read should be written in the right variables
			if (invertedLines)
			{
				//store top-to-bottom
				in.read((char*) &image(j, image.height - 1 - i), 3 * sizeof(uint8_t));
			}
			else
			{
				//store bottom-to-top
				in.read((char*) &image(j, i), 3 * sizeof(uint8_t));
			}
		}
		if (line_padding > 0)
		{
			in.read((char*) padding, line_padding);
		}
	}
	//okay we're done
	return in;
}
