/*
 * easy_image.h
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
#ifndef EASY_IMAGE_INCLUDED
#define EASY_IMAGE_INCLUDED
#include <stdint.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "vector3d.h"
#include "ZBuffer.h"

struct Point2D {

	Point2D()= default;

	Point2D(double ix, double ij) {
		x = ix;
		y = ij;
	}

    double x,y;
};

struct Color {
	double red, green, blue;

	Color()
	{
		red = 0;
		blue = 0;
		green = 0;
	}
	void ini(const std::vector<double> &vec) {
		red = vec[0];
		green = vec[1];
		blue = vec[2];
	}

	void iniDouble(const double &r, const double &g, const double &b) {
		red = r;
		blue = b;
		green = g;
	}

	void iniColor(Color &color) {
		red = color.red;
		blue = color.blue;
		green = color.green;
	}
	void operator+=(const std::vector<double> &right)
	{
		red += right[0];
		blue += right[2];
		green += right[1];
	}
	void operator*=(const std::vector<double> &right)
	{
		red *= right[0];
		blue *= right[2];
		green *= right[1];
	}
	void operator*=(const double &d)
	{
		red *= d;
		blue *= d;
		green *= d;
	}
	std::vector<double> asVector()
	{
		std::vector<double> temp = {red, green, blue};
		return temp;
	}
};


inline int roundToInt(double d){
	return static_cast<int>(std::round(d));
}

/**
 * \brief The namespace of the EasyImage class
 */

class Light;

namespace img
{
	/**
	 * \brief This class represents the color of a pixel in an img::EasyImage object
	 */
	class Color
	{
		//a safety warning to all of you: Do *NOT* rearrange the 'color components' in this class
		//easyimage expects these three fields to be the *first* fields in the class AND expects
		//them to be in *this* order
		//if you alter the arrangement, the generated BMP files will contain garbage
		public:
			/**
			 * \brief The intensity of the blue color component
			 */
			uint8_t blue;

			/**
			 * \brief The intensity of the green color component
			 */
			uint8_t green;

			/**
			 * \brief The intensity of the red color component
			 */
			uint8_t red;

			/**
			 * \brief Default Constructor
			 */
			Color();

			/**
			 * \brief Constructs a Color with the given intensities
			 *
			 * \param r	The red color component
			 * \param g	The green color component
			 * \param b	The blue color component
			 *
			 */
			Color(uint8_t r, uint8_t g, uint8_t b);


			/**
			 * Destructor
			 */
			~Color();
	};

	/**
	 * \brief The exception that is thrown when an error occurs while trying to read an img::EasyImage from an input stream
	 */
	class UnsupportedFileTypeException: public std::exception
	{
		private:

			/**
			 * \brief Message explaining what went wrong
			 */
			std::string message;

		public:
			/**
			 * \brief Construct an exception with the given message
			 *
			 * \param msg	The message explaining what went wrong
			 *
			 */
			UnsupportedFileTypeException(std::string const& msg);

                        /**
                         * \brief Copy Constructor
                         *
                         * \param original	The exception to be copied into this object
                         */
			UnsupportedFileTypeException(const UnsupportedFileTypeException &original);


			/**
			 * \brief Destructor
			 */
			virtual ~UnsupportedFileTypeException() throw ();

			/**
			 * \brief Assignment operator
			 *
			 * \param original	The original exception to be assigned to this one
			 */
			UnsupportedFileTypeException& operator=(const UnsupportedFileTypeException &original);

			/**
                         * \brief Returns a description of the error hat occurred.
                         *
                         * \return A description of the error hat occurred.
                         */
			virtual const char *what() const throw ();
	};


	/**
	 * \brief This class implements a 'minor' image-library that supports basic operations such as setting and retrieving a pixel, and drawing a line.
	 */
	class EasyImage
	{
		public:
			/**
			 * \brief Default Constructor. Creates a zero-pixel image
			 */
			EasyImage();

			/**
			 * \brief Constructor: creates a new EasyImage of the specified width and height
			 *
			 * \param width		the width of the image
			 * \param height	the height of the image
			 * \param color		(optional) the background color of the image
			 */
			EasyImage(unsigned int width, unsigned int height, Color color = Color());

			/**
			 * \brief Copy Constructor
			 *
			 * \param img		the image to be copied
			 */
			EasyImage(EasyImage const& img);

			/**
			 * \brief Destructor
			 */
			virtual ~EasyImage();

			/**
			 * \brief Assignment operator. Allows an easyImage to be assigned to another easyImage
			 *
			 * \param img	The image to be assigned to this image
			 */
			EasyImage& operator=(EasyImage const& img);

			/**
			 * \brief Returns the width of the image
			 *
			 * \return the width of the image
			 */
			unsigned int get_width() const;

			/**
			 * \brief Returns the height of the image
			 * \return the height of the image
			 */
			unsigned int get_height() const;

			/**
			 * \brief Function operator. This operator returns a reference to a particular pixel of the image.
			 *
			 * \param x	the x coordinate of the pixel
			 * \param y	the y coordinate of the pixel
			 *
			 * These assertions apply:
			 * 	assert(x>=0 && x < getWidth())
			 * 	assert(y>=0 && y < getHeight())
			 */
			Color& operator()(unsigned int x, unsigned int y);

			/**
			 * \brief Function operator. This operator returns a const reference to a particular pixel of the image.
			 *
			 * \param x	the x coordinate of the pixel
			 * \param y	the y coordinate of the pixel
			 *
			 * These assertions apply:
			 * 	assert(x>=0 && x < getWidth())
			 * 	assert(y>=0 && y < getHeight())
			 */
			Color const& operator()(unsigned int x, unsigned int y) const;

			/**
			 * \brief Fills the image with a background of a specified color. Defaults to black
			 *
			 * \param color		The color to be assigned to each pixel
			 */
			void clear(Color color = Color());

			/**
			 * \brief Draws a line from pixel (x0,y0) to pixel (x1,y1) in the specified color
			 *
			 * \param x0	the x coordinate of the first pixel
			 * \param y0	the y coordinate of the first pixel
			 * \param x1	the x coordinate of the second pixel
			 * \param y1	the y coordinate of the second pixel
			 * \param color	the color of the line
			 *
			 * These assertions apply:
			 *	assert(x0 < getWidth())
			 * 	assert(y0 < getHeight())
			 * 	assert(x1 < getWidth())
			 * 	assert(y1 < getHeight())
			 */
			void draw_line(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, Color color);

			void draw_zbuf_line(
				ZBuffer &zBuf,
				unsigned int x0, unsigned int y0,
				double z0,
				unsigned int x1, unsigned int y1,
				double z1,
				Color color);

			static void calculateXlXr(Point2D &P, Point2D &Q, double &xl, double &xr, double y);

            void draw_zbuf_triangle(ZBuffer& zBuf,
                             Vector3D const& A,
                             Vector3D const& B,
                             Vector3D const& C,
                             Vector3D const &O,


                             double d,

                             double dx,
                             double dy,

							std::vector<double> ambientReflection,
							std::vector<double> diffuseReflection,
							std::vector<double> specularReflection,

							double reflectionCoeff,
							std::vector<Light>& lights,
							const Matrix &eyePointTrans, bool enableShadows);

        	static void draw_zbuf_triangle_colorless(ZBuffer& zBuf,
                                    Vector3D const& A,
                                    Vector3D const& B,
                                    Vector3D const& C,

                                    double d,

                                    double dx,
                                    double dy);

		void draw_zbuf_triangle_textures(ZBuffer& zBuf,
														 Vector3D const& A,
														 Vector3D const& B,
														 Vector3D const& C,
														 Vector3D const& O,

														 double d,

														 double dx,
														 double dy,

														 Matrix &eyepoint,
														 const img::EasyImage &texture);



		private:
			friend std::istream& operator>>(std::istream& in, EasyImage & image);
			/**
			 * \brief the width of the image
			 */
			unsigned int width;
			/**
			 * \brief the height of the image
			 */
			unsigned int height;
			/**
			 * \brief the vector containing all pixels
			 */
			std::vector<Color> bitmap;
	};

	/**
	 * \brief Writes an img::EasyImage to an output stream in the BMP file format
	 *
	 * \param out		the std::ostream to write the BMP file to.
	 * \param image		the img::EasyImage to be written to the output stream
	 *
	 * \return		a reference to the output stream the image was written to
	 */
	std::ostream& operator<<(std::ostream& out, EasyImage const& image);
	/**
	 * \brief Reads an img::EasyImage from an input stream.
	 *
	 * Please note: at this point only a limited subset of BMP-file format is supported.
	 * In order to correctly read a BMP file it must:
	 * 	- Be an uncompressed bitmap
	 *	- Only contain one plane
	 * 	- Use 24bits/pixel
	 * If the BMP file-format is not supported an UnsupportedFileTypeException is thrown
	 *
	 * \param in		the input stream to read the bitmap from
	 * \param image		the EasyImage object in which the bitmap must be stored
	 *
	 * \return		a reference to the input stream from which the bitmap was read
	 */
	std::istream& operator>>(std::istream& in, EasyImage& image);

}


class Light
{
public:
	Color ambientLight;
	Color diffuseLight;
	Color specularLight;

	bool difLight = false;
	bool amLight = false;
	bool specLight = false;
	bool tex = false;


	Vector3D ldVector;
	bool infinity;

	//shaduw
	ZBuffer shadowMask;
	Matrix eye;

	double d;
	double dx;
	double dy;

	img::EasyImage texture;
};

#endif /* EASY_IMAGE_INCLUDED */
