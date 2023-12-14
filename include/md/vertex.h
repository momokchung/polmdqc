/*********************************************************************************
 *	The Vertex class
 *********************************************************************************/
#pragma once
#include <bitset>
#include <cmath> 
#include <vector>

/*********************************************************************************
  Class that characterizes each vertex of the Delaunay/Alpha complex
 *********************************************************************************/

class Vertex {
	public:
		double Radius;
		double Coordinates[3];
		double Weight;
		double CoefS, CoefV, CoefM, CoefG;
		double gamma;

		std::bitset<8> info;
		bool status;

		Vertex() {
		}

		Vertex(double x, double y, double z, double radius, double aspS,
			double aspV, double aspM, double aspG);
		~Vertex();

	private:
		double truncate_real(double x, int ndigit);
};
