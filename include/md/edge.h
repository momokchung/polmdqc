/*********************************************************************************
 *	The Edge class
 *********************************************************************************/

#pragma once

/*********************************************************************************
  Edge class
 *********************************************************************************/

class Edge {
public:
    int Vertices[2];
    double gamma;
    double sigma;
    double Length, Surf, Vol; 
    double CoefM1, CoefM2, CoefG1, CoefG2;
    double dsurf, dvol, dmean, dgauss;

    Edge() {
    }

    Edge(int i, int j);

    ~Edge();
};
