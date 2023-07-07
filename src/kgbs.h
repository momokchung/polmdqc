///////////////////
//               //
//  kgbs header  //
//               //
///////////////////


#pragma once
#include "init.h"
#include <map>
#include <string>
#include <vector>

namespace gbs
{
class AoBasis
{
public:
    //setters
    void setElementName(std::string name);
    void setShell(std::vector<std::string>& vector);
    void setShellN(int shellN);
    void setContraction(std::vector<int>& vector);
    void setScale(std::vector<real>& vector);
    void setPrimExp(std::vector<std::vector<real>>& vector);
    void setContractionCoeff(std::vector<std::vector<real>>& vector);

    // getters
    std::string getElementName();
    std::vector<std::string>& getShell();
    int getShellN();
    std::vector<int>& getContraction();
    std::vector<real>& getScale();
    std::vector<std::vector<real>>& getPrimExp();
    std::vector<std::vector<real>>& getContractionCoeff();

private:
    // Example of hydrogen 3-21.gbs
    // elementName = "H"
    // shell = {"S", "S"}
    // shellN = 2
    // contraction = {2, 1}
    // scale = {1.0, 1.0}
    // primExp = {{5.4471780, 0.8245470}, {0.1831920}}
    // contractionCoeff= {{0.1562850, 0.9046910}, {1.0000000}}

    std::string elementName;
    std::vector<std::string> shell;
    int shellN;
    std::vector<int> contraction;
    std::vector<real> scale;
    std::vector<std::vector<real>> primExp;
    std::vector<std::vector<real>> contractionCoeff;
};

enum class BasisType
{
    cartesian,
    spherical
};

extern std::string basisName;
extern BasisType basisType;
extern std::map<std::string, AoBasis> AoBasisMap;

void readgbs(std::string fileName);
}
