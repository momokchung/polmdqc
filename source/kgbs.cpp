///////////////////////////////////////////
//                                       //
//  kgbs.cpp  --  assign gbs parameters  //
//                                       //
///////////////////////////////////////////


#include "init.h"
#include "ioUtils.h"
#include "kgbs.h"
#include "stringUtils.h"
#include <algorithm>
#include <iostream>

namespace gbs
{
// basisName     name of basis set
// basisType     cartesian or spherical
// AoBasisMap    map from atom to AoBasis

std::string basisName;
BasisType basisType;
std::map<std::string, AoBasis> AoBasisMap;

// AoBasis class methods
//setters
void AoBasis::setElementName(std::string name)
{
    elementName = name;
}
void AoBasis::setShell(std::vector<std::string>& vector)
{
    shell = vector;
}
void AoBasis::setShellN(int integer)
{
    shellN = integer;
}
void AoBasis::setContraction(std::vector<int>& vector)
{
    contraction = vector;
}
void AoBasis::setScale(std::vector<real>& vector)
{
    scale = vector;
}
void AoBasis::setPrimExp(std::vector< std::vector<real> >& vector)
{
    primExp = vector;
}
void AoBasis::setContractionCoeff(std::vector< std::vector<real> >& vector)
{
    contractionCoeff = vector;
}

// getters
std::string AoBasis::getElementName()
{
    return elementName;
}
std::vector<std::string>& AoBasis::getShell()
{
    return shell;
}
int AoBasis::getShellN()
{
    return shellN;
}
std::vector<int>& AoBasis::getContraction()
{
    return contraction;
}
std::vector<real>& AoBasis::getScale()
{
    return scale;
}
std::vector< std::vector<real> >& AoBasis::getPrimExp()
{
    return primExp;
}
std::vector< std::vector<real> >& AoBasis::getContractionCoeff()
{
    return contractionCoeff;
}

/////////////////////////////////////////
//                                     //
//  void readgbs  --  gbs file reader  //
//                                     //
/////////////////////////////////////////

void readgbs(std::string basisName)
{
    std::string fileName = init::cwd;
    fileName.append("/basis/");
    std::transform(basisName.begin(), basisName.end(), basisName.begin(), tolower);
    fileName.append(basisName);
    fileName.append(".gbs");

    // check if file exists
    ioUtils::fileExists(fileName);

    // get number of lines
    int lineN = ioUtils::lineNumbers(fileName);

    // read line by line
    std::vector< std::string > lines(lineN);
    ioUtils::readlines(fileName, lines);

    // parse file
    for (int i = 0; i < lineN; ++i)
    {
        std::vector< std::string > words = stringUtils::split(lines[i]);
        if (words[0] == "CARTESIAN")
            basisType = BasisType::cartesian;
        else if (words[0] == "SPHERICAL")
            basisType = BasisType::spherical;
        else if (words[0] == "****")
        {
            i += 1;
            if (i == lineN)
                break;
            words = stringUtils::split(lines[i]);
            std::string atom = words[0];
            AoBasis aobasis;
            aobasis.setElementName(atom);
            std::vector<std::string> shellVec;
            std::vector<int> contractionVec;
            std::vector<real> scaleVec;
            std::vector< std::vector<real> > primExpVec;
            std::vector< std::vector<real> > contractionCoeffVec;
            i += 1;
            words = stringUtils::split(lines[i]);
            std::string shell;
            int contractions;
            while (words[0] != "****")
            {
                std::string word = words[0];
                std::replace_if(word.begin(), word.end(), [](char ch) {return (ch == 'D');}, 'E');
                real realNum;
                bool isDouble = stringUtils::checkIsDouble(words[0], realNum);
                if (!isDouble)
                {
                    shell = words[0];
                    contractions = std::stoi(words[1]);
                    real scale = std::stod(words[2]);
                    if (shell == "SP")
                    {
                        shellVec.push_back("S");
                        shellVec.push_back("P");
                        contractionVec.push_back(contractions);
                        contractionVec.push_back(contractions);
                        scaleVec.push_back(scale);
                        scaleVec.push_back(scale);
                    }
                    else
                    {
                        shellVec.push_back(shell);
                        contractionVec.push_back(contractions);
                        scaleVec.push_back(scale);
                    }
                    i += 1;
                    words = stringUtils::split(lines[i]);
                }
                else
                {
                    std::vector<real> primExpSubVec1;
                    std::vector<real> contractionCoeffSubVec1;
                    std::vector<real> primExpSubVec2;
                    std::vector<real> contractionCoeffSubVec2;
                    for (int j = 0; j < contractions; ++j)
                    {
                        words = stringUtils::split(lines[i]);
                        std::replace_if(words[0].begin(), words[0].end(), [](char ch) {return (ch == 'D');}, 'E');
                        std::replace_if(words[1].begin(), words[1].end(), [](char ch) {return (ch == 'D');}, 'E');
                        primExpSubVec1.push_back(std::stod(words[0]));
                        contractionCoeffSubVec1.push_back(std::stod(words[1]));
                        if (shell == "SP")
                        {
                            std::replace_if(words[2].begin(), words[2].end(), [](char ch) {return (ch == 'D');}, 'E');
                            primExpSubVec2.push_back(std::stod(words[0]));
                            contractionCoeffSubVec2.push_back(std::stod(words[2]));
                        }
                        i += 1;
                    }
                    words = stringUtils::split(lines[i]);
                    primExpVec.push_back(primExpSubVec1);
                    contractionCoeffVec.push_back(contractionCoeffSubVec1);
                    if (shell == "SP")
                    {
                        primExpVec.push_back(primExpSubVec2);
                        contractionCoeffVec.push_back(contractionCoeffSubVec2);
                    }
                }
            }
            aobasis.setShell(shellVec);
            aobasis.setShellN(shellVec.size());
            aobasis.setContraction(contractionVec);
            aobasis.setScale(scaleVec);
            aobasis.setPrimExp(primExpVec);
            aobasis.setContractionCoeff(contractionCoeffVec);
            AoBasisMap.insert( {atom, aobasis} );
            i -= 1;
        }
    }

    // // uncomment below to check readgbs
    // AoBasis& H = gbs::AoBasisMap["O"];
    // std::string elementName = H.getElementName();
    // std::vector<std::string>& shell = H.getShell();
    // std::vector<int>& contraction = H.getContraction();
    // std::vector<real>& scale = H.getScale();
    // std::vector< std::vector<real> >& primExp = H.getPrimExp();
    // std::vector< std::vector<real> >& contractionCoeff = H.getContractionCoeff();

    // std::cout << "element Name: " << elementName << std::endl;
    // for (int i = 0; i < shell.size(); ++i)
    // {
    //     std::cout << i << std::endl;
    //     std::cout << "shell " << shell[i] << std::endl;
    //     std::cout << "contraction " << contraction[i] << std::endl;
    //     std::cout << "scale " << scale[i] << std::endl;
    //     std::cout << "primExp" << std::endl;
    //     for (int j = 0; j < contraction[i]; ++j)
    //     {
    //         printf("%15.6f", primExp[i][j]);
    //     }
    //     std::cout << std::endl;
    //     std::cout << "contractionCoeff" << std::endl;
    //     for (int j = 0; j < contraction[i]; ++j)
    //     {
    //         printf("%15.6f", contractionCoeff[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
}

}
