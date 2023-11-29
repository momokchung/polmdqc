////////////////////////////////////////////////
//                                            //
//  katomsqm.cpp  --  assign atom parameters  //
//                                            //
////////////////////////////////////////////////


#include "config.h"
#include "ioUtils.h"
#include "katomsqm.h"
#include "kbasis.h"
#include "kgbs.h"
#include "stringUtils.h"
#include "unitsqm.h"
#include <iostream>
#include <string>
#include <vector>

namespace atoms
{
// n               total number of atoms
// nElec           total number of electrons
// charge          total charge of system
// multiplicity    sping multiplicity
// com             center of mass
// reorient        reorient
// core            nuclear charge of atoms
// coordx          list of x coordinate of atoms
// coordy          list of y coordinate of atoms
// coordz          list of z coordinate of atoms
// name            name of system
// atom            list of atoms
// lengtUnit       length unit
// symmetry        point group symmetry

int n;
int nElec;
int charge;
int multiplicity;
bool com = false;
bool reorient = false;
std::vector<int> core;
std::vector<real> coordx;
std::vector<real> coordy;
std::vector<real> coordz;
std::string name;
std::vector<std::string> atom;
LengthUnit lengthUnit = LengthUnit::angstrom;
Symmetry symmetry = Symmetry::C1;

std::map<std::string, int> atomicNumber = {
    { "H" ,   1 },
    { "HE",   2 },
    { "LI",   3 },
    { "BE",   4 },
    { "B" ,   5 },
    { "C" ,   6 },
    { "N" ,   7 },
    { "O" ,   8 },
    { "F" ,   9 },
    { "NE",  10 },
    { "NA",  11 },
    { "MG",  12 },
    { "AL",  13 },
    { "SI",  14 },
    { "P" ,  15 },
    { "S" ,  16 },
    { "CL",  17 },
    { "AR",  18 },
    { "K" ,  19 },
    { "CA",  20 },
    { "SC",  21 },
    { "TI",  22 },
    { "V" ,  23 },
    { "CR",  24 },
    { "MN",  25 },
    { "FE",  26 },
    { "CO",  27 },
    { "NI",  28 },
    { "CU",  29 },
    { "ZN",  30 },
    { "GA",  31 },
    { "GE",  32 },
    { "AS",  33 },
    { "SE",  34 },
    { "BR",  35 },
    { "KR",  36 },
    { "RB",  37 },
    { "SR",  38 },
    { "Y" ,  39 },
    { "ZR",  40 },
    { "NB",  41 },
    { "MO",  42 },
    { "TC",  43 },
    { "RU",  44 },
    { "RH",  45 },
    { "PD",  46 },
    { "AG",  47 },
    { "CD",  48 },
    { "IN",  49 },
    { "SN",  50 },
    { "SB",  51 },
    { "TE",  52 },
    { "I" ,  53 },
    { "XE",  54 },
    { "CS",  55 },
    { "BA",  56 },
    { "LA",  57 },
    { "CE",  58 },
    { "PR",  59 },
    { "ND",  60 },
    { "PM",  61 },
    { "SM",  62 },
    { "EU",  63 },
    { "GD",  64 },
    { "TB",  65 },
    { "DY",  66 },
    { "HO",  67 },
    { "ER",  68 },
    { "TM",  69 },
    { "YB",  70 },
    { "LU",  71 },
    { "HF",  72 },
    { "TA",  73 },
    { "W" ,  74 },
    { "RE",  75 },
    { "OS",  76 },
    { "IR",  77 },
    { "PT",  78 },
    { "AU",  79 },
    { "HG",  80 },
    { "TL",  81 },
    { "PB",  82 },
    { "BI",  83 },
    { "PO",  84 },
    { "AT",  85 },
    { "RN",  86 },
    { "FR",  87 },
    { "RA",  88 },
    { "AC",  89 },
    { "TH",  90 },
    { "PA",  91 },
    { "U" ,  92 },
    { "NP",  93 },
    { "PU",  94 },
    { "AM",  95 },
    { "CM",  96 },
    { "BK",  97 },
    { "CF",  98 },
    { "ES",  99 },
    { "FM", 100 },
    { "MD", 101 },
    { "NO", 102 },
    { "LR", 103 },
    { "RF", 104 },
    { "DB", 105 },
    { "SG", 106 },
    { "BH", 107 },
    { "HS", 108 },
    { "MT", 109 },
    { "DS", 110 },
    { "RG", 111 },
};

std::map<std::string, real> atomicWeight = {
    { "H" , 1.00797 },
    { "HE", 4.0026 },
    { "LI", 6.941 },
    { "BE", 9.01218 },
    { "B" , 10.81 },
    { "C" , 12.011 },
    { "N" , 14.0067 },
    { "O" , 15.9994 },
    { "F" , 18.998403 },
    { "NE", 20.179 },
    { "NA", 22.98977 },
    { "MG", 24.305 },
    { "AL", 26.98154 },
    { "SI", 28.0855 },
    { "P" , 30.97376 },
    { "S" , 32.06 },
    { "CL", 35.453 },
    { "AR", 39.948 },
    { "K" , 39.0983 },
    { "CA", 40.08 },
    { "SC", 44.9559 },
    { "TI", 47.9 },
    { "V" , 50.9415 },
    { "CR", 51.996 },
    { "MN", 54.938 },
    { "FE", 55.847 },
    { "CO", 58.9332 },
    { "NI", 58.7 },
    { "CU", 63.546 },
    { "ZN", 65.38 },
    { "GA", 69.72 },
    { "GE", 72.59 },
    { "AS", 74.9216 },
    { "SE", 78.96 },
    { "BR", 79.904 },
    { "KR", 83.8 },
    { "RB", 85.4678 },
    { "SR", 87.62 },
    { "Y" , 88.9059 },
    { "ZR", 91.22 },
    { "NB", 92.9064 },
    { "MO", 95.94 },
    { "TC", 98.0 },
    { "RU", 101.07 },
    { "RH", 102.9055 },
    { "PD", 106.4 },
    { "AG", 107.868 },
    { "CD", 112.41 },
    { "IN", 114.82 },
    { "SN", 118.69 },
    { "SB", 121.75 },
    { "TE", 127.6 },
    { "I" , 126.9045 },
    { "XE", 131.3 },
    { "CS", 132.9054 },
    { "BA", 137.33 },
    { "LA", 138.9055 },
    { "CE", 140.12 },
    { "PR", 140.9077 },
    { "ND", 144.24 },
    { "PM", 145.0 },
    { "SM", 150.4 },
    { "EU", 151.96 },
    { "GD", 157.25 },
    { "TB", 158.9254 },
    { "DY", 162.5 },
    { "HO", 164.9304 },
    { "ER", 167.26 },
    { "TM", 168.9342 },
    { "YB", 173.04 },
    { "LU", 174.967 },
    { "HF", 178.49 },
    { "TA", 180.9479 },
    { "W" , 183.85 },
    { "RE", 186.207 },
    { "OS", 190.2 },
    { "IR", 192.22 },
    { "PT", 195.09 },
    { "AU", 196.9665 },
    { "HG", 200.59 },
    { "TL", 204.37 },
    { "PB", 207.2 },
    { "BI", 208.9804 },
    { "PO", 209.0 },
    { "AT", 210.0 },
    { "RN", 222.0 },
    { "FR", 223.0 },
    { "RA", 226.0254 },
    { "AC", 227.0278 },
    { "TH", 232.0381 },
    { "PA", 231.0359 },
    { "U" , 238.029 },
    { "NP", 237.0482 },
    { "PU", 242.0 },
    { "AM", 243.0 },
    { "CM", 247.0 },
    { "BK", 247.0 },
    { "CF", 251.0 },
    { "ES", 252.0 },
    { "FM", 257.0 },
    { "MD", 258.0 },
    { "NO", 250.0 },
    { "LR", 260.0 },
    { "RF", 261.0 },
    { "DB", 262.0 },
    { "SG", 263.0 },
    { "BH", 262.0 },
    { "HS", 255.0 },
    { "MT", 256.0 },
    { "DS", 269.0 },
    { "RG", 272.0 },
};

//////////////////////////////////////////
//                                      //
//  void readxyz  --  read in xyz file  //
//                                      //
//////////////////////////////////////////

void readxyz(std::string fileName)
{
    // initialize
    core.resize(0);
    coordx.resize(0);
    coordy.resize(0);
    coordz.resize(0);
    atom.resize(0);

    // check if file exists
    ioUtils::fileExists(fileName);

    // get number of lines
    int lineN = ioUtils::lineNumbers(fileName);

    // read line by line
    std::vector<std::string> lines(lineN);
    ioUtils::readlines(fileName, lines);

    // parse file
    for (int i = 0; i < lineN; ++i)
    {
        std::vector<std::string> words = stringUtils::split(lines[i]);

        // get machine memory
        if (words[0] == "MEMORY")
        {
            if (words.size() == 2)
            {
                config::memory = std::stoi(words[1]);
            }
            else if (words.size() > 2)
            {
                if (words[2] == "MB")
                    config::memory = std::stoi(words[1]);
                else if (words[2] == "GB")
                    config::memory = std::stoi(words[1]) * 1024;
            }
        }
        // get geometric configuration
        else if (words[0] == "MOLECULE")
        {
            if (words.size() != 3)
            {
                std::cerr << "Opening line of xyz file format must be \"molecule {name of moelcule} {\"\n"; 
                std::exit(1);
            }
            if (words[2] != "{")
            {
                std::cerr << "Opening line of xyz file format must be \"molecule {name of moelcule} {\"\n"; 
                std::exit(1);
            }
            name = words[1];
            i += 1;
            words = stringUtils::split(lines[i]);
            charge = std::stoi(words[0]);
            multiplicity = std::stoi(words[1]);
            i += 1;
            words = stringUtils::split(lines[i]);
            bool l = false;
            while (!l)
            {
                for (auto word: words)
                {
                    if (word.find("}") != std::string::npos)
                    {
                        l = true;
                        words.back() = words.back().substr(0, words.back().size() - 1);
                    }
                }
                if (words.size() == 4)
                {
                    atom.push_back(words[0]);
                    core.push_back(atomicNumber[words[0]]);
                    coordx.push_back(std::stod(words[1]));
                    coordy.push_back(std::stod(words[2]));
                    coordz.push_back(std::stod(words[3]));
                }
                else
                {
                    if (words[0] == "UNITS")
                    {
                        if (words[1] == "ANGSTROM")
                            lengthUnit = LengthUnit::angstrom;
                        else if (words[2] == "BOHR")
                            lengthUnit = LengthUnit::bohr;
                    }
                    else if (words[0] == "SYMMETRY")
                    {
                        if (words[1] == "C1")
                            symmetry = Symmetry::C1;
                        else if (words[1] == "C2")
                            symmetry = Symmetry::C2;
                        else if (words[1] == "Ci")
                            symmetry = Symmetry::Ci;
                        else if (words[1] == "Cs")
                            symmetry = Symmetry::Cs;
                        else if (words[1] == "D2")
                            symmetry = Symmetry::D2;
                        else if (words[1] == "C2h")
                            symmetry = Symmetry::C2h;
                        else if (words[1] == "C2v")
                            symmetry = Symmetry::C2v;
                        else if (words[1] == "D2h")
                            symmetry = Symmetry::D2h;
                    }
                    else if (words[0] == "COM")
                        com = true;
                    else if (words[0] == "NO_COM")
                        com = false;
                    else if (words[0] == "REORIENT")
                        reorient = true;
                    else if (words[0] == "NO_REORIENT")
                        reorient = false;
                }
                i += 1;
                words = stringUtils::split(lines[i]);
            }
            i -= 1;
        }
        // get basis set name
        else if (words[0] == "BASIS")
        {
            gbs::basisName = words[1];
        }
    }
    // get total number of atoms and electrons
    n = atom.size();
    nElec = 0;
    for (int i = 0; i < n; ++i)
    {
        nElec += core[i];
    }
    nElec -= charge;
    // convert to bohrs
    if (lengthUnit == LengthUnit::angstrom)
    {
        for (int i = 0; i < n; ++i)
        {
            coordx[i] = coordx[i] * unitsqm::bohr;
            coordy[i] = coordy[i] * unitsqm::bohr;
            coordz[i] = coordz[i] * unitsqm::bohr;
        }
    }
    // // uncomment below to check readxyz
    // std::cout << name << std::endl;
    // for (int i = 0; i < n; ++i)
    // {
    //     std::cout << atom[i] << std::endl;
    //     printf("core: %2i\n", core[i]);
    //     printf("coordx: %9.6f\n", coordx[i]);
    //     printf("coordy: %9.6f\n", coordy[i]);
    //     printf("coordz: %9.6f\n", coordz[i]);
    // }
    // if (lengthUnit == LengthUnit::angstrom)
    //     printf("angstrom\n");
    // else if (lengthUnit == LengthUnit::bohr)
    //     printf("bohr\n");
    // if (symmetry == Symmetry::C1)
    //     printf("C1\n");
}
}
