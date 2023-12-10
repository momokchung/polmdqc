// Author: Moses KJ Chung
// Year:   2023

#include "readjson.h"
#include <iostream>
#include <fstream>
#include <sstream>

//////////////////////////////////////////////
//                                          //
//  readjson  --  read json to std::vector  //
//                                          //
//////////////////////////////////////////////

// "readjson" reads in the json file and outputs vector

std::vector<double> readVectorFromJson(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<double> vector;
    std::string line;

    while (std::getline(file, line)) {
        if (line.find_first_not_of("[\n]\n \t\n\v\f\r") == std::string::npos) {
            continue;
        }
        std::istringstream iss(line);
        double value;

        char peek;
        while (iss.get(peek) && peek != '[') {}

        while (iss >> value) {
            vector.push_back(value);
            peek = iss.peek();
            if (peek == ',' or peek == ']' or peek == '[') {
                iss.ignore();
            }
        }
    }

    return vector;
}

std::vector<std::vector<double>> readMatrixFromJson(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<std::vector<double>> matrix;
    std::string line;

    while (std::getline(file, line)) {
        if (line.find_first_not_of("[\n]\n \t\n\v\f\r") == std::string::npos) {
            continue;
        }
        std::vector<double> row;
        std::istringstream iss(line);
        double value;

        char peek;
        while (iss.get(peek) && peek != '[') {}

        while (iss >> value) {
            row.push_back(value);
            peek = iss.peek();
            if (peek == ',' or peek == ']' or peek == '[') {
                iss.ignore();
            }
        }
        matrix.push_back(row);
    }

    return matrix;
}
