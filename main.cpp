#include <iostream>
#include <ostream>
#include "vawt/aerofoil.hpp"
#include "external/csv-parser/single_include/csv.hpp"

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;


    csv::CSVFormat format;
    format.no_header().delimiter(',');
    csv::CSVReader reader("examples/NACA0018/NACA0018Re0040.data", format);
    for (csv::CSVRow& row: reader) {
        for (csv::CSVField& field: row ) {
            std::cout << field.get<>();
        }
        std::cout << std::endl;
    }

    std::cout << "bye" << std::endl;
    
    return 0;
}
