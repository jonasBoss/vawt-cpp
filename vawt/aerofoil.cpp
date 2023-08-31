#include "aerofoil.hpp"
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <cassert>
#include <csv.hpp>
#include <iostream>
#include <locale>
#include <memory>
#include <string>

using namespace vawt;
using namespace csv;
using namespace std;

const double TO_RAD = boost::math::double_constants::pi / 180;

double to_double(string s) {
    boost::algorithm::trim(s);
    return boost::lexical_cast<double>(s);
}

AerofoilBuilder& AerofoilBuilder::load_data(string_view file, double re) {
    CSVFormat format;
    format.no_header().delimiter(',');
    unique_ptr<CSVReader> reader(new CSVReader(file, format));

    for (CSVRow& row: *reader) {
        assert(row.size() == 3);
        auto iter =  row.begin();
        this->alpha.push_back(to_double(iter->get()) * TO_RAD);
        iter++;
        this->cl.push_back(to_double(iter->get()));
        iter++;
        this->cd.push_back(to_double(iter->get()));
        this->re.push_back(re);
    }
    return *this;
}
