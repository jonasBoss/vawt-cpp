#include "aerofoil.hpp"
#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/combine.hpp>
#include <cassert>
#include <csv.hpp>
#include <iomanip>
#include <iostream>
#include <locale>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

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

shared_ptr<Aerofoil> AerofoilBuilder::build() {
    vector<double> alpha;
    vector<double> re;
    vector<double> cl;
    vector<double> cd;


    if (this->_update_aspect_ratio) {
        throw "not yet implemented";
    } else {
        alpha = this->alpha;
        re = this->re;
        cl = this->cl;
        cd = this->cd;
    }

    cout << alpha.size()<< " " << re.size()<< " "  << cl.size()<< " "  << cd.size() << endl;
    assert(alpha.size() == re.size() && re.size()  == cl.size() && cl.size() == cd.size());
    // we know that re contains contiguous ranges of repeating values
    // find the rages for the minimum and maximum value
    auto min_max = minmax_element(re.begin(), re.end());
    auto min_range = equal_range(re.begin(), re.end(), *min_max.first);
    auto max_range = equal_range(re.begin(), re.end(), *min_max.second);

    // extrapolate the data for re = 0 to re = inf by
    // Appending values from each vector within their respective ranges.
    for (auto it = min_range.first; it != min_range.second; ++it) {
        size_t index = std::distance(re.begin(), it);
        alpha.push_back(alpha[index]);
        cl.push_back(cl[index]);
        cd.push_back(cd[index]);
    }

    for (auto it = max_range.first; it != max_range.second; ++it) {
        size_t index = std::distance(re.begin(), it);
        alpha.push_back(alpha[index]);
        cl.push_back(cl[index]);
        cd.push_back(cd[index]);
    }

    // Calculate the lengths of the ranges and append to 're'.
    size_t min_range_length = std::distance(min_range.first, min_range.second);
    size_t max_range_length = std::distance(max_range.first, max_range.second);

    for (size_t i = 0; i < min_range_length; ++i) {
        re.push_back(0.0); 
    }

    for (size_t i = 0; i < max_range_length; ++i) {
        re.push_back(std::numeric_limits<double>::infinity());
    }

    for (auto tup: boost::range::combine(alpha, re, cl, cd)) {
        double alpha, re, cl, cd;
        boost::tie(alpha, re, cl, cd) = tup;
        cout << setw(10) << alpha << setw(10) << re << setw(10) << cl << setw(10) << cd << endl;
    }

    return shared_ptr<Aerofoil>(new Aerofoil(alpha,re,cl,cd,this->_symmetric));
}