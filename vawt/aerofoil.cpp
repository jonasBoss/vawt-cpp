#include "aerofoil.hpp"
#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/detail/combine_cxx11.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>
#include <cassert>
#include <csv.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
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

void sort_by_re(vector<double>& alpha, vector<double>& re,vector<double>& cl,vector<double>& cd) {
    // Sort 're' in ascending order and rearrange other vectors accordingly.
    std::vector<std::pair<double, std::pair<double, std::pair<double, double>>>> combined;
    for (size_t i = 0; i < re.size(); ++i) {
        combined.push_back({re[i], {alpha[i], {cl[i], cd[i]}}});
    }
    std::stable_sort(combined.begin(), combined.end());
    for (size_t i = 0; i < re.size(); ++i) {
        re[i] = combined[i].first;
        alpha[i] = combined[i].second.first;
        cl[i] = combined[i].second.second.first;
        cd[i] = combined[i].second.second.second;
    }
}

void AerofoilBuilder::add_data(tuple<double, vector<double>, vector<double>, vector<double>> data){
    double re = get<0>(data);
    auto it = lower_bound(this->data.begin(), this->data.end(), re, [](const std::tuple<double, std::vector<double>, std::vector<double>, std::vector<double>>& tuple, double value) {
        return std::get<0>(tuple) < value;
    });
    this->data.insert(it, data);
}

AerofoilBuilder& AerofoilBuilder::load_data(string_view file, double re) {
    CSVFormat format;
    format.no_header().delimiter(',');
    unique_ptr<CSVReader> reader(new CSVReader(file, format));
    if (this->contains_re(re)) {
        throw "data is already loaded";
    }
    
    vector<double> alpha, cl, cd;

    for (CSVRow& row: *reader) {
        assert(row.size() == 3);
        auto iter =  row.begin();
        alpha.push_back(to_double(iter->get()) * TO_RAD);
        iter++;
        cl.push_back(to_double(iter->get()));
        iter++;
        cd.push_back(to_double(iter->get()));
    }
    this->add_data(tuple(re, alpha, cl, cd));
    return *this;
}

shared_ptr<Aerofoil> AerofoilBuilder::build() {
    std::list<std::tuple<double, std::vector<double>, std::vector<double>, std::vector<double>>> data;
    if (this->_update_aspect_ratio) {
        throw "not yet implemented";
    } else {
        data = this->data;
    }

    // duplicate highest and lowest values for extrapolation over re
    auto lowest = data.front();
    get<0>(lowest) = 0.0;
    data.push_front(lowest);
    auto highest = data.back();
    get<0>(highest) = numeric_limits<double>::max();
    data.push_back(highest);

    // collect everything into coniguous vectors for the interpolator
    vector<double> alpha, re, cl, cd;
    for (auto dataset: data){
        for (auto datapoint: boost::range::combine(get<1>(dataset), get<2>(dataset), get<3>(dataset))) {
            double _alpha, _cl, _cd;
            boost::tie(_alpha, _cl, _cd) = datapoint;
            re.push_back(get<0>(dataset));
            alpha.push_back(_alpha);
            cl.push_back(_cl);
            cd.push_back(_cd);
        }
    }
    return shared_ptr<Aerofoil>(new Aerofoil(alpha, re, cl, cd,this->_symmetric));
}