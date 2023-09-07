#include "aerofoil.hpp"
#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/detail/combine_cxx11.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>
#include <cassert>
#include <cmath>
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
const double PI = boost::math::double_constants::pi;

/**
 * @brief cast a string to double
 * 
 * @param s 
 * @return double 
 */
double to_double(string s) {
    boost::algorithm::trim(s);
    return boost::lexical_cast<double>(s);
}

struct DataPoint {
    double& alpha;
    double& cl;
    double& cd;
};

/**
 * @brief update alpha and cd for aspect ratio correction
 * 
 * @param p - the datapoint
 * @param ar - aspect ratio
 */
void lanchester_prandtl(DataPoint p, double ar) {
    p.cd = p.cd + p.cl * p.cl / (PI * ar);
    p.alpha = p.alpha + p.cl / (PI *ar);
}

/**
 * @brief update cl and cd values inplace for aspect ratio correction
 * 
 * @param p - the datapoint
 * @param stall - the data at the stall point
 * @param ar - aspect ratio
 */
void viterna_corrigan(DataPoint p, DataPoint stall, double ar) {
    double cd_max;
    if (ar > 50.0){
        cd_max = 2.01;
    } else {
        cd_max = 1.1 + 0.018 * ar;
    }
    double kd = (stall.cd - cd_max * pow(sin(stall.alpha), 2)) / cos(stall.alpha);
    double kl = (stall.cl - cd_max * sin(stall.alpha) * cos(stall.alpha)) * sin(stall.alpha) / pow(cos(stall.alpha), 2);
    p.cl = cd_max / 2.0 * sin(2.0 * p.alpha) + kl * pow(cos(p.alpha), 2) / sin(p.alpha);
    p.cd = cd_max * pow(sin(p.alpha), 2) + kd * cos(p.alpha);
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