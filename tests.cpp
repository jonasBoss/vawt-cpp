#include <boost/algorithm/string/trim.hpp>
#include <cassert>
#include <cmath>
#include <memory>
#include <vawt.hpp>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <ostream>
#include <vector>
#include <csv.hpp>

using namespace vawt;
using namespace csv;
using namespace std;

const double TO_RAD = boost::math::double_constants::pi / 180;

double to_double(string s) {
    boost::algorithm::trim(s);
    return boost::lexical_cast<double>(s);
}

bool rel_eq(double a, double b, double rel, double epsilon){
    double abs_diff = fabs(a - b);
    double rel_diff = (abs_diff / max(fabs(a), fabs(b)));
    if (rel_diff <= rel){
        return true;
    }
    return abs_diff <= epsilon;
}

class MatlabSolution{
public:
    vector<double> theta, a, w, alpha, re, CtubeThru, Ctan, Cnorm; 
    MatlabSolution(){
        CSVFormat format;
        format.delimiter('\t');
        unique_ptr<CSVReader> reader(new CSVReader("examples/matlab_NACA0018_tsr-3.25.txt", format));

        for (CSVRow& row: *reader) {
            assert(row.size() == 8);
            auto iter = row.begin();

            theta.push_back(to_double(iter->get()) * TO_RAD);
            iter++;
            a.push_back(to_double(iter->get()));
            iter++;
            w.push_back(to_double(iter->get()));
            iter++;
            alpha.push_back(to_double(iter->get()) * TO_RAD);
            iter++;
            re.push_back(to_double(iter->get()));
            iter++;
            CtubeThru.push_back(to_double(iter->get()));
            iter++;
            Ctan.push_back(to_double(iter->get()));
            iter++;
            Cnorm.push_back(to_double(iter->get()));
        }
    }

    uint n_streamtubes(){
        return this->theta.size();
    }
};


int main(int argc, char** argv) {
    std::cout << "Preparing Test Enviroment" << std::endl;

    vawt::AerofoilBuilder* builder = new vawt::AerofoilBuilder;
    auto aerofoil =
        builder->load_data("examples/NACA0018/NACA0018Re0080.data", 80000.0)
            .load_data("examples/NACA0018/NACA0018Re0040.data", 40000.0)
            .load_data("examples/NACA0018/NACA0018Re0160.data", 160000.0)
            .set_aspect_ratio(12.8)
            .update_aspect_ratio(true)
            .symmetric(true)
            .build();
    
    std::cout << "Loading Matlab solution" << std::endl;
    MatlabSolution* matlab = new MatlabSolution();

    std::cout << "Solving Turbine" << std::endl;
    auto testresult = VAWTSolver(aerofoil)
        .re(31'000.0)
        .solidity(0.3525)
        .n_streamtubes(matlab->n_streamtubes())
        .tsr(3.25)
        .solve(0.0);

    std::cout << "Checking results" << std::endl;
    for (int i=0; i< matlab->n_streamtubes(); i++){
        double theta = matlab->theta[i];
        double a = matlab->a[i];
        assert(rel_eq(a, testresult.a(theta), 0.01, testresult.epsilon() * 2));
    }
    std::cout << "Ok!" << std::endl;
    return 0;
}