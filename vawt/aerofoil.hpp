#ifndef AEROFOIL_HPP
#define AEROFOIL_HPP

#include <Interpolators/_2D/BilinearInterpolator.hpp>
#include <list>
#include <tuple>
#include <vector>

namespace vawt {
class AerofoilBuilder;
class Aerofoil;
using DataRow = std::tuple<double, std::vector<double>, std::vector<double>,
                           std::vector<double>>;
using DataSet = std::list<DataRow>;

/**
 * @brief Aerofoil coefficients of lift and drag
 *
 */
class ClCd {
    friend Aerofoil;

  private:
    double _cl;
    double _cd;
    ClCd(double cl, double cd) : _cl(cl), _cd(cd) {}

  public:
    /**
     * @brief Coefficient of Lift
     *
     * @return double
     */
    double cl() { return this->_cl; };

    /**
     * @brief Coefficient of Drag
     *
     * @return double
     */
    double cd() { return this->_cd; };

    /**
     * @brief Convert coefficients to normal and tangential turbine coordinates
     *
     * @param alpha - the foil angle of attac in radians
     * @param beta - the pitch agle from turbine tangent to wing chord in
     * radians
     * @return std::pair<double, double> - (normal coefficinent, tangential
     * coefficient)
     */
    std::pair<double, double> to_tangential(double alpha, double beta);

    /**
     * @brief Convert coefficients to global xy direction
     *
     * @param alpha - the foil angle of attac in radians
     * @param beta - the pitch agle from turbine tangent to wing chord in
     * radians
     * @param theta - the position angle at the turbine in radians
     * @return std::pair<double, double> - (x-coefficient, y-coefficient)
     */
    std::pair<double, double> to_global(double alpha, double beta,
                                        double theta);
};

class Aerofoil {
    friend AerofoilBuilder;

  private:
    bool symmetric;
    _2D::BilinearInterpolator<double> cl;
    _2D::BilinearInterpolator<double> cd;
    Aerofoil(std::vector<double> alpha, std::vector<double> re,
             std::vector<double> cl, std::vector<double> cd, bool symmetric) {
        this->cl.setData(re, alpha, cl);
        this->cd.setData(re, alpha, cd);
        this->symmetric = symmetric;
    }

  public:
    /**
     * @brief lift and drag coefficients
     *
     * @param alpha
     * @param re
     * @return ClCd
     */
    ClCd cl_cd(double alpha, double re) {
        double cl, cd;
        if (this->symmetric) {
            double sgn = (alpha >= 0) ? 1 : -1;
            cl = this->cl(re, abs(alpha)) * sgn;
            cd = this->cd(re, abs(alpha));
        } else {
            cl = this->cl(re, alpha);
            cd = this->cd(re, alpha);
        }
        return ClCd(cl, cd);
    }
};

class AerofoilBuilder {
  private:
    DataSet data;
    bool _symmetric = false;
    bool _update_aspect_ratio = false;
    double aspect_ratio = std::numeric_limits<double>::infinity();

    /**
     * @brief is data for the reynodlsnumber available?
     *
     * @param re
     * @return true
     * @return false
     */
    bool contains_re(double re) {
        return find_if(this->data.begin(), this->data.end(),
                       [re](const DataRow& tuple) {
                           return std::get<0>(tuple) == re;
                       }) != this->data.end();
    }

    /**
     * @brief adds a datarow to `this.data`
     *
     * The datarow is inserted such that `this.data` is sorted by `re`
     *
     * @param data
     */
    void add_data(DataRow data);

    /**
     * @brief create a new dataset with transformed datapoints for aspect ratio
     * correction
     *
     * @return DataSet
     */
    DataSet transformed_set();

    /**
     * @brief transform the datapoints in the row for aspect ratio correction
     *
     * @param row
     * @return DataRow
     */
    void transform_row(DataRow& row);

  public:
    /**
     * @brief load aerofoil data for a given reynolds number from a file
     *
     * The file is expected to be in csv format delimited by `,` without a
     * header. It should contain 3 columns: alpha (in degrees), cl, cd
     *
     * @param file - the file path
     * @param re - the reynoldsnumber for the data
     * @return AerofoilBuilder&
     */
    AerofoilBuilder& load_data(std::string_view file, double re);

    /**
     * @brief is the aerofoil profile symmetric
     *
     * @param symmetric
     * @return AerofoilBuilder&
     */
    AerofoilBuilder& symmetric(bool yes) {
        this->_symmetric = yes;
        return *this;
    }

    /**
     * @brief Set the aspect ratio of the areofoil
     *
     * When the data does not reflect this aspect ratio, but instead
     * is profile data for an infinte aspect ratio set [`update_aspect_ratio`]
     *
     * @param ar
     * @return AerofoilBuilder&
     */
    AerofoilBuilder& set_aspect_ratio(double ar) {
        this->aspect_ratio = ar;
        return *this;
    }

    /**
     * @brief Assume the provided data to be for a infinite aspect ratio and
     * change it.
     *
     * When the Aerofoil is built, the data will be
     * Updtated with the Lanchester-Prandtl model below the stalling angle
     * and with the Viterna-Corrigan above the stall angle
     *
     * @param yes
     * @return AerofoilBuilder&
     */
    AerofoilBuilder& update_aspect_ratio(bool yes) {
        this->_update_aspect_ratio = yes;
        return *this;
    }

    /**
     * @brief build the Aerofoil
     *
     * @return std::shared_ptr<Aerofoil>
     */
    std::shared_ptr<Aerofoil> build();
};

} // namespace vawt

#endif // AEROFOIL_HPP