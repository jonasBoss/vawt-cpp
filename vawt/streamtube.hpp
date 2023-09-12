#pragma once

namespace vawt {

class StreamTube {
  private:
    double a_0;
    double theta;
    double beta;

  public:
    /**
     * @brief Construct a new Stream Tube object
     * 
     * @param theta - stramtube position in turbine (radians)
     * @param beta - foil pitch angle (radians)
     * @param a_0 - upstream induction factor when `theta < PI` this is probably 0
     */
    StreamTube(double theta, double beta, double a_0) {
        this->a_0 = a_0;
        this->beta = beta;
        this->theta = theta;
    }

};

} // namespace vawt