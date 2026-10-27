#ifndef INTERPOLATE_HPP
#define INTERPOLATE_HPP

#include <vector>

namespace WaveformUtilities {
  
  std::vector<double> Interpolate(const std::vector<double>& X1, const std::vector<double>& Y1, const std::vector<double>& X2);
  void Interpolate(const std::vector<double>& X1, const std::vector<double>& Y1, const std::vector<double>& X2, std::vector<double>& Y2);
  std::vector<double> Interpolate(const std::vector<double>& X1, const std::vector<double>& Y1, const std::vector<double>& X2, const double ExtrapVal);
  void Interpolate(const std::vector<double>& X1, const std::vector<double>& Y1, const std::vector<double>& X2, std::vector<double>& Y2, const double ExtrapVal);
  double Interpolate(const std::vector<double>& X1, const std::vector<double>& Y1, const double& X2);
  
  std::vector<double> SplineIntegral(const std::vector<double>& X1, const std::vector<double>& Y1);
  std::vector<double> SplineIntegral(const std::vector<double>& X1, const std::vector<double>& Y1, const std::vector<double>& X2);
  double SplineCumulativeIntegral(const std::vector<double>& X1, const std::vector<double>& Y1);
  
  class Interpolator {
  protected:
    int n, mm, jsav, cor, dj;
    const std::vector<double> &xx, &yy;
    
  public:
    Interpolator(const std::vector<double>& x, const std::vector<double>& y, int m);
    virtual ~Interpolator() { }//  delete[] xx; delete[] yy; }
    
    double interp(double x) {
      int jlo = cor ? hunt(x) : locate(x);
      return rawinterp(jlo,x);
    }
    
    int locate(const double x);
    int hunt(const double x);
    
    virtual double rawinterp(int jlo, double x) = 0;
  };
  
  
  class PolynomialInterpolator : public Interpolator {
    double dy;
    
  public:
    PolynomialInterpolator(const std::vector<double>& xv, const std::vector<double>& yv, int m)
      : Interpolator(xv,yv,m), dy(0.) { }
    double rawinterp(int jl, double x);
  };
  
  class SplineInterpolator : public Interpolator {
  protected:
    std::vector<double> y2;
    
  public:
    SplineInterpolator(const std::vector<double>& xv, const std::vector<double>& yv, double yp1=1.e99, double ypn=1.e99)
      : Interpolator(xv,yv,2), y2(xv.size())
    { sety2(xv,yv,yp1,ypn); }
    
    void sety2(const std::vector<double>& xv, const std::vector<double>& yv, double yp1, double ypn);
    double rawinterp(int jl, double xv);
  };
  
  class SplineIntegrator : public SplineInterpolator {
  private:
    std::vector<double> IntegrationConstants;
    std::vector<double> IntegrationCoefficients1, IntegrationCoefficients2,
      IntegrationCoefficients3, IntegrationCoefficients4;
    void SetUpIntegrationCoefficients();
    void SetUpIntegrationConstants();
    
  public:
    SplineIntegrator(const std::vector<double>& xv, const std::vector<double>& yv,
		     double yp1=1.e99, double ypn=1.e99)
      : SplineInterpolator(xv, yv, yp1, ypn),
	IntegrationConstants(xv.size()),
	IntegrationCoefficients1(xv.size()),
	IntegrationCoefficients2(xv.size()),
	IntegrationCoefficients3(xv.size()),
	IntegrationCoefficients4(xv.size())
    {
      SetUpIntegrationCoefficients();
      SetUpIntegrationConstants();
    }
    std::vector<double> operator()(const std::vector<double>& x); // Return integral at selected points
    double operator()(const double x); // Return integral at selected points
    std::vector<double> operator()(); // Return integral at all original data points
    double CumulativeIntegral(); // Return the total integral over all original data points
  };
  
} // namespace WaveformUtilities

#endif // INTERPOLATE_HPP
