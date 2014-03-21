#ifndef WAVEFORMATAPOINTFT_HPP
#define WAVEFORMATAPOINTFT_HPP

#include "WaveformAtAPoint.hpp"

namespace WaveformObjects {
  
  /// The WaveformAtAPointFT class is a derived class, constructed
  /// from waveforms evaluated at a point, using the given complex
  /// detector response (F+ + i*Fx) -- or more particularly, its
  /// amplitude and phase.
  class WaveformAtAPointFT : public WaveformAtAPoint {
  private:  // Member data
    bool Normalized;
    
  public:  // Constructors and Destructor
    WaveformAtAPointFT();
    WaveformAtAPointFT(const WaveformAtAPoint& W, const unsigned int WindowNCycles=1,
		       const double DetectorResponseAmp=1.0, const double DetectorResponsePhase=0.0);
    ~WaveformAtAPointFT() { }
    
  public: // Access functions
    inline const double F(const unsigned int i) const { return T(i); }
    inline const std::vector<double>& F() const { return T(); }
    
  public:  // Member functions
    double InnerProduct(const WaveformAtAPointFT& B, const std::vector<double>& InversePSD) const;
    double SNR(const std::vector<double>& InversePSD) const;
    double Match(const WaveformAtAPointFT& B, const std::vector<double>& InversePSD) const;
    double Match(const WaveformAtAPointFT& B, const std::string& Detector="AdvLIGO_ZeroDet_HighP") const;
    void Match(const WaveformAtAPointFT& B, const std::vector<double>& InversePSD, double& timeOffset, double& phaseOffset, double& match) const;
    void Match(const WaveformAtAPointFT& B, double& timeOffset, double& phaseOffset, double& match, const std::string& Detector="AdvLIGO_ZeroDet_HighP") const;
    WaveformAtAPointFT& Normalize(const std::vector<double>& InversePSD);
    WaveformAtAPointFT& ZeroAbove(const double Frequency);
    WaveformAtAPointFT operator-(const WaveformAtAPointFT& b) const;
    WaveformAtAPointFT operator*(const double b) const;
  }; // class
  
} // namespace WaveformObjects

std::ostream& operator<<(std::ostream& os, const WaveformObjects::WaveformAtAPointFT& a);

#endif // WAVEFORMATAPOINTFT_HPP
