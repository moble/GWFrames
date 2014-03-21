#ifndef WAVEFORMATAPOINTFT_HPP
#define WAVEFORMATAPOINTFT_HPP

#include <vector>
#include <string>

namespace GWFrames {
  class Waveform;
}

namespace GWFrames {


  /// The WaveformAtAPointFT class is a derived class, constructed
  /// from waveforms evaluated at a point, using the given complex
  /// detector response (F+ + i*Fx) -- or more particularly, its
  /// amplitude and phase.
  class WaveformAtAPointFT {
  private:  // Member data
    double mDt, mVartheta, mVarphi;
    std::vector<double> mRealF, mImagF, mFreqs;
    bool mNormalized;

  public:  // Constructors and Destructor
    WaveformAtAPointFT(const GWFrames::Waveform& W,
                       const double Dt,
                       const double Vartheta,
                       const double Varphi,
                       const double TotalMass, // In solar masses
                       const unsigned int WindowNCycles=1,
                       const double DetectorResponseAmp=1.0,
                       const double DetectorResponsePhase=0.0);

  public: // Access functions
    /// Returns the physical frequencies in Hz
    const std::vector<double>& F() const { return mFreqs; }
    const double& F(const unsigned int f) const { return mFreqs[f]; }
    unsigned int NFreq() const { return mFreqs.size(); }
    bool IsNormalized() const { return mNormalized; }

    const std::vector<double>& Re() const { return mRealF; }
    const std::vector<double>& Im() const { return mImagF; }
    const double& Re(const unsigned int f) const { return mRealF[f]; }
    const double& Im(const unsigned int f) const { return mImagF[f]; }

  public:  // Member functions
    std::vector<double> InversePSD(const std::string& Detector="AdvLIGO_ZeroDet_HighP") const;
    double SNR(const std::vector<double>& InversePSD) const;
    double SNR(const std::string& Detector="AdvLIGO_ZeroDet_HighP") const;
    double InnerProduct(const WaveformAtAPointFT& B,
                        const std::vector<double>& InversePSD) const;
    double Match(const WaveformAtAPointFT& B,
                 const std::vector<double>& InversePSD) const;
    double Match(const WaveformAtAPointFT& B,
                 const std::string& Detector="AdvLIGO_ZeroDet_HighP") const;
    void Match(const WaveformAtAPointFT& B,
               const std::vector<double>& InversePSD,
               double& timeOffset, double& phaseOffset, double& match) const;
    void Match(const WaveformAtAPointFT& B, double& timeOffset,
               double& phaseOffset, double& match,
               const std::string& Detector="AdvLIGO_ZeroDet_HighP") const;
  public:
    WaveformAtAPointFT& Normalize(const std::vector<double>& InversePSD);
    WaveformAtAPointFT& Normalize(const std::string& Detector="AdvLIGO_ZeroDet_HighP");
    WaveformAtAPointFT& ZeroAbove(const double Frequency);
  }; // class

} // namespace GWFrames

#endif // WAVEFORMATAPOINTFT_HPP
