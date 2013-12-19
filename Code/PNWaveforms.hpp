// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#ifndef PNWAVEFORMS_HPP
#define PNWAVEFORMS_HPP

#include "Waveforms.hpp"

#define PNWaveforms_ellMax 8

namespace GWFrames {

  /// Object for calculating a post-Newtonian Waveform with (optional) precession
  class PNWaveform : public Waveform {

  public:  // Constructors and Destructor
    PNWaveform();
    PNWaveform(const PNWaveform& W);
    PNWaveform(const std::string& Approximant, const double delta, const std::vector<double>& chi1_i, const std::vector<double>& chi2_i,
	       const double Omega_orb_i, const Quaternions::Quaternion& R_frame_i=Quaternions::Quaternion(1,0,0,0),
	       const double PNOrder=4.0, double v_0=-1.0);
    ~PNWaveform() { }

  private:  // Member data
    // std::stringstream history;           // inherited from Waveform
    // std::vector<double> t;               // inherited from Waveform
    // std::vector<Quaternion> frame;       // inherited from Waveform
    // std::vector<std::vector<int> > lm;   // inherited from Waveform
    // MatrixC data;                        // inherited from Waveform
    std::vector<std::vector<double> > mchi1;
    std::vector<std::vector<double> > mchi2;
    std::vector<std::vector<double> > mOmega_orb;
    std::vector<std::vector<double> > mOmega_prec;
    std::vector<std::vector<double> > mL;
    std::vector<double> mPhi_orb;

  public:  // Data access functions
    // Vector a specific time index
    inline const std::vector<double>& chi1(const unsigned int iTime) const { return mchi1[iTime]; }
    inline const std::vector<double>& chi2(const unsigned int iTime) const { return mchi2[iTime]; }
    inline const std::vector<double>& Omega_orb(const unsigned int iTime) const { return mOmega_orb[iTime]; }
    inline const std::vector<double>& Omega_prec(const unsigned int iTime) const { return mOmega_prec[iTime]; }
    std::vector<double> Omega_tot(const unsigned int iTime) const;
    inline const std::vector<double>& L(const unsigned int iTime) const { return mL[iTime]; }
    // Magnitude at a specific time index
    inline double chi1Mag(const unsigned int iTime) const { return GWFrames::abs(mchi1[iTime]); }
    inline double chi2Mag(const unsigned int iTime) const { return GWFrames::abs(mchi2[iTime]); }
    inline double Omega_orbMag(const unsigned int iTime) const { return GWFrames::abs(mOmega_orb[iTime]); }
    inline double Omega_precMag(const unsigned int iTime) const { return GWFrames::abs(mOmega_prec[iTime]); }
    inline double Omega_totMag(const unsigned int iTime) const { return GWFrames::abs(mOmega_orb[iTime]+mOmega_prec[iTime]); }
    inline double LMag(const unsigned int iTime) const { return GWFrames::abs(mL[iTime]); }
    // Direction at a specific time index
    inline std::vector<double> chiHat1(const unsigned int iTime) const { return mchi1[iTime]/GWFrames::abs(mchi1[iTime]); }
    inline std::vector<double> chiHat2(const unsigned int iTime) const { return mchi2[iTime]/GWFrames::abs(mchi2[iTime]); }
    inline std::vector<double> OmegaHat_orb(const unsigned int iTime) const { return mOmega_orb[iTime]/GWFrames::abs(mOmega_orb[iTime]); }
    inline std::vector<double> OmegaHat_prec(const unsigned int iTime) const { return mOmega_prec[iTime]/GWFrames::abs(mOmega_prec[iTime]); }
    inline std::vector<double> OmegaHat_tot(const unsigned int iTime) const { return Omega_tot(iTime)/GWFrames::abs(Omega_tot(iTime)); }
    inline std::vector<double> LHat(const unsigned int iTime) const { return mL[iTime]/GWFrames::abs(mL[iTime]); }
    // Vector at all times
    inline const std::vector<std::vector<double> >& chi1() const { return mchi1; }
    inline const std::vector<std::vector<double> >& chi2() const { return mchi2; }
    inline const std::vector<std::vector<double> >& Omega_orb() const { return mOmega_orb; }
    inline const std::vector<std::vector<double> >& Omega_prec() const { return mOmega_prec; }
    std::vector<std::vector<double> > Omega_tot() const;
    inline const std::vector<std::vector<double> >& L() const { return mL; }
    // Magnitude at all times
    std::vector<double> Omega_orbMag() const;
    std::vector<double> Omega_precMag() const;
    std::vector<double> Omega_totMag() const;
    std::vector<double> LMag() const;
    // Direction at all times
    inline std::vector<std::vector<double> > chiHat1() const { return mchi1/GWFrames::abs(mchi1); }
    inline std::vector<std::vector<double> > chiHat2() const { return mchi2/GWFrames::abs(mchi2); }
    inline std::vector<std::vector<double> > OmegaHat_orb() const { return mOmega_orb/Omega_orbMag(); }
    inline std::vector<std::vector<double> > OmegaHat_prec() const { return mOmega_prec/Omega_precMag(); }
    inline std::vector<std::vector<double> > OmegaHat_tot() const { return Omega_tot()/Omega_totMag(); }
    inline std::vector<std::vector<double> > LHat() const { return mL/LMag(); }
    // Phase
    inline double Phi_orb(const unsigned int iTime) const { return mPhi_orb[iTime]; }
    inline const std::vector<double>& Phi_orb() const { return mPhi_orb; }

  }; // class PNWaveform


  /// Construct PN waveform given dynamics, assuming standard orientation
  class PNWaveformModes {
  private:
    // These constants will be fixed once the mass-difference is known
    const double delta, nu, pownu2, pownu3;
    const double HhatL2M0Rev0, HhatL2M1Imv1, HhatL2M1Imv3, HhatL2M1Imv4, HhatL2M1Imv5, HhatL2M1Rev6, HhatL2M1Imv6, HhatL2M2Rev0, HhatL2M2Rev2, HhatL2M2Rev5, HhatL2M2Imv5, HhatL2M2Rev6, HhatL2M2Imv6, HhatL2M2Rev6lnv, HhatL2M2Rev7, HhatL2M2Imv7, HhatL3M0Imv5, HhatL3M1Imv1, HhatL3M1Imv3, HhatL3M1Imv4, HhatL3M1Imv5, HhatL3M1Rev6, HhatL3M1Imv6, HhatL3M2Rev2, HhatL3M2Rev4, HhatL3M2Rev5, HhatL3M2Imv5, HhatL3M2Rev6, HhatL3M3Imv1, HhatL3M3Imv3, HhatL3M3Imv4, HhatL3M3Imv5, HhatL3M3Rev6, HhatL3M3Imv6, HhatL4M0Rev0, HhatL4M1Imv3, HhatL4M1Imv5, HhatL4M1Rev6, HhatL4M1Imv6, HhatL4M2Rev2, HhatL4M2Rev4, HhatL4M2Rev5, HhatL4M2Imv5, HhatL4M2Rev6, HhatL4M3Imv3, HhatL4M3Imv5, HhatL4M3Rev6, HhatL4M3Imv6, HhatL4M4Rev2, HhatL4M4Rev4, HhatL4M4Rev5, HhatL4M4Imv5, HhatL4M4Rev6, HhatL5M1Imv3, HhatL5M1Imv5, HhatL5M1Rev6, HhatL5M1Imv6, HhatL5M2Rev4, HhatL5M2Rev6, HhatL5M3Imv3, HhatL5M3Imv5, HhatL5M3Rev6, HhatL5M3Imv6, HhatL5M4Rev4, HhatL5M4Rev6, HhatL5M5Imv3, HhatL5M5Imv5, HhatL5M5Rev6, HhatL5M5Imv6, HhatL6M1Imv5, HhatL6M2Rev4, HhatL6M2Rev6, HhatL6M3Imv5, HhatL6M4Rev4, HhatL6M4Rev6, HhatL6M5Imv5, HhatL6M6Rev4, HhatL6M6Rev6, HhatL7M1Imv5, HhatL7M2Rev6, HhatL7M3Imv5, HhatL7M4Rev6, HhatL7M5Imv5, HhatL7M6Rev6, HhatL7M7Imv5, HhatL8M2Rev6, HhatL8M4Rev6, HhatL8M6Rev6, HhatL8M8Rev6;
  public:
    PNWaveformModes(const double idelta);
    void operator()(const double v, const double chis, const double chia, std::vector<std::complex<double> >& modes);
  }; // class PNWaveformModes

} // namespace GWFrames

#endif // PNWAVEFORMS_HPP
