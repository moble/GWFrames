// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#ifndef PNWAVEFORMS_HPP
#define PNWAVEFORMS_HPP

#include "Waveforms.hpp"

namespace GWFrames {
  
  /// Object for calculating a post-Newtonian Waveform with (optional) precession
  class PNWaveform : public Waveform {
    
  public:  // Constructors and Destructor
    PNWaveform();
    PNWaveform(const PNWaveform& W);
    PNWaveform(const double delta, const std::vector<double>& chi1_0, const std::vector<double>& chi2_0, const double Omega_orb_0,
	       const GWFrames::Quaternion& R_0=GWFrames::Quaternion(1,0,0,0));
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
  
} // namespace GWFrames

#endif // PNWAVEFORMS_HPP
