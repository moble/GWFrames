//////////////////////////////////////
//// Read in the PNWaveform class ////
//////////////////////////////////////
//// Ignore things that don't translate well...
// %ignore operator<<;
// %ignore GWFrames::Waveform::operator=;
//// Get output from PNWaveformModes
%apply std::vector<std::complex<double> >& OUTPUT { std::vector<std::complex<double> >& modes };
//// These will convert the output data to numpy.ndarray for easier use
#ifndef SWIGPYTHON_BUILTIN
%feature("pythonappend") GWFrames::PNWaveform::chi1() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::chi2() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_orb() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_prec() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_tot() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::L() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_orbMag() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_precMag() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_totMag() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::LMag() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Phi_orb() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::chiHat1() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::chiHat2() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::OmegaHat_orb() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::OmegaHat_prec() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::OmegaHat_tot() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::LHat() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
#endif
//// Parse the header file to generate wrappers
%include "../PNWaveforms.hpp"
