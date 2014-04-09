
////////////////////////////////////
//// Read in the Waveform class ////
////////////////////////////////////
// Rename C++ Waveform objects to be _Waveform, so that we can
// subclass Waveform objects in python.
%rename(_Waveform) Waveform;
//// Ignore things that don't translate well...
%ignore operator<<;
%ignore GWFrames::Waveform::operator=;
%ignore GWFrames::Waveforms::operator[];
%rename(__getitem__) GWFrames::Waveforms::operator[] const;
// //// Allow us to get Quaternions returned naturally when passed by reference
// %typemap(in,numinputs=0) Quaternions::Quaternion& OUTPUT (Quaternions::Quaternion temp) { $1 = &temp; }
// %typemap(argout) Quaternions::Quaternion& OUTPUT {
//   $result = SWIG_Python_AppendOutput($result, SWIG_NewPointerObj((new Quaternions::Quaternion(static_cast< const Quaternions::Quaternion& >(temp$argnum))), SWIGTYPE_p_Quaternions__Quaternion, SWIG_POINTER_OWN |  0 ));
// }
//// These will convert the output data to numpy.ndarray for easier use
#ifndef SWIGPYTHON_BUILTIN
%feature("pythonappend") GWFrames::Waveform::T() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Frame() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::LM() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Re() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Im() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Abs() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Arg() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::ArgUnwrapped() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Data() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::DataDot() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Norm() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::LdtVector() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::LLMatrix() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::OShaughnessyEtAlVector() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::AngularVelocityVector() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::CorotatingFrame() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::PNEquivalentOrbitalAV() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::PNEquivalentPrecessionalAV() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
#endif
//// Allow us to extract the outputs naturally in python
%apply double& OUTPUT { double& deltat };
// %apply Quaternions::Quaternion& OUTPUT { Quaternions::Quaternion& R_delta };
//// Parse the header file to generate wrappers
%include "../Waveforms.hpp"
//// Make sure vectors of Waveform are understood
namespace std {
  %template(_vectorW) vector<GWFrames::Waveform>;
};
//// Make any additions to the Waveform class here
%extend GWFrames::Waveform {
  //// This function is called when printing the Waveform object
  std::string __str__() {
    std::stringstream S;
    S << ($self->HistoryStr()) << "#\n"
      << "# <in_python>\n"
      << "import GWFrames\n"
      << "print(this)\n"
      << "# </in_python>" << std::endl << std::setprecision(14);
    for(unsigned int t=0; t<$self->NTimes(); ++t) {
      S << $self->T(t) << " ";
      for(unsigned int mode=0; mode<$self->NModes(); ++mode) {
        S << $self->Re(mode, t) << " " << $self->Im(mode, t) << " ";
      }
      S << std::endl;
    }
    return S.str();
  }
  //// Allow Waveform objects to be pickled
  %insert("python") %{
    def __getstate__(self) :
      return (self.HistoryStr(),
              self.T(),
              self.Frame(),
              self.FrameType(),
              self.DataType(),
              self.RIsScaledOut(),
              self.MIsScaledOut(),
              self.LM(),
              self.Data()
              )
    __safe_for_unpickling__ = True
    def __reduce__(self) :
        return (Waveform, (), self.__getstate__())
    def __setstate__(self, data) :
        self.SetHistory(data[0])
        self.SetTime(data[1])
        self.SetFrame(data[2])
        self.SetFrameType(data[3])
        self.SetDataType(data[4])
        self.SetRIsScaledOut(data[5])
        self.SetMIsScaledOut(data[6])
        self.SetLM(data[7].tolist())
        self.SetData(data[8])
  %}
 };
// Note the 's' on 'Waveforms' below!
%extend GWFrames::Waveforms {
  void __setitem__(int i, const GWFrames::Waveform& W) {
    $self->operator[](i) = W;
    return;
  }
 };
