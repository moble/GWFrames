///////////////////////////////////////
//// Import the quantities on scri ////
///////////////////////////////////////
%rename(__getitem__) GWFrames::DataGrid::operator [](unsigned int const) const;
%rename(__setitem__) GWFrames::DataGrid::operator [](unsigned int const);
%rename(__getitem__) GWFrames::Modes::operator [](unsigned int const) const;
%rename(__setitem__) GWFrames::Modes::operator [](unsigned int const);
%rename(__getitem__) GWFrames::SliceModes::operator [](unsigned int const) const;
%rename(__setitem__) GWFrames::SliceModes::operator [](unsigned int const);
%rename(__getitem__) GWFrames::SliceGrid::operator [](unsigned int const) const;
%rename(__setitem__) GWFrames::SliceGrid::operator [](unsigned int const);
%rename(__getitem__) GWFrames::Scri::operator [](unsigned int const) const;
%rename(__setitem__) GWFrames::Scri::operator [](unsigned int const);
%include "../Scri.hpp"
namespace GWFrames {
  %template(SliceOfScriGrid) SliceOfScri<DataGrid>;
  %template(SliceOfScriModes) SliceOfScri<Modes>;
}
%extend GWFrames::DataGrid { // None of the above seem to work, so...
  const std::complex<double> __getitem__(const unsigned int i) const { return $self->operator[](i); }
  void __setitem__(const unsigned int i, const std::complex<double>& a) { $self->operator[](i)=a; }
};
%extend GWFrames::Modes { // None of the above seem to work, so...
  const std::complex<double> __getitem__(const unsigned int i) const { return $self->operator[](i); }
  void __setitem__(const unsigned int i, const std::complex<double>& a) { $self->operator[](i)=a; }
};
%extend GWFrames::SliceModes { // None of the above seem to work, so...
  const GWFrames::Modes& __getitem__(const unsigned int i) const { return $self->operator[](i); }
  void __setitem__(const unsigned int i, const GWFrames::Modes& a) { $self->operator[](i)=a; }
};
%extend GWFrames::SliceGrid { // None of the above seem to work, so...
  const GWFrames::DataGrid& __getitem__(const unsigned int i) const { return $self->operator[](i); }
  void __setitem__(const unsigned int i, const GWFrames::DataGrid& a) { $self->operator[](i)=a; }
};
%extend GWFrames::Scri { // None of the above seem to work, so...
  const GWFrames::SliceModes __getitem__(const unsigned int i) const { return $self->operator[](i); }
  void __setitem__(const unsigned int i, const GWFrames::SliceModes& a) { $self->operator[](i) = a; }
};
%extend GWFrames::SuperMomenta { // None of the above seem to work, so...
  const GWFrames::Modes __getitem__(const unsigned int i) const { return $self->operator[](i); }
  void __setitem__(const unsigned int i, const GWFrames::Modes& a) { $self->operator[](i) = a; }
};
