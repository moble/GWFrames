//////////////////////////
//// Import utilities ////
//////////////////////////
%ignore GWFrames::Matrix::operator=;
%ignore GWFrames::Matrix::operator[];
%ignore GWFrames::operator+;
%ignore GWFrames::operator-;
%ignore GWFrames::operator*;
%ignore GWFrames::operator/;
%ignore GWFrames::abs;
%ignore GWFrames::pow;
%include "../Utilities.hpp"
namespace std {
  %template(_vectorM) vector<GWFrames::Matrix>;
  %template() pair<string,string>;
  %template() vector<pair<string,string> >;
};
%extend GWFrames::Matrix {
  // Print the Matrix nicely at the prompt
  %pythoncode{
    def __repr__(self):
        return "".join(
                       [  'Matrix(['+repr([self(0,c) for c in range(self.ncols())])+',\n']
                       + ['        '+repr([self(r,c) for c in range(self.ncols())])+',\n' for r in range(1,self.nrows()-1)]
                       + ['        '+repr([self(r,c) for c in range(self.ncols()) for r in [self.nrows()-1]])+'])'] )
  };
 };
