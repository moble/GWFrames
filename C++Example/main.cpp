#include "Errors.hpp"
#include "Utilities.hpp"
#include "Quaternions.hpp"
#include "Waveforms.hpp"
#include <cmath>
using namespace std;

// To build this program, change any necessary paths in the
// accompanying Makefile, and run 'make'.

// NOTE: This file is a simple example to demonstrate basic usage of a
// few GWFrames commands and direct compilation of the C++ code.  More
// extensive examples are given in the IPython notebook
// 'Documentation.ipynb'.

int main() {
  // Quaternion examples
  // ===================
  // A few random angles we will use below
  const double vartheta = 0.2;
  const double varphi = 0.4;
  const double alpha = 0.3;
  const double beta = 0.5;
  const double gamma = 0.7;
  // Construct some Quaternion objects
  GWFrames::Quaternion R1(vartheta, varphi); // Construct rotor from polar coordinates
  GWFrames::Quaternion R2(alpha, beta, gamma); // Construct rotor from Euler angles
  GWFrames::Quaternion Q(1.2, 3.4, 5.6, 7.8); // Some random Quaternion
  GWFrames::Quaternion R3 = Q.normalized(); // That Quaternion normalized
  // Construct basis vectors as Quaternion objects
  GWFrames::Quaternion x(0.0,1.0,0.0,0.0);
  GWFrames::Quaternion y(0.0,0.0,1.0,0.0);
  GWFrames::Quaternion z(0.0,0.0,0.0,1.0);
  // Print the Quaternions and show what they do
  cout << "\n\n\nQuaternion examples\n===================\n"
       << "x = " << x << "\n"
       << "y = " << y << "\n"
       << "z = " << z << "\n"
       << "\nR1 = " << R1 << " rotates the z basis vector to the polar coordinates (vartheta, varphi)=(" << vartheta << ", " << varphi << "):\n"
       << "R1 * z * R1.conjugate() = \n\t" << R1 * z * R1.conjugate() << "\n"
       << "[0, cos(varphi)*sin(vartheta), sin(varphi)*sin(vartheta), cos(vartheta)] = \n\t[0, "
       << cos(varphi)*sin(vartheta) << ", " << sin(varphi)*sin(vartheta) << ", " << cos(vartheta) << "]\n"
       << "\nR2 = " << R2 << " rotates any vector by Euler angles (alpha, beta, gamma)=(" << alpha << ", " << beta << ", " << gamma << "):\n"
       << "R2 * x * R2.conjugate() = " << R2 * x * R2.conjugate() << "\n"
       << "R2 * y * R2.conjugate() = " << R2 * y * R2.conjugate() << "\n"
       << "R2 * z * R2.conjugate() = " << R2 * z * R2.conjugate() << "\n"
       << "\nMany more uses for Quaternion objects are available.  See Quaternions.hpp or 'Documentation.ipynb' for more possibilities." << endl;
  
  
  
  // Waveform examples
  // =================
  cout << "\n\n\nWaveform examples\n=================\n";
  // Read a Waveform object
  const string InFileName = "../ExampleData/rMPsi4_HighlyPrecessing.dat";
  cout << "Reading Waveform object from '" << InFileName << "'" << endl;
  GWFrames::Waveform W1(InFileName, "MagArg");
  cout << "Finished!\n"
       << "W1.NTimes() = " << W1.NTimes()
       << "\nW1.NModes() = " << W1.NModes()
       << endl;
  // Make a copy of that object
  cout << "\nCopying W1 to W2" << endl;
  GWFrames::Waveform W2(W1);
  cout << "Finished!\n"
       << "W2.NTimes() = " << W2.NTimes()
       << "\nW2.NModes() = " << W2.NModes()
       << endl;
  // Transform W2 to the corotating frame
  vector<int> Lmodes(3);
  Lmodes[0] = 2;
  Lmodes[1] = 3;
  Lmodes[2] = 4;
  cout << "\nTransforming W2 to the corotating frame using modes with ell in {2,3,4}" << endl;
  W2.TransformToCorotatingFrame(Lmodes);
  cout << "Finished!\n";
  // Output the data to a file
  const string OutFileName = "rMPsi4_HighlyPrecessing_Corotating.dat";
  cout << "\nSaving W2 to '" << OutFileName << "'" << endl;
  W2.Output(OutFileName);
  cout << "Finished!\n" << endl;
  
  return 0;
}
