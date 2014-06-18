%feature("docstring") SliceModes::SuperMomentum """
Find the Moreschi supermomentum.
================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Modes
  
  Description
  -----------
    $\\psi = \\psi_2 + \\sigma \\dot{\\bar{\\sigma}} + \\eth^2 \\bar{\\sigma}$
  
"""

%feature("docstring") GWFrames::Waveform::Re """


  Parameters
  ----------
    const unsigned int Mode
    const unsigned int TimeIndex
  
  Returns
  -------
    double
  

Return vector of real parts of a given mode as function of time.
================================================================
  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    vector<double>
  

Return vector of vector of real parts of all modes as function of time.
=======================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") WaveformUtilities::SplineCumulativeIntegral """


  Parameters
  ----------
    const vector<double>& X1
    const vector<double>& Y1
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::DataGrid::operator[] """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    const complex<double>&
  



  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    complex<double>&
  
"""

%feature("docstring") GWFrames::SuperMomenta::BMSTransform """
Return value of Psi on u'=const slice centered at delta[0].
===========================================================
  Parameters
  ----------
    const Modes& OneOverK
    const Modes& delta
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") InverseConformalFactorFunctor::InverseConformalFactorFunctor """


  Parameters
  ----------
    const ThreeVector& vi
  
  Returns
  -------
    InverseConformalFactorFunctor
  
"""

%feature("docstring") WaveformUtilities::TimeToFrequency """
This function creates a frequency vector in (0 -> positives -> negatives -> 0) order.
=====================================================================================
  Parameters
  ----------
    const vector<double>& Time
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::SliceOfScri<D>::SliceOfScri """


  Parameters
  ----------
    const SliceOfScri& S
  
  Returns
  -------
    SliceOfScri
  
"""

%feature("docstring") GWFrames::Waveform::FrameTypeString """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    string
  
"""

%feature("docstring") GWFrames::Waveform::TransformToAngularVelocityFrame """
Transform Waveform to frame aligned with angular-velocity vector.
=================================================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    This function combines the steps required to obtain the Waveform in the
    frame aligned with the angular-velocity vector. Note that this frame is not
    the corotating frame; this frame has its z axis aligned with the
    angular-velocity vector.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") GWFrames::Waveform::EllMax """
Return greatest ell value present in the data.
==============================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") ellmax """


  Parameters
  ----------
    const int maxindex
  
  Returns
  -------
    int
  
"""

%feature("docstring") Modes::edth2edthbar2 """
The operator edth^2 bar{edth}^2.
================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") GWFrames::PNWaveform::Phi_orb """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    double
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<double>&
  
"""

%feature("docstring") GWFrames::Modes::SetSpin """


  Parameters
  ----------
    const int ess
  
  Returns
  -------
    Modes&
  
"""

%feature("docstring") GWFrames::PNWaveform::Omega_precMag """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    double
  

Vector of magnitudes of Omega_prec at each instant of time.
===========================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::NTimes """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") GWFrames::Waveform::SliceOfTimeIndicesWithEll2 """
Copy of the Waveform between indices i_t_a and i_t_b, only ell=2 modes.
=======================================================================
  Parameters
  ----------
    const unsigned int i_t_a
    unsigned int i_t_b = 0
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    i_t_a and i_t_b should hold the indices pointing to the first time in t
    after t_a, and the first time in t after t_b (or one-past-the-end of t if
    necessary)
  
"""

%feature("docstring") SplineInterpolator::sety2 """


  Parameters
  ----------
    const vector<double>& xv
    const vector<double>& yv
    double yp1
    double ypn
  
  Returns
  -------
    void
  
"""

%feature("docstring") WaveformUtilities::Interpolate """


  Parameters
  ----------
    const vector<double>& X1
    const vector<double>& Y1
    const vector<double>& X2
    vector<double>& Y2
  
  Returns
  -------
    void
  



  Parameters
  ----------
    const vector<double>& X1
    const vector<double>& Y1
    const vector<double>& X2
    vector<double>& Y2
    const double ExtrapVal
  
  Returns
  -------
    void
  



  Parameters
  ----------
    const vector<double>& X1
    const vector<double>& Y1
    const double& X2
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::Re """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<double>&
  



  Parameters
  ----------
    const unsigned int f
  
  Returns
  -------
    const double&
  
"""

%feature("docstring") GWFrames::AlignWaveforms """
Do everything necessary to align two waveform objects.
======================================================
  Parameters
  ----------
    Waveform& A
    Waveform& B
    const double t_1
      Beginning of alignment interval
    const double t_2
      End of alignment interval
    unsigned int InitialEvaluations = 0
      Number of evaluations for dumb initial optimization
    vector<double> nHat_A = vector<double>(0)
      Approximate nHat vector at (t_1+t_2)/2. [optional]
  
  Returns
  -------
    void
  
  Description
  -----------
    This function aligns the frame to the waveform modes for both input
    Waveform objects at time t_mid = (t_1+t_2)/2. It also optimizes the
    alignment of W_B by adjusting its time and overall orientations to align
    with W_A as well as possible. While doing so, it re-adjusts the frame
    alignment to the modes for W_B to account for the changing meaning of t_mid.
    
    Note that t_1 and t_2 refer to fixed times with respect to the time axis of
    W_A.
    
    The input waveforms are transformed to their co-rotating frames if they are
    in the inertial frame. Otherwise, they must already be in the co-rotating
    frame. (E.g., the co-orbital frame is an error.)
    
    The nHat quantity is just the approximate direction for that vector
    (pointing from black hole A to black hole B) in the systems, used to set
    the direction of the x axis for the rotating frame. Only the value at t_mid
    for W_A is needed.
    
    The alignment algorithm assumes that the waveforms are already reasonably
    well aligned in time. In particular, the final value of t_mid+deltat for
    W_B must lie somewhere in the interval (t_1, t_2) at least, and after the
    time shift, W_B must have data over all of that interval.
    
    As long as this last condition is satisfied, and the waveforms are even
    remotely well sampled, and your data does not contain other grievous
    errors, I (Mike Boyle) do hereby guarantee that this function will find the
    optimal alignment. Or your money back.
  
"""

%feature("docstring") GWFrames::Modes::Spin """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::Scri::T """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::SliceOfTimesWithEll2 """
Copy of the Waveform between t_a and t_b, only ell=2 modes.
===========================================================
  Parameters
  ----------
    const double t_a = -1e300
    const double t_b = 1e300
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::NFreq """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") GWFrames::Waveform::HistoryStream """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    stringstream&
  
"""

%feature("docstring") GWFrames::PNWaveform::OmegaHat_orb """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::Waveform::Frame """


  Parameters
  ----------
    const unsigned int TimeIndex
  
  Returns
  -------
    Quaternions::Quaternion
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<Quaternions::Quaternion>&
  
"""

%feature("docstring") GWFrames::PNWaveform::Omega_orbMag """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    double
  

Vector of magnitudes of Omega_orb at each instant of time.
==========================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") WaveformUtilities::SplineIntegrator::SplineIntegrator """


  Parameters
  ----------
    const vector<double>& xv
    const vector<double>& yv
    double yp1 = 1.e99
    double ypn = 1.e99
  
  Returns
  -------
    SplineIntegrator
  
"""

%feature("docstring") GWFrames::Waveform::SetRIsScaledOut """


  Parameters
  ----------
    const bool Scaled
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") SQR """


  Parameters
  ----------
    const double a
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Waveform::RotatePhysicalSystem """
Rotate the physical content of the Waveform by a constant rotor.
================================================================
  Parameters
  ----------
    const Quaternions::Quaternion& R_phys
  
  Returns
  -------
    Waveform&
  

Rotate the physical content of the Waveform.
============================================
  Parameters
  ----------
    vector<Quaternions::Quaternion> R_phys
      Vector of Quaternions by which to rotate
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    This rotates the physical system, leaving the coordinates in place.
    
    The Waveform's frame data records the rotors needed to rotate the standard
    (x,y,z) basis into the (X,Y,Z) basis with respect to which the Waveform
    modes are decomposed. If this is not the first rotation of the frame, we
    need to be careful about how we record the total rotation. Here, we are
    rotating the physical system, while leaving fixed the basis with respect to
    which the modes are decomposed. Therefore, the new frame must be the
    original frame data times $\\bar{R}_{phys}$.
    
    Note that this function does not change the frameType; this is left to the
    calling function.
  
"""

%feature("docstring") GWFrames::PNWaveform """
class GWFrames::PNWaveform
==========================
  Object for calculating a post-Newtonian Waveform with (optional) precession.
  
  Member variables
  ----------------
    vector<vector<double>> mchi1
    vector<vector<double>> mchi2
    vector<vector<double>> mOmega_orb
    vector<vector<double>> mOmega_prec
    vector<vector<double>> mL
    vector<double> mPhi_orb
  
"""

%feature("docstring") GWFrames::Waveform::Data """


  Parameters
  ----------
    const unsigned int Mode
    const unsigned int TimeIndex
  
  Returns
  -------
    complex<double>
  

Return vector of complex data of a given mode as function of time.
==================================================================
  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    vector<complex<double>>
  

Return vector of vector of complex data of all modes as function of time.
=========================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<complex<double>>>
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::SNR """


  Parameters
  ----------
    const vector<double>& InversePSD
      Noise spectrum
  
  Returns
  -------
    double
  



  Parameters
  ----------
    const string& Detector = 'AdvLIGO_ZeroDet_HighP'
      Noise spectrum from this detector
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Waveform::AngularVelocityVector """
Calculate the angular velocity of the Waveform.
===============================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<vector<double>>
  
  Description
  -----------
    This returns the angular velocity of the Waveform, as defined in Sec. II of
    'Angular velocity of gravitational radiation and the corotating frame'.
    Note that the returned vector is relative to the inertial frame.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") SliceModes::Mass """
Calculate the mass of the system from the four-momentum.
========================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Modes """
class GWFrames::Modes
=====================
  Member variables
  ----------------
    int s

            This object holds complex spin-weighted data on the sphere in a
      spin-weighted spherical-harmonic mode representation. All modes are
      present, even for $\\ell<|s|$. The modes are also assumed to be stored in
      order, as $(\\ell,m) = (0,0), (1,-1), (1,0), (1,1), (2,-2), \\ldots$.    int ellMax
    vector<complex<double>> data
  
"""

%feature("docstring") Modes::edth """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Modes
  
  Description
  -----------
    This operator is the one defined by Geroch et al. (1973). It raises the
    spin weight of any field on the sphere by 1, while leaving the boost weight
    unchanged.
    
    This operator is very similar to the basic Newman-Penrose edth operator,
    except that it preserves boost weights. Its effect in this implementation
    is identical (up to a factor of $\\sqrt{2}$) to the NP edth. There is an
    additional term in the definition of the GHP operator, but its value is
    zero. (It transforms nontrivially, though.) In this context, we have
    NPEdth() = sqrt(2)*GHPEdth().
    
    The complex shear $\\sigma$ has spin weight +2 and boost weight +1. The
    radial coordinate $r$ has boost weight -1, and the derivative with respect
    to time $d/du$ has boost weight -1. The asymptotic metric shear $r\\, h$
    has spin weight -2 and boost weight -1. In particular, it seems that $r\\,
    h = r^2\\, \\bar{\\sigma}$.
    
    The Newman-Penrose scalars $\\Psi_i$ have spin weight and boost weight
    equal to $2-i$. (E.g., $\\Psi_4$ has $s = b = -2$.) However, when these are
    multiplied by the appropriate factors of $r$ to find the leading-order
    terms, they acquire boost weights. In particular, we need to multiply
    $\\Psi_i$ by $r^{5-i}$ to get nonzero values at scri, which adds $i-5$ to
    the boost weight, so that the asymptotic NP scalars all have boost weight
    -3.
  
"""

%feature("docstring") Interpolator::locate """


  Parameters
  ----------
    const double x
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames """
namespace GWFrames
==================
"""

%feature("docstring") TransitionFunction_Smooth """
Local utility function.
=======================
  Parameters
  ----------
    const double x
  
  Returns
  -------
    double
  
  Description
  -----------
    This smoothly transitions from 0.0 for x<=0.0 to 1.0 for x>=1.0. The
    function is just the usual transition function based on the familiar smooth
    but non-analytic function.
  
"""

%feature("docstring") zHat """


  Parameters
  ----------
    0 
    0 
    0 
    1 
  
  Returns
  -------
    const Quaternion
  
"""

%feature("docstring") Modes::operator- """


  Parameters
  ----------
    const Modes& M
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") Modes::operator/ """


  Parameters
  ----------
    const Modes& M
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") GWFrames::Waveform::swap """
Efficiently swap data between two Waveform objects.
===================================================
  Parameters
  ----------
    Waveform& b
  
  Returns
  -------
    void
  
  Description
  -----------
    This function uses the std::vector method 'swap' which simply swaps
    pointers to data, for efficiency.
  
"""

%feature("docstring") GWFrames::Modes::operator= """


  Parameters
  ----------
    const Modes& B
  
  Returns
  -------
    Modes&
  
"""

%feature("docstring") Modes::operator* """


  Parameters
  ----------
    const Modes& M
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") Modes::operator+ """


  Parameters
  ----------
    const Modes& M
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") AdvLIGO_ZeroDet_HighP """


  Parameters
  ----------
    const vector<double>& F
    const bool Invert = false
    const double NoiseFloor = 0.0
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::IsNormalized """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    bool
  
"""

%feature("docstring") WaveformUtilities """
namespace WaveformUtilities
===========================
"""

%feature("docstring") AdvLIGO_ZeroDet_LowP """


  Parameters
  ----------
    const vector<double>& F
    const bool Invert = false
    const double NoiseFloor = 0.0
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform """
class GWFrames::Waveform
========================
  Object storing data and other information for a single waveform.
  
  Member variables
  ----------------
    int spinweight
    int boostweight
    stringstream history
    vector<double> t
    vector<Quaternions::Quaternion> frame
    WaveformFrameType frameType
    WaveformDataType dataType
    bool rIsScaledOut
    bool mIsScaledOut
    vector<vector<int>> lm
    MatrixC data
  
  Non-public member functions
  ---------------------------
    Waveform& TransformModesToRotatedFrame
      Rotate modes of the Waveform object.  
"""

%feature("docstring") GWFrames::PNWaveform::chi2Mag """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Waveform::SetFrameType """


  Parameters
  ----------
    const WaveformFrameType Type
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") Scri::BMSTransformation """
Apply a (constant) BMS transformation to data on null infinity.
===============================================================
  Parameters
  ----------
    const double& u0
      Initial time slice to transform
    const ThreeVector& v
      Three-vector of the boost relative to the current frame
    const Modes& delta
      Spherical-harmonic modes of the supertranslation
  
  Returns
  -------
    SliceModes
  
  Description
  -----------
    A general BMS transformation is expressed as a conformal transformation of
    the sphere (which encompasses rotations and boosts) and a supertranslation
    (which affects only the time coordinate on null infinity). Here, the
    conformal transformation of the sphere is assumed to be a simple boost,
    described by the velocity vector v. (The other freedom in the conformal
    group is a simple rotation, which we assume is zero.) The supertranslation
    is decomposed into (scalar) spherical harmonics. The input data is assumed
    to satisfy conditions on the nonzero data such that the values everywhere
    on the sphere are real.
    
    The work to be done by this function includes (1) evaluating the data on
    the appropriate equi-angular grid of the final frame, while simultaneously
    transforming the data at each point by the appropriate spin factor, boost
    factor, and mixing; (2) interpolating those data to the appropriate values
    of retarded time of the final frame; and (3) transforming back to spectral
    space to store the data in their usual representation.
    
    The relation between the new and old time coordinates is $u' =
    K(u-\\delta)$, where $K$ is the conformal factor (which is a function of
    angle). So we need to interpolate at each point to the original time
    coordinate $u = u'/K + \\delta$, which again depends on angle. That is, we
    have to interpolate to a different time for each grid point. In this case,
    we arbitrarily set $u' = 0$, because any other choice can be absorbed into
    a time- and space-translation. This does not matter, of course, because
    that choice is not stored in any way.
  
"""

%feature("docstring") GWFrames::Waveform::SetDataType """


  Parameters
  ----------
    const WaveformDataType Type
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::Waveform::MaxNormIndex """
Return the data index corresponding to the time of the largest norm.
====================================================================
  Parameters
  ----------
    const unsigned int SkipFraction = 4
      Integer fraction of data to skip before looking
  
  Returns
  -------
    unsigned int
  
  Description
  -----------
    The default value of SkipFraction is 4, meaning that we start looking for
    the maximum after 1/4th of the data, so as to cut out junk radiation. Note
    that this is integer division, so an argument of NTimes()+1 will look
    through all of the data.
    
    Norm()
    
    MaxNormTime()
  
"""

%feature("docstring") SplineIntegrator::operator() """


  Parameters
  ----------
    const vector<double>& x
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const double x
  
  Returns
  -------
    double
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::Waveform """
Default constructor for an empty object.
========================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  

Copy constructor.
=================
  Parameters
  ----------
    const Waveform& W
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    Simply copies all fields in the input object to the constructed object,
    including history
  

Constructor from data file.
===========================
  Parameters
  ----------
    const string& FileName
      Relative path to data file
    const string& DataFormat
      Either 'ReIm' or 'MagArg'
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    NOTE: This function assumes that the data are stored as (ell,m) modes,
    starting with (2,-2), incrementing m, then incrementing ell and starting
    again at m=-ell. If this is not how the modes are stored, the 'lm' data of
    this object needs to be reset or bad things will happen when trying to find
    the angular-momentum vector or rotate the waveform.
  

Explicit constructor from data.
===============================
  Parameters
  ----------
    const vector<double>& T
    const vector<vector<int>>& LM
    const vector<vector<complex<double>>>& Data
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    Arguments are T, LM, Data, which consist of the explicit data.
  
"""

%feature("docstring") GWFrames::Waveform::SetTime """


  Parameters
  ----------
    const vector<double>& a
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::Waveform::SliceOfTimeIndices """
Copy the Waveform between indices i_t_a and i_t_b.
==================================================
  Parameters
  ----------
    const unsigned int i_t_a
      Index of initial time
    unsigned int t_b = 0
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    i_t_a and i_t_b should hold the indices pointing to the first time in t
    after t_a, and the first time in t after t_b (or one-past-the-end of t if
    necessary)
  
"""

%feature("docstring") WaveformUtilities::TimeToPositiveFrequencies """
This function returns the positive half of the frequencies, so returned size is 1/2 input size + 1.
===================================================================================================
  Parameters
  ----------
    const vector<double>& Time
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::SuperMomenta::MoreschiIteration """
Transform to given slice with given BMS transformation, and return next step in Moreschi algorithm.
===================================================================================================
  Parameters
  ----------
    Modes& OneOverK
      Inverse conformal factor (input/output)
    Modes& delta
      Supertranslation (input/output)
  
  Returns
  -------
    void
  
  Description
  -----------
    This function first transforms Psi to a u'=constant slice centered at the
    value delta[0] with the given BMS transformation. It then replaces the
    values of that BMS transformation with the next step in the Moreschi
    algorithm.
  
"""

%feature("docstring") GWFrames::Waveform::DataDot """
Return time derivative of data.
===============================
  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    vector<complex<double>>
  
"""

%feature("docstring") GWFrames::Waveform::Differentiate """
Differentiate the waveform as a function of time.
=================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") InverseConformalFactorFunctor::operator() """


  Parameters
  ----------
    const Quaternions::Quaternion& R
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::F """
Returns the physical frequencies in Hz.
=======================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<double>&
  



  Parameters
  ----------
    const unsigned int f
  
  Returns
  -------
    const double&
  
"""

%feature("docstring") GWFrames::Waveform::ZParityViolationNormalized """
Measure the relative magnitude of the violation of parity in the z direction.
=============================================================================
  Parameters
  ----------
    vector<int> Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    This function measures the violation of invariance under z-parity
    (reflection across the x-y plane). Nonprecessing systems in a suitable
    frame should have zero violation. Precessing systems in any frame and
    nonprecessing systems in the wrong frame will show violations. This
    quantity can be minimized over orientation to show the presence or absence
    of a plane of symmetry.
    
    The quantity is normalized by the overall norm of the data at each instant,
    and the square-root of that ratio is returned.
  
"""

%feature("docstring") Modes::Modes """


  Parameters
  ----------
    const int spin
    const vector<complex<double>>& Data
  
  Returns
  -------
    Modes
  



  Parameters
  ----------
    DataGrid D
    const int L = -1
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") GWFrames::SliceModes::BMSTransformationOnSlice """
Execute a BMS transformation except for the supertranslation of points.
=======================================================================
  Parameters
  ----------
    const double u
    const ThreeVector& v
    const Modes& delta
  
  Returns
  -------
    SliceGrid
  
  Description
  -----------
    A full BMS transformation is only possible using information from multiple
    slices due to the supertranslation moving points 'between slices'. This
    function simply transforms the data within the slice by accounting for the
    change of grid at each point, and the change of grid points themselves. The
    returned object is a DataGrid object, each point of which can then be used
    to interpolate to the supertranslated time.
  
"""

%feature("docstring") GWFrames::Waveform::BinaryOp """


  Parameters
  ----------
    typename Op 
    const Waveform& b
  
  Returns
  -------
    typename Op
  
"""

%feature("docstring") Modes::bar """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") SliceModes::MoreschiIteration """
Find the next iteration of the BMS transformation via Moreschi's algorithm.
===========================================================================
  Parameters
  ----------
    Modes& OneOverK_ip1
      Inverse conformal factor for the next step
    Modes& delta_ip1
      Supertranslation for the next step
  
  Returns
  -------
    void
  
  Description
  -----------
    This member function applies to a SliceModes object that has already been
    transformed by the BMS transformation represented by $K_i$ and $\\delta_i$.
    This then takes the data on that slice and computes the values of $K_{i+1}$
    and $\\delta_{i+1}$, returning them by reference.
  
"""

%feature("docstring") GWFrames::Waveform::SpinWeight """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") SplineIntegrator::CumulativeIntegral """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Modes::Modes """


  Parameters
  ----------
    const int size = 0
  
  Returns
  -------
    Modes
  



  Parameters
  ----------
    const Modes& A
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") StringForm """


  Parameters
  ----------
    const vector<int>& Lmodes
  
  Returns
  -------
    string
  
"""

%feature("docstring") tolower """


  Parameters
  ----------
    const string& A
  
  Returns
  -------
    string
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::Im """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<double>&
  



  Parameters
  ----------
    const unsigned int f
  
  Returns
  -------
    const double&
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::WaveformAtAPointFT """


  Parameters
  ----------
    const Waveform& W
    const double Dt
    const double Vartheta
    const double Varphi
    const double TotalMass
    const unsigned int WindowNCycles = 1
    const double DetectorResponseAmp = 1.0
    const double DetectorResponsePhase = 0.0
  
  Returns
  -------
    WaveformAtAPointFT
  
"""

%feature("docstring") WaveformUtilities::realdft """


  Parameters
  ----------
    vector<double>& data
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::PNWaveform::PNWaveform """
Default constructor for an empty object.
========================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    PNWaveform
  

Copy constructor.
=================
  Parameters
  ----------
    const PNWaveform& W
  
  Returns
  -------
    PNWaveform
  
  Description
  -----------
    Simply copies all fields in the input object to the constructed object,
    including history
  

Constructor of PN waveform from parameters.
===========================================
  Parameters
  ----------
    const string& Approximant
      'TaylorT1'|'TaylorT4'|'TaylorT5'
    const double delta
      Normalized BH mass difference (M1-M2)/(M1+M2)
    const vector<double>& chi1_i
      Initial dimensionless spin vector of BH1
    const vector<double>& chi2_i
      Initial dimensionless spin vector of BH2
    const double Omega_orb_i
      Initial orbital angular frequency
    double Omega_orb_0 = -1.0
      Earliest orbital angular frequency to compute (optional)
    const Quaternions::Quaternion& R_frame_i = Quaternions::Quaternion(1, 0, 0, 0)
      Overall rotation of the system (optional)
    const double PNOrder = 4.0
      PN order at which to compute all quantities (default: 4.0)
  
  Returns
  -------
    PNWaveform
  
  Description
  -----------
    The PN system is initialized having the BHs along the x axis, with the
    orbital angular velocity along the positive z axis, having magnitude
    Omega_orb_i. The input spin vectors must be defined with respect to this
    basis. Note that the optional parameter Omega_orb_0 may be given, in which
    case the PN system is also evolved backwards to that point.
    
    The TaylorTn system is first integrated to compute the dynamics of the
    binary. The evolved spin vectors chi1 and chi2, orbital angular-velocity
    vector Omega_orb, and orbital phase Phi_orb are stored. Simultaneously, the
    minimal-rotation frame of the angular-velocity vector is computed, then
    rotated about the z' axis by Phi_orb, resulting in the binary's frame. Once
    this step is completed, the information is used to construct the waveform
    in the minimal-rotation frame. (That is, the waveform will be essentially
    corotating.)
    
    Note that, to get the PNWaveform in an inertial frame, you must first apply
    the method TransformToInertialFrame().
  
"""

%feature("docstring") GWFrames::Waveform::operator() """


  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    const complex<double> *
  
"""

%feature("docstring") GWFrames::Waveform::Arg """


  Parameters
  ----------
    const unsigned int Mode
    const unsigned int TimeIndex
  
  Returns
  -------
    double
  

Return vector of arg of a given mode as function of time.
=========================================================
  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    Note that this quantity is not 'unwrapped'. That is, the arg is between -pi
    and +pi. To get a smooth, continuous phase in python, use numpy.unwrap.
  

Return vector of vector of arg of all modes as function of time.
================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::Waveform::DescriptorString """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    string
  
"""

%feature("docstring") SliceModes::EllMax """
Find largest ell value in the data on this slice.
=================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") SliceOfScri::SliceOfScri """
Empty constructor with reserved storage.
========================================
  Parameters
  ----------
    const int size = 0
  
  Returns
  -------
    SliceOfScri
  
"""

%feature("docstring") GWFrames::Waveform::LM """


  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    const vector<int>&
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<vector<int>>&
  
"""

%feature("docstring") GWFrames::Waveform::InterpolateToPoint """
Evaluate Waveform at a particular sky location and an instant of time.
======================================================================
  Parameters
  ----------
    const double vartheta
      Polar angle of detector
    const double varphi
      Azimuthal angle of detector
    const double t_i
      New time to interpolate to
  
  Returns
  -------
    complex<double>
  
  Description
  -----------
    Note that the input angle parameters are measured relative to the binary's
    coordinate system. In particular, this will make no sense if the frame type
    is something other than inertial, and will fail if the FrameType is neither
    UnknownFrameType nor Inertial.
  
"""

%feature("docstring") GWFrames::InverseConformalFactorGrid """
Construct a grid with the conformal factor at each point.
=========================================================
  Parameters
  ----------
    const ThreeVector& v
    const int n_theta
    const int n_phi
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::Waveform::Hybridize """
Hybridize this Waveform with another.
=====================================
  Parameters
  ----------
    const Waveform& B
      Second Waveform to hybridize with
    const double t1
      Beginning of time over which to transition
    const double t2
      End of time over which to transition
    const double tMinStep = 0.005
      Lower limit on time step appearing in the output
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This function simply takes two Waveforms and blends them together. In
    particular, it does not align the Waveforms; that is assumed to have been
    done already.
    
    The transition function is a $C^\\infty$ function, meaning that the output
    data has exactly this Waveform's data before t1, exactly Waveform B's data
    after t2, and a smooth blend in between.
    
    Note that this function does NOT operate in place; a new Waveform object is
    constructed and returned.
  
"""

%feature("docstring") GWFrames::Modes::size """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") GWFrames::PNWaveform::LMag """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    double
  

Vector of magnitudes of angular momentum L at each instant of time.
===================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::operator= """
Assignment operator.
====================
  Parameters
  ----------
    const Waveform& 
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::ScriFunctor::operator() """


  Parameters
  ----------
    const Quaternions::Quaternion& 
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::SuperMomenta::SuperMomenta """


  Parameters
  ----------
    const unsigned int size
  
  Returns
  -------
    SuperMomenta
  



  Parameters
  ----------
    const SuperMomenta& S
  
  Returns
  -------
    SuperMomenta
  



  Parameters
  ----------
    const vector<double>& T
    const vector<Modes>& psi
  
  Returns
  -------
    SuperMomenta
  



  Parameters
  ----------
    const Scri& scri
  
  Returns
  -------
    SuperMomenta
  
"""

%feature("docstring") GWFrames::Waveform::BoostWeight """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::Waveform::Abs """


  Parameters
  ----------
    const unsigned int Mode
    const unsigned int TimeIndex
  
  Returns
  -------
    double
  

Return vector of absolute value of a given mode as function of time.
====================================================================
  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    vector<double>
  

Return vector of vector of absolute value of all modes as function of time.
===========================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::Waveform::NModes """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") GWFrames::Scri """
class GWFrames::Scri
====================
  Member variables
  ----------------
    vector<double> t

            A Scri object contains all the information needed to express the
      asymptotic geometry of null infinity, and symmetry transformations
      thereof. This geometry is encoded in the NewmanPenrose curvature scalars
      $\\psi_0, \\ldots, \\psi_4$ and the complex shear $\\sigma$ of outgoing
      null rays (strain in gravitational-wave detectors). The general symmetry
      transformation is an element of the BondiMetznerSachs (BMS) group, which
      transforms the data contained by Scri among itself.    vector<SliceModes> slices
  
"""

%feature("docstring") Interpolator::Interpolator """


  Parameters
  ----------
    const vector<double>& x
    const vector<double>& y
    int m
  
  Returns
  -------
    Interpolator
  
"""

%feature("docstring") GWFrames::Waveform::DataTypeString """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    string
  
"""

%feature("docstring") GWFrames::Waveform::operator- """


  Parameters
  ----------
    const Waveform& B
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") GWFrames::Waveform::operator/ """


  Parameters
  ----------
    const Waveform& B
  
  Returns
  -------
    Waveform
  



  Parameters
  ----------
    const double b
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") GWFrames::Waveform::operator* """


  Parameters
  ----------
    const Waveform& B
  
  Returns
  -------
    Waveform
  



  Parameters
  ----------
    const double b
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") GWFrames::Waveform::TransformToInertialFrame """
Transform Waveform to an inertial frame.
========================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    This function uses the stored frame information to transform from whatever
    rotating frame the waveform is currently in, to a stationary, inertial
    frame. This is the usual frame of scri^+, and is the frame in which GW
    observations should be made.
  
"""

%feature("docstring") GWFrames::Waveform::operator+ """


  Parameters
  ----------
    const Waveform& B
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") WaveformUtilities::Interpolator """
class WaveformUtilities::Interpolator
=====================================
  Member variables
  ----------------
    int n
    int mm
    int jsav
    int cor
    int dj
    const vector<double>& xx
    const vector<double>& yy
  
"""

%feature("docstring") Modes::edthbar """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Modes
  
  Description
  -----------
    This operator is the one defined by Geroch et al. (1973). It lowers the
    spin weight of any field on the sphere by 1, while leaving the boost weight
    unchanged.
    
    This operator is very similar to the basic Newman-Penrose edth operator,
    except that it preserves boost weights. Its effect in this implementation
    is identical (up to a factor of $\\sqrt{2}$) to the NP edth. There is an
    additional term in the definition of the GHP operator, but its value is
    zero. (It transforms nontrivially, though.) In this context, we have
    NPEdthBar() = sqrt(2)*GHPEdthBar().
    
    The complex shear $\\sigma$ has spin weight +2 and boost weight +1. The
    radial coordinate $r$ has boost weight -1, and the derivative with respect
    to time $d/du$ has boost weight -1. The asymptotic metric shear $r\\, h$
    has spin weight -2 and boost weight -1. In particular, it seems that $r\\,
    h = r^2\\, \\bar{\\sigma}$.
    
    The Newman-Penrose scalars $\\Psi_i$ have spin weight and boost weight
    equal to $2-i$. (E.g., $\\Psi_4$ has $s = b = -2$.) However, when these are
    multiplied by the appropriate factors of $r$ to find the leading-order
    terms, they acquire boost weights. In particular, we need to multiply
    $\\Psi_i$ by $r^{5-i}$ to get nonzero values at scri, which adds $i-5$ to
    the boost weight, so that the asymptotic NP scalars all have boost weight
    -3.
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT """
class GWFrames::WaveformAtAPointFT
==================================
  The WaveformAtAPointFT class is a derived class, constructed from waveforms
  evaluated at a point, using the given complex detector response (F+ + i*Fx) 
  or more particularly, its amplitude and phase.
  
  Member variables
  ----------------
    double mDt
    double mVartheta
    double mVarphi
    vector<double> mRealF
    vector<double> mImagF
    vector<double> mFreqs
    bool mNormalized
  
"""

%feature("docstring") GWFrames::Waveform::SetHistory """


  Parameters
  ----------
    const string& Hist
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::SuperMomenta::NTimes """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::PNWaveform::OmegaHat_tot """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::ConformalFactorGrid """
Construct a grid with the conformal factor at each point.
=========================================================
  Parameters
  ----------
    const ThreeVector& v
    const int n_theta
    const int n_phi
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") WaveformUtilities::dft """


  Parameters
  ----------
    vector<double>& data
  
  Returns
  -------
    void
  
  Description
  -----------
    The following call the Numerical Recipes fft routines (from fourier.h) Note
    that the returned quantities represent the bare fft sum, with no
    normalization constants
  
"""

%feature("docstring") GWFrames::Waveform::SliceOfTimesWithoutModes """
Copy of the Waveform between t_a and t_b without mode data.
===========================================================
  Parameters
  ----------
    const double t_a = -1e300
    const double t_b = 1e300
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") WaveformUtilities::WrapVecDoub """
class WaveformUtilities::WrapVecDoub
====================================
  This helper class is useful for untangling NR's fft routines.
  
  Member variables
  ----------------
    vector<double> vvec
    vector<double>& v
    int n
    int mask
  
"""

%feature("docstring") GWFrames::Modes::SetEllMax """


  Parameters
  ----------
    const int ell
  
  Returns
  -------
    Modes&
  
"""

%feature("docstring") GWFrames::Waveform::TransformToCorotatingFrame """
Transform Waveform to corotating frame.
=======================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    This function combines the steps required to obtain the Waveform in the
    corotating frame. Note that this leaves an integration constant unset. To
    set it, the modes should be rotated so that they are aligned with the frame
    using AlignModesToFrame.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") GWFrames::Waveform::SetT """


  Parameters
  ----------
    const vector<double>& a
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::Boost """


  Parameters
  ----------
    ThreeVector v
    ThreeVector n
  
  Returns
  -------
    Quaternions::Quaternion
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::Match """
Compute the match between two WaveformAtAPointFT.
=================================================
  Parameters
  ----------
    const WaveformAtAPointFT& B
    const vector<double>& InversePSD
  
  Returns
  -------
    double
  

Compute the match between two WaveformAtAPointFT.
=================================================
  Parameters
  ----------
    const WaveformAtAPointFT& B
    const string& Detector = 'AdvLIGO_ZeroDet_HighP'
  
  Returns
  -------
    double
  

Compute the match between two WaveformAtAPointFT.
=================================================
  Parameters
  ----------
    const WaveformAtAPointFT& B
      WaveformAtAPointFT to compute match with
    const vector<double>& InversePSD
      Spectrum used to weight contributions by frequencies to match
    double& timeOffset
      Time offset used between the waveforms
    double& phaseOffset
      Phase offset used between the waveforms
    double& match
      Match between the two waveforms
  
  Returns
  -------
    void
  
  Description
  -----------
    The return from ifft is just the bare FFT sum, so we multiply by df to get
    the continuum-analog FT. This is correct because the input data (re,im) are
    the continuum-analog data, rather than just the return from the bare FFT
    sum. See, e.g., Eq. (A.33) [rather than Eq. (A.35)] of
    http://etd.caltech.edu/etd/available/etd-01122009-143851/
  

Compute the match between two WaveformAtAPointFT.
=================================================
  Parameters
  ----------
    const WaveformAtAPointFT& B
    double& timeOffset
    double& phaseOffset
    double& match
    const string& Detector = 'AdvLIGO_ZeroDet_HighP'
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::FindModeIndex """
Find index of mode with given (l,m) data.
=========================================
  Parameters
  ----------
    const int L
    const int M
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") WaveformUtilities::WrapVecDoub::imag """


  Parameters
  ----------
    int i
  
  Returns
  -------
    double&
  
"""

%feature("docstring") GWFrames::DataGrid::SetNPhi """


  Parameters
  ----------
    const int N_phi
  
  Returns
  -------
    DataGrid&
  
"""

%feature("docstring") GWFrames::Waveform::~Waveform """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    ~Waveform
  
"""

%feature("docstring") PolynomialInterpolator::rawinterp """


  Parameters
  ----------
    int jl
    double x
  
  Returns
  -------
    double
  
"""

%feature("docstring") WaveformUtilities::idft """


  Parameters
  ----------
    vector<double>& data
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::PNWaveform::chiHat2 """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::PNWaveform::chiHat1 """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::SuperMomenta::operator[] """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    const Modes
  



  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    Modes&
  
"""

%feature("docstring") GWFrames::Waveform::ArgUnwrapped """


  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    vector<double>
  

Return vector of vector of arg of all modes as function of time.
================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::Waveform::InterpolateInPlace """
Interpolate the Waveform to a new set of time instants.
=======================================================
  Parameters
  ----------
    const vector<double>& NewTime
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::operator+ """


  Parameters
  ----------
    const double& a
    const DataGrid& b
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::Waveform::T """


  Parameters
  ----------
    const unsigned int TimeIndex
  
  Returns
  -------
    double
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<double>&
  
"""

%feature("docstring") GWFrames::Waveform::KeepOnlyEll2 """
Remove data relating to all but the ell=2 modes.
================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::Waveform::operator() """


  Parameters
  ----------
    const unsigned int Mode
    const unsigned int TimeIndex
  
  Returns
  -------
    complex<double>
  
"""

%feature("docstring") GWFrames::SliceOfScri<D>::operator[] """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    const D&
  



  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    D&
  
"""

%feature("docstring") GWFrames::operator* """


  Parameters
  ----------
    const double& a
    const DataGrid& b
  
  Returns
  -------
    DataGrid
  



  Parameters
  ----------
    const double b
    const Waveform& A
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") GWFrames::Waveform::DropTimesOutside """
Remove all data relating to times outside of the given range.
=============================================================
  Parameters
  ----------
    const double ta
    const double tb
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::PNWaveform::Omega_orb """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    const vector<double>&
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<vector<double>>&
  
"""

%feature("docstring") WaveformUtilities::NoiseCurve """


  Parameters
  ----------
    const vector<double>& F
    const string& Detector = 'AdvLIGO_ZeroDet_HighP'
    const bool Invert = false
    const double NoiseFloor = 0.0
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    These functions return various noise curves as described e.g., in
    LIGO-T0900288-v3, and fit by Collin Capano. Noise curves implemented thus
    far are: IniLIGO_Approx  The Initial LIGO design goal AdvLIGO_ZeroDet_HighP
     LIGO-T0900288-v3 AdvLIGO_ZeroDet_LowP  LIGO-T0900288-v3
    AdvLIGO_NSNSOptimal  Collin Capano's fit for the NS-NS optimized noise curve
  
"""

%feature("docstring") AdvLIGO_NSNSOptimal """


  Parameters
  ----------
    const vector<double>& F
    const bool Invert = false
    const double NoiseFloor = 0.0
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::PNWaveform::Omega_tot """
Total angular velocity of PN binary at an instant of time.
==========================================================
  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    vector<double>
  

Total angular velocity of PN binary at each instant of time.
============================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::SuperMomenta """
class GWFrames::SuperMomenta
============================
  Member variables
  ----------------
    vector<double> t
    vector<Modes> Psi
  
"""

%feature("docstring") GWFrames::Scri::Scri """


  Parameters
  ----------
    const Scri& S
  
  Returns
  -------
    Scri
  
"""

%feature("docstring") GWFrames::Waveform::SetMIsScaledOut """


  Parameters
  ----------
    const bool Scaled
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") WaveformUtilities::InverseNoiseCurve """


  Parameters
  ----------
    const vector<double>& F
    const string& Detector = 'AdvLIGO_ZeroDet_HighP'
    const double NoiseFloor = 0.0
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::LdtVector """
Calculate the $<L \\partial_t>$ quantity defined in the paper.
==============================================================
  Parameters
  ----------
    vector<int> Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<vector<double>>
  
  Description
  -----------
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
    
    $<L \\partial_t>^a = \\sum_{\\ell,m,m'} \\Im [ \\bar{f}^{\\ell,m'}
    <\\ell,m' | L_a | \\ell,m> \\dot{f}^{\\ell,m} ]$
  
"""

%feature("docstring") DataGrid::operator/ """


  Parameters
  ----------
    const DataGrid& 
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") WaveformUtilities::SplineIntegral """


  Parameters
  ----------
    const vector<double>& X1
    const vector<double>& Y1
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const vector<double>& X1
    const vector<double>& Y1
    const vector<double>& X2
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") DataGrid::operator- """


  Parameters
  ----------
    const DataGrid& 
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::DataGrid::N_phi """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") DataGrid::operator+ """


  Parameters
  ----------
    const DataGrid& 
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") DataGrid::operator* """


  Parameters
  ----------
    const DataGrid& 
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::Waveform::TransformToCoprecessingFrame """
Transform Waveform to co-precessing frame.
==========================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    This function combines the steps required to obtain the Waveform in a frame
    with its z axis aligned with the vector V_h defined by O'Shaughnessy et al.
    [2011], and its alignment about the z axis satisfying the minimal-rotation
    condition, as given by Boyle, Owen, Pfeiffer [2011].
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") WaveformUtilities::WrapVecDoub::validate """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    void
  
"""

%feature("docstring") InverseConformalFactorFunctor """
class InverseConformalFactorFunctor
===================================
  Member variables
  ----------------
    double gamma
    Quaternion v
  
"""

%feature("docstring") GWFrames::Waveform::Interpolate """
Interpolate the Waveform to a new set of time instants.
=======================================================
  Parameters
  ----------
    const vector<double>& NewTime
      New vector of times to which this interpolates
    const bool AllowTimesOutsideCurrentDomain = false
      [Default: false]
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    If AllowTimesOutsideCurrentDomain is true, the values of all modes will be
    set to 0.0 for times outside the current set of time data. If false, and
    such times are requested, an error will be thrown.
  
"""

%feature("docstring") GWFrames::Waveform::HistoryStr """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    string
  
"""

%feature("docstring") WaveformUtilities::WrapVecDoub::operator[] """


  Parameters
  ----------
    int i
  
  Returns
  -------
    complex<double>&
  
"""

%feature("docstring") GWFrames::Scri::operator[] """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    const SliceModes
  



  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    SliceModes&
  
"""

%feature("docstring") Scri::Scri """


  Parameters
  ----------
    const Waveform& psi0
    const Waveform& psi1
    const Waveform& psi2
    const Waveform& psi3
    const Waveform& psi4
    const Waveform& sigma
  
  Returns
  -------
    Scri
  
"""

%feature("docstring") GWFrames::vFromOneOverK """
Derive three-velocity from the inverse conformal metric.
========================================================
  Parameters
  ----------
    const Modes& OneOverK
  
  Returns
  -------
    ThreeVector
  
"""

%feature("docstring") GWFrames::Waveform::SetData """


  Parameters
  ----------
    const vector<vector<complex<double>>>& a
  
  Returns
  -------
    Waveform&
  



  Parameters
  ----------
    const unsigned int i_Mode
    const unsigned int i_Time
    const complex<double>& a
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") complexi """


  Parameters
  ----------
    0. 0
    1. 0
  
  Returns
  -------
    const complex<double>
  
"""

%feature("docstring") VectorStringForm """


  Parameters
  ----------
    const vector<double>& V
  
  Returns
  -------
    string
  
"""

%feature("docstring") Interpolator::hunt """


  Parameters
  ----------
    const double x
  
  Returns
  -------
    int
  
"""

%feature("docstring") WaveformUtilities::Interpolate """


  Parameters
  ----------
    const vector<double>& X1
    const vector<double>& Y1
    const vector<double>& X2
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const vector<double>& X1
    const vector<double>& Y1
    const vector<double>& X2
    const double ExtrapVal
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") WaveformUtilities::SplineIntegrator """
class WaveformUtilities::SplineIntegrator
=========================================
  Member variables
  ----------------
    vector<double> IntegrationConstants
    vector<double> IntegrationCoefficients1
    vector<double> IntegrationCoefficients2
    vector<double> IntegrationCoefficients3
    vector<double> IntegrationCoefficients4
  
  Non-public member functions
  ---------------------------
    void SetUpIntegrationCoefficients
    void SetUpIntegrationConstants
  
"""

%feature("docstring") WaveformUtilities::WrapVecDoub::real """


  Parameters
  ----------
    int i
  
  Returns
  -------
    double&
  
"""

%feature("docstring") GWFrames::ScriFunctor """
class GWFrames::ScriFunctor
===========================
"""

%feature("docstring") GWFrames::Waveform::ZParityViolationMinimized """
Measure the relative magnitude of the violation of parity in the z direction.
=============================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    This function measures the violation of invariance under z-parity
    (reflection across the x-y plane). Nonprecessing systems in a suitable
    frame should have zero violation. Precessing systems in any frame and
    nonprecessing systems in the wrong frame will show violations. The quantity
    is normalized by the overall norm of the data at each instant, and the
    square-root of that ratio is taken.
    
    This is performed iteratively at each time step, as the system is rotated,
    and the parity violation is minimized.
  
"""

%feature("docstring") GWFrames::Waveform::ResizeData """


  Parameters
  ----------
    const unsigned int NModes
    const unsigned int NTimes
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") std """
namespace std
=============
  STL namespace.
  
"""

%feature("docstring") zero """


  Parameters
  ----------
    0. 0
    0. 0
  
  Returns
  -------
    const complex<double>
  
"""

%feature("docstring") GWFrames::DataGrid::Data """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<complex<double>>
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::InnerProduct """


  Parameters
  ----------
    const WaveformAtAPointFT& B
    const vector<double>& InversePSD
  
  Returns
  -------
    double
  
"""

%feature("docstring") virtualWaveformUtilities::Interpolator::~Interpolator """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    ~Interpolator
  
"""

%feature("docstring") IniLIGO_Approx """


  Parameters
  ----------
    const vector<double>& F
    const bool Invert = false
    const double NoiseFloor = 0.0
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::ZeroAbove """


  Parameters
  ----------
    const double Frequency
  
  Returns
  -------
    WaveformAtAPointFT&
  
"""

%feature("docstring") GWFrames::Modes::EllMax """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::Scri::NTimes """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::SliceModes::SliceModes """


  Parameters
  ----------
    const SliceModes& S
  
  Returns
  -------
    SliceModes
  
"""

%feature("docstring") DataGrid::pow """


  Parameters
  ----------
    const int p
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::DataGrid::N_theta """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::Waveform::NormalizedAntisymmetry """
Return the normalized asymmetry as a function of time.
======================================================
  Parameters
  ----------
    vector<int> LModesForAsymmetry = vector<int>(0)
      $\\ell$ modes to use when calculating numerator
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    This function just returns the value of the normalized asymmetry $a$
    defined by Boyle et al. (2014), which is the difference between the field
    at a point and its conjugate at the antipodal point, the amplitude squared
    and integrated over the sphere. The normalization is just the Norm of the
    waveformits overall power at each instant.
    
    By default, all ell modes in the data are used for both the asymmetry and
    the normalization factor. If an argument is input, only modes with ell
    values in that argument will be used to calculate the asymmetry. All modes
    will always be used to calculate the normalization.
  
"""

%feature("docstring") GWFrames::Waveform::AlignDecompositionFrameToModes """
Fix the orientation of the corotating frame.
============================================
  Parameters
  ----------
    const double t_fid
      Fiducial time at which the alignment should happen
    const Quaternions::Quaternion& nHat_t_fid = Quaternions::xHat
      The approximate direction of nHat at t_fid
    const vector<int>& Lmodes = vector<int>(0)
      Lmodes to use in computing $<LL>$
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    The corotating frame is only defined up to some constant rotor R_eps; if
    R_corot is corotating, then so is R_corot*R_eps. This function uses that
    freedom to ensure that the frame is aligned with the Waveform modes at the
    fiducial time. In particular, it ensures that the Z axis of the frame in
    which the decomposition is done is along the dominant eigenvector of $<LL>$
    (suggested by O'Shaughnessy et al.), and the phase of the (2,2) mode is
    zero.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") GWFrames::Waveform::Output """
Output Waveform object to data file.
====================================
  Parameters
  ----------
    const string& FileName
    const unsigned int precision = 14
  
  Returns
  -------
    const Waveform&
  
"""

%feature("docstring") GWFrames::Waveform::MaxNormTime """


  Parameters
  ----------
    const unsigned int SkipFraction = 4
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Waveform::DataType """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::Waveform::RIsScaledOut """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    bool
  
"""

%feature("docstring") GWFrames::SliceOfScri """
class GWFrames::SliceOfScri
===========================
  Member variables
  ----------------
    D psi0

            This class holds all the necessary objects needed to understand the
      geometry of a given slice of null infinity.    D psi1
    D psi2
    D psi3
    D psi4
    D sigma
    D sigmadot
  
"""

%feature("docstring") GWFrames::SliceModes """
class GWFrames::SliceModes
==========================
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::InversePSD """


  Parameters
  ----------
    const string& Detector = 'AdvLIGO_ZeroDet_HighP'
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") dotproduct """


  Parameters
  ----------
    const double * a
    const double * b
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Waveform::SetBoostWeight """


  Parameters
  ----------
    const int NewBoostWeight
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") SplineInterpolator::rawinterp """


  Parameters
  ----------
    int jl
    double xv
  
  Returns
  -------
    double
  
"""

%feature("docstring") SliceModes::SliceModes """
Constructor from ellMax.
========================
  Parameters
  ----------
    const int ellMax = 0
  
  Returns
  -------
    SliceModes
  
"""

%feature("docstring") GWFrames::Waveform::FindModeIndexWithoutError """
Find index of mode with given (l,m) data without the chance of throwing an exception.
=====================================================================================
  Parameters
  ----------
    const int L
    const int M
  
  Returns
  -------
    unsigned int
  
  Description
  -----------
    If the requested mode is not present, the returned index is 1 beyond the
    end of the mode vector.
  
"""

%feature("docstring") GWFrames::PNWaveform::chi2 """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    const vector<double>&
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<vector<double>>&
  
"""

%feature("docstring") GWFrames::PNWaveform::chi1 """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    const vector<double>&
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<vector<double>>&
  
"""

%feature("docstring") GWFrames::operator/ """


  Parameters
  ----------
    const double& a
    const DataGrid& b
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::operator- """


  Parameters
  ----------
    const double& a
    const DataGrid& b
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::Waveform::EvaluateAtPoint """
Evaluate Waveform at a particular sky location.
===============================================
  Parameters
  ----------
    const double vartheta
      Polar angle of detector
    const double varphi
      Azimuthal angle of detector
  
  Returns
  -------
    vector<complex<double>>
  
  Description
  -----------
    Note that the input angle parameters are measured relative to the binary's
    coordinate system. In particular, this will make no sense if the frame type
    is something other than inertial, and will fail if the FrameType is neither
    UnknownFrameType nor Inertial.
  

Evaluate Waveform at a particular sky location and an instant of time.
======================================================================
  Parameters
  ----------
    const double vartheta
      Polar angle of detector
    const double varphi
      Azimuthal angle of detector
    const unsigned int i_t
      Index of time at which to evaluate
  
  Returns
  -------
    complex<double>
  
  Description
  -----------
    Note that the input angle parameters are measured relative to the binary's
    coordinate system. In particular, this will make no sense if the frame type
    is something other than inertial, and will fail if the FrameType is neither
    UnknownFrameType nor Inertial.
  
"""

%feature("docstring") GWFrames::PNWaveform::Omega_totMag """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    double
  

Vector of magnitudes of Omega_tot at each instant of time.
==========================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::DataGrid::Spin """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") WaveformUtilities::convlv """


  Parameters
  ----------
    const vector<double>& data
    const vector<double>& respns
    const int isign
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::CorotatingFrame """
Frame in which the rotation is minimal.
=======================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<Quaternions::Quaternion>
  
  Description
  -----------
    This function combines the steps required to obtain the corotating frame.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") GWFrames::DataGrid::SetSpin """


  Parameters
  ----------
    const int ess
  
  Returns
  -------
    DataGrid&
  
"""

%feature("docstring") ImaginaryI """


  Parameters
  ----------
    0. 0
    1. 0
  
  Returns
  -------
    const complex<double>
  
"""

%feature("docstring") GWFrames::Waveform::DropEllModes """
Remove data relating to the given ell modes.
============================================
  Parameters
  ----------
    const vector<unsigned int>& EllModesToDrop
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") WaveformUtilities::Interpolator::rawinterp """


  Parameters
  ----------
    int jlo
    double x
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::DataGrid::SetNTheta """


  Parameters
  ----------
    const int N_theta
  
  Returns
  -------
    DataGrid&
  
"""

%feature("docstring") GWFrames::Waveform::ZParityViolationSquared """
Measure the absolute magnitude of the violation of parity in the z direction.
=============================================================================
  Parameters
  ----------
    vector<int> Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    This function measures the violation of invariance under z-parity
    (reflection across the x-y plane). Nonprecessing systems in a suitable
    frame should have zero violation. Precessing systems in any frame and
    nonprecessing systems in the wrong frame will show violations. This
    quantity can be minimized over orientation to show the presence or absence
    of a plane of symmetry.
  
"""

%feature("docstring") WaveformUtilities::PolynomialInterpolator::PolynomialInterpolator """


  Parameters
  ----------
    const vector<double>& xv
    const vector<double>& yv
    int m
  
  Returns
  -------
    PolynomialInterpolator
  
"""

%feature("docstring") WaveformUtilities::Interpolator::interp """


  Parameters
  ----------
    double x
  
  Returns
  -------
    double
  
"""

%feature("docstring") DataGrid::DataGrid """


  Parameters
  ----------
    const int Spin
    const int N_theta
    const int N_phi
    const vector<complex<double>>& D
  
  Returns
  -------
    DataGrid
  



  Parameters
  ----------
    Modes M
    const int N_theta = 0
    const int N_phi = 0
  
  Returns
  -------
    DataGrid
  



  Parameters
  ----------
    const Modes& M
    const ThreeVector& v
    const int N_theta = 0
    const int N_phi = 0
  
  Returns
  -------
    DataGrid
  

Constructor on boosted grid by means of functor.
================================================
  Parameters
  ----------
    const int Spin
      Integer spin weight
    const int N_theta
      Number of points in output grid in theta
    const int N_phi
      Number of points in output grid in phi
    const ThreeVector& v
      Three-vector velocity of boosted frame relative to current frame
    const ScriFunctor& f
      Functor operating on a Quaternion object
  
  Returns
  -------
    DataGrid
  
  Description
  -----------
    The functor takes a Quaternion argument, which describes the location and
    orientation of the point to be evaluated. In particular, the rotor takes
    the $\\hat{z}$ vector into the point at which the field is to be measured,
    and takes $\\hat{x} + i \\hat{y}$ into the $m$ vector (within
    normalization) needed for spin-weighted fields.
  
"""

%feature("docstring") GWFrames::Modes::pow """


  Parameters
  ----------
    const int p
  
  Returns
  -------
    Modes
  
"""

%feature("docstring") GWFrames::PNWaveform::OmegaHat_prec """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::Waveform::DataTypeLaTeXString """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    string
  
"""

%feature("docstring") GWFrames::Waveform::SliceOfTimeIndicesWithoutModes """
Copy of the Waveform between indices i_t_a and i_t_b without mode data.
=======================================================================
  Parameters
  ----------
    const unsigned int i_t_a
    unsigned int i_t_b = 0
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    i_t_a and i_t_b should hold the indices pointing to the first time in t
    after t_a, and the first time in t after t_b (or one-past-the-end of t if
    necessary)
  
"""

%feature("docstring") GWFrames::Waveform::DipoleMoment """
Evaluate the dipole moment of the waveform.
===========================================
  Parameters
  ----------
    int ellMax = 0
      Maximum ell mode to include [default: all]
  
  Returns
  -------
    vector<vector<double>>
  
  Description
  -----------
    This function evaluates the dipole moment of the waveform's magnitude,
    defined as $\\vec{d} = \\int \\hat{n} \\lvert f \\rvert^2 d\\Omega$. Up to
    a geometric factor, this function applied to $\\dot{h}$ is the rate of
    emission of momentum in gravitational waves.
  
"""

%feature("docstring") GWFrames::Waveform::BinaryOp """
Pointwise multiply this object by another Waveform object.
==========================================================
  Parameters
  ----------
    typename Op 
    const Waveform& B
  
  Returns
  -------
    typename Op
  
"""

%feature("docstring") GWFrames::PNWaveform::Omega_prec """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    const vector<double>&
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<vector<double>>&
  
"""

%feature("docstring") GWFrames::PNWaveform::L """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    const vector<double>&
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<vector<double>>&
  
"""

%feature("docstring") GWFrames::Modes::operator[] """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    complex<double>
  



  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    complex<double>&
  
"""

%feature("docstring") GWFrames::DataGrid::size """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") Modes::EvaluateAtPoint """
Evaluate Waveform at a particular sky location.
===============================================
  Parameters
  ----------
    const double vartheta
      Polar angle of detector
    const double varphi
      Azimuthal angle of detector
  
  Returns
  -------
    complex<double>
  

Evaluate Waveform at a particular sky location.
===============================================
  Parameters
  ----------
    const Quaternions::Quaternion& R
      Quaternion giving point by rotation of $\\hat{z}$
  
  Returns
  -------
    complex<double>
  
  Description
  -----------
    Note that the argument R might typically be thought of as the rotor taking
    the unit $z$ vector into a point $(\\vartheta, \\varphi)$. However, more
    general arguments are possible; this feature is used to greatly simplify
    the process of constructing a DataGrid object from a Modes object with a
    boost.
  
"""

%feature("docstring") GWFrames::Waveform::SliceOfTimes """
Copy the Waveform between t_a and t_b.
======================================
  Parameters
  ----------
    const double t_a = -1e300
    const double t_b = 1e300
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") GWFrames::Waveform::LLMatrix """
Calculate the $<LL>$ quantity defined in the paper.
===================================================
  Parameters
  ----------
    vector<int> Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<Matrix>
  
  Description
  -----------
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
    
    $<LL>^{ab} = \\sum_{\\ell,m,m'} [\\bar{f}^{\\ell,m'} <\\ell,m' | L_a L_b |
    \\ell,m> f^{\\ell,m} ]$
  
"""

%feature("docstring") GWFrames::PNWaveform::LHat """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::Waveform::Norm """
Return the norm (sum of squares of modes) of the waveform.
==========================================================
  Parameters
  ----------
    const bool TakeSquareRoot = false
      If true, the square root is taken at each instant before returning
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    This returns the norm of the waveform, defined as the sum of the complex
    norms of the modes. Note that we are calling this norm in analogy with the
    c++ std::complex norm, which is the square of the absolute value. However,
    there is also an option to take the square root of the data at each time
    step, which would be the usual L2 norm of the waveform.
    
    MaxNormIndex
    
    MaxNormTime
  
"""

%feature("docstring") GWFrames::WaveformAtAPointFT::Normalize """


  Parameters
  ----------
    const vector<double>& InversePSD
  
  Returns
  -------
    WaveformAtAPointFT&
  



  Parameters
  ----------
    const string& Detector = 'AdvLIGO_ZeroDet_HighP'
  
  Returns
  -------
    WaveformAtAPointFT&
  
"""

%feature("docstring") GWFrames::Waveform::LLDominantEigenvector """
Calculate the principal axis of the LL matrix, as prescribed by O'Shaughnessy et al.
====================================================================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<vector<double>>
  
  Description
  -----------
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") GWFrames::DataGrid::DataGrid """


  Parameters
  ----------
    const int size = 0
  
  Returns
  -------
    DataGrid
  



  Parameters
  ----------
    const DataGrid& A
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::Waveform::MIsScaledOut """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    bool
  
"""

%feature("docstring") GWFrames::PNWaveform::~PNWaveform """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    ~PNWaveform
  
"""

%feature("docstring") WaveformUtilities::SplineInterpolator """
class WaveformUtilities::SplineInterpolator
===========================================
  Member variables
  ----------------
    vector<double> y2
  
"""

%feature("docstring") GWFrames::Waveform::KeepOnlyEllModes """
Remove data relating to all but the given ell modes.
====================================================
  Parameters
  ----------
    const vector<unsigned int>& EllModesToKeep
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") WaveformUtilities::WrapVecDoub::WrapVecDoub """


  Parameters
  ----------
    const int nn
  
  Returns
  -------
    WrapVecDoub
  



  Parameters
  ----------
    vector<double>& vec
  
  Returns
  -------
    WrapVecDoub
  
"""

%feature("docstring") GWFrames::Waveform::SetSpinWeight """


  Parameters
  ----------
    const int NewSpinWeight
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::Waveform::AppendHistory """


  Parameters
  ----------
    const string& Hist
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::Modes::Data """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<complex<double>>
  
"""

%feature("docstring") GWFrames::Waveform::Compare """
Return a Waveform with differences between the two inputs.
==========================================================
  Parameters
  ----------
    const Waveform& B
    const double MinTimeStep = 0.005
    const double MinTime = -3.0e300
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This function simply subtracts the data in this Waveform from the data in
    Waveform A, and finds the rotation needed to take this frame into frame A.
    Note that the waveform data are stored as complex numbers, rather than as
    modulus and phase.
  
"""

%feature("docstring") GWFrames::Waveform::Im """


  Parameters
  ----------
    const unsigned int Mode
    const unsigned int TimeIndex
  
  Returns
  -------
    double
  

Return vector of imaginary parts of a given mode as function of time.
=====================================================================
  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    vector<double>
  

Return vector of vector of imaginary parts of all modes as function of time.
============================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::SuperMomenta::T """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::GetAlignmentOfDecompositionFrameToModes """
Find the appropriate rotation to fix the orientation of the corotating frame.
=============================================================================
  Parameters
  ----------
    const double t_fid
      Fiducial time at which the alignment should happen
    const Quaternions::Quaternion& nHat_t_fid = Quaternions::xHat
      The approximate direction of nHat at t_fid
    const vector<int>& Lmodes = vector<int>(0)
      Lmodes to use in computing $<LL>$
  
  Returns
  -------
    Quaternions::Quaternion
  
  Description
  -----------
    This function simply finds the rotation necessary to align the corotating
    frame to the waveform at the fiducial time, rather than applying it. This
    is called by AlignDecompositionFrameToModes and probably does not need to
    be called directly; see that function's documentation for more details.
    
    AlignDecompositionFrameToModes
  
"""

%feature("docstring") GWFrames::Waveform::SetFrame """


  Parameters
  ----------
    const vector<Quaternions::Quaternion>& a
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") WaveformUtilities::PolynomialInterpolator """
class WaveformUtilities::PolynomialInterpolator
===============================================
  Member variables
  ----------------
    double dy
  
"""

%feature("docstring") SliceModes::FourMomentum """
Calculate the four-momentum of the system from the supermomentum.
=================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    FourVector
  
  Description
  -----------
    The (Bondi) four-momentum is given by the ell=0 and ell=1 modes of the
    supermomentum.
  
"""

%feature("docstring") GWFrames::Waveform::FrameType """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") ComplexI """


  Parameters
  ----------
    0. 0
    1. 0
  
  Returns
  -------
    const complex<double>
  
"""

%feature("docstring") GWFrames::PNWaveform::chi1Mag """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::DataGrid """
class GWFrames::DataGrid
========================
  Member variables
  ----------------
    int s

            This object holds complex spin-weighted data on the sphere in an
      equi-angular representation. That is, given integers n_theta and n_phi,
      the data are recorded on a 'rectangular' grid with n_theta*n_phi points,
      including n_phi points at each of the poles. The purpose of these objects
      is primarily to serve as a computational tool for pointwise
      multiplication. Otherwise, Modes is expected to be the preferable
      representation.    int n_theta
    int n_phi
    vector<complex<double>> data
  
"""

%feature("docstring") GWFrames::Waveform::SetLM """


  Parameters
  ----------
    const vector<vector<int>>& a
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::Waveform::RotateDecompositionBasis """
Rotate the basis in which this Waveform is measured by a constant rotor.
========================================================================
  Parameters
  ----------
    const Quaternions::Quaternion& R_frame
  
  Returns
  -------
    Waveform&
  

Rotate the basis in which this Waveform is measured.
====================================================
  Parameters
  ----------
    const vector<Quaternions::Quaternion>& R_frame
      Vector of Quaternions by which to rotate
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    This rotates the coordinate basis, leaving the physical system in place.
    
    The Waveform's frame data records the rotors needed to rotate the standard
    (x,y,z) basis into the (X,Y,Z) basis with respect to which the Waveform
    modes are decomposed. If this is not the first rotation of the frame, we
    need to be careful about how we record the total rotation. Here, we are
    just composing rotations, so we need to store R_frame times the original
    frame data.
    
    Note that this function does not change the frameType; this is left to the
    calling function.
  
"""

%feature("docstring") GWFrames::InverseConformalFactorBoostedGrid """
Construct a boosted grid with the conformal factor at each point.
=================================================================
  Parameters
  ----------
    const ThreeVector& v
    const int n_theta
    const int n_phi
  
  Returns
  -------
    DataGrid
  
"""

%feature("docstring") GWFrames::Waveform::GetAlignmentsOfDecompositionFrameToModes """
Find the appropriate rotations to fix the orientation of the corotating frame.
==============================================================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      Lmodes to use in computing $<LL>$
  
  Returns
  -------
    vector<Quaternions::Quaternion>
  
  Description
  -----------
    This function finds the appropriate pre-multiplied rotation
    $R_{\\varepsilon}$ so that the decomposition frame is aligned to the
    waveform. This particular version finds the appropriate $R_{\\varepsilon}$
    at each time in the input Waveform. This is useful in cases where we need
    to try many such alignments, because the setup for interpolation is very
    slow.
    
    Note that this function has no option to choose the direction of X based on
    some nHat vector, as other similar functions have. That issue is assumed to
    be handled elsewhere.
  
"""

%feature("docstring") GWFrames::Waveform::CopyWithoutData """
Copy the Waveform, except for the data (t, frame, lm, data)
===========================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") WaveformUtilities::SplineInterpolator::SplineInterpolator """


  Parameters
  ----------
    const vector<double>& xv
    const vector<double>& yv
    double yp1 = 1.e99
    double ypn = 1.e99
  
  Returns
  -------
    SplineInterpolator
  
"""

