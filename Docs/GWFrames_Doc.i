%feature("docstring") GWFrames::Waveform::TransformUncertaintiesToInertialFrame """
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

%feature("docstring") GWFrames::SWSH::spin """


  Parameters
  ----------
    s 
  
  Returns
  -------
    spin
  
"""

%feature("docstring") GWFrames::Quaternion """
class GWFrames::Quaternion
==========================
  Object representing an individual quaternion.
  
  Member variables
  ----------------
    double w
    double x
    double y
    double z
  
"""

%feature("docstring") func """


  Parameters
  ----------
    double t
    const double y
    double dydt
    void * params
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::FrameFromPrescribedRotation """
Construct minimal-rotation frame from Z basis vector of that frame.
===================================================================
  Parameters
  ----------
    const vector<Quaternion>& omega
      Vector of Quaternions
    const vector<double>& T
      Vector of corresponding times
    const unsigned int NIterations = 5
      Number of refinements [default: 5]
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    The input vector of Quaternions represent the angular-velocity vector
    (omega) of the frame at each instant of time. The returned vector of rotors
    will rotate the stationary frame's (x,y,z) vectors into the new frame's
    (X,Y,Z) vectors, where Z is parallel to omega, and the X and Y vectors are
    deduced by enforcing the condition that the instantaneous rotation of the
    frame about Z is |omega|. Note that this leaves an unfixed initial rotation
    in the XY plane.
  
"""

%feature("docstring") MatrixC::MatrixC """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    MatrixC
  



  Parameters
  ----------
    int n
    int m
  
  Returns
  -------
    MatrixC
  



  Parameters
  ----------
    int n
    int m
    const complex<double>& a
  
  Returns
  -------
    MatrixC
  



  Parameters
  ----------
    int n
    int m
    const complex<double> * a
  
  Returns
  -------
    MatrixC
  



  Parameters
  ----------
    const vector<vector<complex<double>>>& DataIn
  
  Returns
  -------
    MatrixC
  



  Parameters
  ----------
    const MatrixC& rhs
  
  Returns
  -------
    MatrixC
  
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

%feature("docstring") MatrixC::assign """


  Parameters
  ----------
    int newn
    int newm
    const complex<double>& a
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::SpinWeight """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
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

%feature("docstring") GWFrames::Quaternion::Quaternion """
Empty constructor  initialized to 0s.
=====================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  

Copy constructor.
=================
  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  

Constructor from spherical coordinates.
=======================================
  Parameters
  ----------
    const double vartheta
      Float representing the polar angle
    const double varphi
      Float representing the azimuthal angle
  
  Returns
  -------
    Quaternion
  
  Description
  -----------
    The unit Quaternion constructed in this way rotates the z axis onto the
    point given by the coordinates (vartheta, varphi).
  

Constructor from Euler angles.
==============================
  Parameters
  ----------
    const double alpha
      First Euler angle
    const double beta
      Second Euler angle
    const double gamma
      Third Euler angle
  
  Returns
  -------
    Quaternion
  
  Description
  -----------
    The unit Quaternion constructed in this way corresponds to a rotation by
    the given Euler angles. The convention used here is the z-y-z convention.
    That is, the rotations occur about the fixed axes: first a rotation by
    gamma about the z axis, then a rotation by beta about the y axis, and
    finally a rotation by alpha about the z axis.
  

Constructor by components.
==========================
  Parameters
  ----------
    const double w0
      Scalar component of Quaternion
    const double x0
      First vector component of Quaternion
    const double y0
      Second vector component of Quaternion
    const double z0
      Third vector component of Quaternion
  
  Returns
  -------
    Quaternion
  

Constructor from vector.
========================
  Parameters
  ----------
    const vector<double>& q
      Vector containing three or four components
  
  Returns
  -------
    Quaternion
  
  Description
  -----------
    If the input vector has three components, they are assumed to represent the
    vector components of the Quaternion, and the scalar component is set to
    zero. If the input vector has four components, they are assumed to
    represent the four components of the Quaternion, with the 0 component being
    the scalar part.
  

Constructor from axis-angle.
============================
  Parameters
  ----------
    const double angle
      Single number giving the rotation angle
    const vector<double>& axis
      Three-component vector (assumed to be normalized) giving the axis
  
  Returns
  -------
    Quaternion
  
  Description
  -----------
    This constructs a rotor (assuming 'axis' is normalized) corresponding to
    rotation about the given axis through the given angle.
  
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

%feature("docstring") GWFrames::pow """


  Parameters
  ----------
    const Quaternion& Q
    const double x
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const Quaternion& Q
    const Quaternion& P
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const double x
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const Quaternion& P
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& Q
    const vector<double>& x
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& Q
    const vector<Quaternion>& P
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const vector<double>& x
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const vector<Quaternion>& P
  
  Returns
  -------
    vector<Quaternion>
  
"""

%feature("docstring") GWFrames::Waveform::HackSpECSignError """
Correct the error in RWZ extraction from older SpEC files.
==========================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform&
  
"""

%feature("docstring") GWFrames::Waveform::PNEquivalentOrbitalAV """
Deduce PN-equivalent orbital angular velocity from Waveform.
============================================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<vector<double>>
  
  Description
  -----------
    This function simply takes the projection of the field's angular-velocity
    vector $\\vec{\\omega}$ along the dominant eigenvector $\\hat{V}_f$ of
    $<LL>$. This should be equivalent to the orbital angular velocity of the PN
    system. Note that the returned vector is relative to the inertial frame.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") GWFrames::Waveform::SchmidtEtAlVector """


  Parameters
  ----------
    const double alpha0Guess = 0.0
    const double beta0Guess = 0.0
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::SWSH::SetRotation """


  Parameters
  ----------
    const Quaternion& iR
  
  Returns
  -------
    SWSH&
  
"""

%feature("docstring") Matrix::Matrix """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Matrix
  



  Parameters
  ----------
    unsigned int rows
    unsigned int cols
  
  Returns
  -------
    Matrix
  



  Parameters
  ----------
    unsigned int rows
    unsigned int cols
    const double a
  
  Returns
  -------
    Matrix
  



  Parameters
  ----------
    const Matrix& rhs
  
  Returns
  -------
    Matrix
  



  Parameters
  ----------
    const vector<vector<double>>& DataIn
  
  Returns
  -------
    Matrix
  
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
    Quaternion
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    const vector<Quaternion>&
  
"""

%feature("docstring") GWFrames::MatrixC::nrows """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
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

%feature("docstring") GWFrames::Waveform::SetTime """


  Parameters
  ----------
    const vector<double>& a
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::ComplexDerivative """
Three-point finite-differencing of vector of complex<double>.
=============================================================
  Parameters
  ----------
    const vector<complex<double>>& f
      Vector of complex<double>.
    const vector<double>& t
      Vector of corresponding time steps.
  
  Returns
  -------
    vector<complex<double>>
  
  Description
  -----------
    Sundquist and Veronis, Tellus XXII (1970), 1
  
"""

%feature("docstring") GWFrames::Waveform::SetRIsScaledOut """


  Parameters
  ----------
    const bool Scaled
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::QuaternionDerivative """
Three-point finite-differencing of vector of Quaternions.
=========================================================
  Parameters
  ----------
    const vector<Quaternion>& f
      Vector of Quaternions.
    const vector<double>& t
      Vector of corresponding time steps.
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    Sundquist and Veronis, Tellus XXII (1970), 1
  
"""

%feature("docstring") GWFrames::Waveform::RotatePhysicalSystem """
Rotate the physical content of the Waveform by a constant rotor.
================================================================
  Parameters
  ----------
    const Quaternion& R_phys
  
  Returns
  -------
    Waveform&
  

Rotate the physical content of the Waveform.
============================================
  Parameters
  ----------
    vector<Quaternion> R_phys
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

%feature("docstring") GWFrames::WignerCoefficientFunctor::WignerCoefficientFunctor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    WignerCoefficientFunctor
  
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

%feature("docstring") SQR """


  Parameters
  ----------
    const double& x
  
  Returns
  -------
    double
  



  Parameters
  ----------
    const double a
  
  Returns
  -------
    double
  



  Parameters
  ----------
    const double a
  
  Returns
  -------
    double
  
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

%feature("docstring") GWFrames::Waveform::GHPEdthBar """
Geroch-Held-Penrose edth operator conjugate.
============================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This operator is the one defined by Geroch et al. (1973). It raises the
    spin weight of any field on the sphere by 1, while leaving the boost weight
    unchanged.
    
    This operator is very similar to the basic Newman-Penrose edth operator,
    except that it preserves boost weights. In certain cases, it reduces to the
    NP edth (up to a factor of $\\sqrt{2}$), but not generally. In those
    certain cases, we have NPEdthBar() = sqrt(2)*GHPEdthBar().
    
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
    
    NPEdth
    
    NPEdthBar
    
    GHPEdth
  
"""

%feature("docstring") WignerDMatrix::SetRotation """
Reset the rotor for this object to the given value.
===================================================
  Parameters
  ----------
    const Quaternion& iR
  
  Returns
  -------
    WignerDMatrix&
  
"""

%feature("docstring") GWFrames::Waveform::NPEdth """
Newman-Penrose edth operator.
=============================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This operator is the one defined by Newman and Penrose (1966) and further
    described by Goldberg et al. (1967). It raises the spin weight of any field
    on the sphere by 1. Note that this operator does not preserve boost weights
    in any nice way  except in special cases. The GHP version does.
    
    Note that the boost weight is set to the value of WeightError, which is
    just meant to be large enough that it will give improbable values if used.
    This is not fool-proof.
    
    NPEdthBar
    
    GHPEdth
    
    GHPEdthBar
  
"""

%feature("docstring") GWFrames::Eigenvectors """


  Parameters
  ----------
    Matrix& M
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::TransformToOShaughnessyEtAlFrame """
Transform Waveform to O'Shaughnessy et al. frame.
=================================================
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
    O'Shaughnessy et al. frame.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") GWFrames::BinomialCoefficientFunctor """
class GWFrames::BinomialCoefficientFunctor
==========================================
  Object for pre-computing and retrieving binomials.
  
  Member variables
  ----------------
    const vector<double> BinomialCoefficientTable
  
"""

%feature("docstring") GWFrames::Waveform::SetLM """


  Parameters
  ----------
    const vector<vector<int>>& a
  
  Returns
  -------
    void
  
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

%feature("docstring") GWFrames::Waveform::AlignDecompositionFrameToModes """
Fix the orientation of the corotating frame.
============================================
  Parameters
  ----------
    const double t_fid
      Fiducial time at which the alignment should happen
    const vector<int>& Lmodes = vector<int>(0)
      Lmodes to use in computing $<LL>$
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    The corotating frame is only defined up to some constant rotor R_c; if
    R_corot is corotating, then so is R_corot*R_c. This function uses that
    freedom to ensure that the frame is aligned with the Waveform modes at the
    fiducial time. In particular, it ensures that the Z axis of the frame in
    which the decomposition is done is along the dominant eigenvector of $<LL>$
    (suggested by O'Shaughnessy et al.), and the phase of the (2,2) mode is
    zero.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
"""

%feature("docstring") GWFrames::Quaternions """


  Parameters
  ----------
    const vector<double>& vartheta
    const vector<double>& varphi
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& alpha
    const vector<double>& beta
    const vector<double>& gamma
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& w0
    const vector<double>& x0
    const vector<double>& y0
    const vector<double>& z0
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<vector<double>>& q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& angle
    const vector<vector<double>>& axis
  
  Returns
  -------
    vector<Quaternion>
  
"""

%feature("docstring") GWFrames::Intersection """
Return the intersection of two time sequences.
==============================================
  Parameters
  ----------
    const vector<double>& t1
    const vector<double>& t2
    const double MinStep = 0.005
    const double MinTime = -1e300
    const double MaxTime = 1e300
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    The time step at each point is the minimum of the time steps in t1 and t2
    at that instant, or MinStep, whichever is greater. The output starts at the
    earliest moment common to t1 and t2, or MinTime, whichever is greater.
    
    The input to this function is assumed to be strictly monotonic.
  
"""

%feature("docstring") GWFrames::Eigenvalues """


  Parameters
  ----------
    Matrix& M
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::WignerCoefficientFunctor """
class GWFrames::WignerCoefficientFunctor
========================================
  Object for pre-computing and retrieving coefficients for the Wigner D
  matrices.
  
  Member variables
  ----------------
    const vector<double> CoefficientTable
  
"""

%feature("docstring") GWFrames::Waveform::TransformToSchmidtEtAlFrame """
Transform Waveform to Schmidt et al. frame.
===========================================
  Parameters
  ----------
    const double alpha0Guess = 0.0
      Initial guess for optimal direction alpha
    const double beta0Guess = 0.0
      Initial guess for optimal direction beta
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    This function combines the steps required to obtain the Waveform in the
    Schmidt et al. frame.
  
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
    vector<Quaternion> frame
    WaveformFrameType frameType
    WaveformDataType dataType
    bool rIsScaledOut
    bool mIsScaledOut
    vector<vector<int>> lm
    MatrixC data
  
  Non-public member functions
  ---------------------------
    Waveform& TransformModesToRotatedFrame
      Rotate modes of the Waveform object.    Waveform& TransformUncertaintiesToRotatedFrame
      Rotate modes of the uncertainty of a Waveform object.  
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
    void
  
"""

%feature("docstring") GWFrames::Quaternion::operator+ """


  Parameters
  ----------
    const double t
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Quaternion::operator* """


  Parameters
  ----------
    const double t
  
  Returns
  -------
    Quaternion
  

Quaternion multiplication.
==========================
  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Matrix::gslobj """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    gsl_matrix *
  
"""

%feature("docstring") GWFrames::Quaternion::operator/ """


  Parameters
  ----------
    const double t
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  
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

%feature("docstring") GWFrames::Quaternion::operator- """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const double t
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  
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

%feature("docstring") GWFrames::FrameFromAngularVelocity """


  Parameters
  ----------
    const vector<Quaternion>& Omega
      Vector of Quaternions.
    const vector<double>& T
      Vector of corresponding times.
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    Note that each element of Omega should be a pure-vector Quaternion,
    corresponding to the angular-velocity vector at the instant of time.
  
"""

%feature("docstring") GWFrames::Waveforms::operator[] """


  Parameters
  ----------
    const int i
  
  Returns
  -------
    const Waveform&
  



  Parameters
  ----------
    const int i
  
  Returns
  -------
    Waveform&
  
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
    done already. The transition function is a smooth
    
    Note that this function does NOT operate in place; a new Waveform object is
    constructed and returned.
  
"""

%feature("docstring") GWFrames::Matrix::~Matrix """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    ~Matrix
  
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

%feature("docstring") VectorStringForm """


  Parameters
  ----------
    const vector<double>& V
  
  Returns
  -------
    string
  
"""

%feature("docstring") GWFrames::Quaternion::operator= """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion&
  
"""

%feature("docstring") GWFrames::exp """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  
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

%feature("docstring") GWFrames::Waveforms::clear """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::SetDataType """


  Parameters
  ----------
    const WaveformDataType Type
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::FrameFromXY """
Construct frame given the X and Y basis vectors of that frame.
==============================================================
  Parameters
  ----------
    const vector<Quaternion>& X
      Vector of Quaternions
    const vector<Quaternion>& Y
      Vector of Quaternions
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    The input parameters are Quaternions, assumed to be pure unit vectors,
    representing the X and Y basis vectors of the frame at each instant of
    time. The returned vector of rotors will rotate the stationary frame's
    (x,y,z) vectors into the new frame's (X,Y,Z) vectors.
  
"""

%feature("docstring") GWFrames::Component """


  Parameters
  ----------
    const vector<Quaternion>& Q
    const unsigned int i
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") FactorialTableCalculator """
Class to create an object returning the factorial of an argument.
=================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    Note that because a double is returned, only values up to 28! will be
    exact; higher values will be accurate to machine precision. Values up to
    170! only are allowed because higher values overflow.
  
"""

%feature("docstring") StringForm """


  Parameters
  ----------
    const vector<int>& Lmodes
  
  Returns
  -------
    string
  
"""

%feature("docstring") GWFrames::inverse """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  
"""

%feature("docstring") GWFrames::WignerDMatrix """
class GWFrames::WignerDMatrix
=============================
  Object for computing the Wigner D matrices as functions of quaternion rotors.
  
  Member variables
  ----------------
    BinomialCoefficientFunctor BinomialCoefficient
    WignerCoefficientFunctor WignerCoefficient
    complex<double> Ra
    complex<double> Rb
    double absRa
    double absRb
    double absRRatioSquared
  
"""

%feature("docstring") GWFrames::Quaternion::pow """


  Parameters
  ----------
    const double t
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") Matrix::swap """


  Parameters
  ----------
    Matrix& b
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::GetAlignmentOfFrame """
Get the rotor needed to align this waveform's frame to the other's at the given time.
=====================================================================================
  Parameters
  ----------
    const Waveform& A
      Fixed Waveform in corotating frame to which this Waveform is aligned
    const double t_fid
      Fiducial time at which to equate frames
    Quaternion& R_delta
      Returned rotor
  
  Returns
  -------
    void
  
  Description
  -----------
    This function simply finds the rotation necessary to align this waveform's
    frame to the other at the fiducial time, rather than applying it. This is
    called by AlignFrame and probably does not need to be called directly; see
    that function's documentation for more details.
    
    AlignFrame
  
"""

%feature("docstring") GWFrames::Quaternion::abs """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    double
  
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
    const double delta
      Normalized BH mass difference (M1-M2)/(M1+M2)
    const vector<double>& chi1_0
      Initial dimensionless spin vector of BH1
    const vector<double>& chi2_0
      Initial dimensionless spin vector of BH2
    const double Omega_orb_0
      Initial orbital angular frequency
    const Quaternion& R_0 = Quaternion(1, 0, 0, 0)
      Overall rotation of the system (optional)
  
  Returns
  -------
    PNWaveform
  
  Description
  -----------
    The PN system is initialized having the BHs along the x axis, with the
    orbital angular velocity along the positive z axis, having magnitude
    Omega_orb_0. The input spin vectors must be defined with respect to this
    basis.
    
    The TaylorT1 system is first integrated to compute the dynamics of the
    binary. The evolved spin vectors chi1 and chi2, orbital angular-velocity
    vector Omega_orb, and orbital phase Phi_orb are stored. Simultaneously, the
    minimal-rotation frame of the angular-velocity vector is computed, then
    rotated about the z' axis by Phi_orb, resulting in the binary's frame. Once
    this step is completed, the information is used to construct the waveform
    in the minimal-rotation frame. (That is, the waveform will be essentially
    corotating.)
    
    Note that, to get the PNWaveform in an inertial frame, you must first apply
    the method TransformToCorotatingFrame().
  
"""

%feature("docstring") GWFrames::Waveform::operator() """


  Parameters
  ----------
    const unsigned int Mode
  
  Returns
  -------
    const complex<double> *
  
"""

%feature("docstring") GWFrames::Quaternion::operator!= """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    bool
  
"""

%feature("docstring") GWFrames::CumulativeScalarIntegral """
Integrate scalar function by simple trapezoidal rule.
=====================================================
  Parameters
  ----------
    const vector<double>& fdot
      Vector of scalars.
    const vector<double>& t
      Vector of corresponding time steps.
  
  Returns
  -------
    double
  
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

%feature("docstring") MatrixC::operator= """


  Parameters
  ----------
    const MatrixC& rhs
  
  Returns
  -------
    MatrixC&
  
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

%feature("docstring") GWFrames::Waveform::GetAlignmentOfTimeAndFrame """
Get time and frame offset for alignment over extended region.
=============================================================
  Parameters
  ----------
    const Waveform& A
      Fixed Waveform in corotating frame to which this Waveform is aligned
    const double t1
      Initial time of region over which differences are minimized
    const double t2
      Final time of region over which differences are minimized
    double& deltat
      Returned time offset
    Quaternion& R_delta
      Returned rotation offset
  
  Returns
  -------
    void
  
  Description
  -----------
    This function simply finds the time and rotation shifts necessary to align
    this waveform to the other at the fiducial time, rather than applying it.
    This is called by AlignTimeAndFrame and probably does not need to be called
    directly; see that function's documentation for more details.
    
    AlignTimeAndFrame
  
"""

%feature("docstring") GWFrames::normsquared """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    double
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveforms::SetCommonTime """
Interpolate to a common set of times.
=====================================
  Parameters
  ----------
    vector<vector<double>>& Radii
    const double MinTimeStep = 0.005
    const double EarliestTime = -3e300
    const double LatestTime = 3e300
  
  Returns
  -------
    void
  
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

%feature("docstring") GWFrames::Waveforms::Extrapolate """
Main extrapolation routine.
===========================
  Parameters
  ----------
    vector<vector<double>>& Radii
      Array of radii for each Waveform (first index) and each time (second
      index)
    const vector<int>& ExtrapolationOrders
      List of integers denote extrapolation orders
    const vector<double>& Omegas = vector<double>(0)
      Optional list of angular frequencies for scaling extrapolation polynomial
  
  Returns
  -------
    Waveforms
  
  Description
  -----------
    The input FiniteRadiusWaveforms are assumed to be properly scaled and
    time-retarded, and interpolated to a uniform set of retarded times. This
    function simply steps through the indices, fitting those data to
    polynomials in 1/radius, and evaluating at 0 (for infinity).
    
    The extrapolation orders can be negative. In this case, the scaled,
    time-retarded waveform at finite radius is given, where N=-1 is the
    outermost Waveform, N=-2 is the second to outermost, etc.
    
    Note that the fitting uses gsl_multifit_linear_usvd, which is GSL's fitting
    function that does NOT use column scaling (specified by the 'u' in front of
    'svd' in the function name). The basic GSL fitting function uses column
    scaling 'to improve
the accuracy of the singular values'. However, for
    convergent series, this scaling can make all the coefficients roughly equal
    (just as the Omegas option does), which defeats the SVD.
  
"""

%feature("docstring") GWFrames::Waveforms::~Waveforms """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    ~Waveforms
  
"""

%feature("docstring") GWFrames::BinomialCoefficientFunctor::operator() """


  Parameters
  ----------
    const unsigned int n
    const unsigned int k
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::BinomialCoefficientFunctor::BinomialCoefficientFunctor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    BinomialCoefficientFunctor
  
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

%feature("docstring") GWFrames::Waveform::AlignTime """
Change this Waveform by aligning to the other at the given time.
================================================================
  Parameters
  ----------
    const Waveform& A
      Fixed Waveform in inertial frame to which this Waveform is aligned
    const double t_fid
      Note that this function operates in place; the Waveform to which it is
      applied will change.
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    As noted above, it is implicitly assumed that both Waveforms are in an
    inertial frame, so that the magnitude of the angular velocity may be
    properly measured. This could be adjusted to account for the angular
    velocity of the frame, but hasn't been yet.
    
    To improve accuracy, the angular velocity of A is interpolated to t_fid.
    The time of B is then interpolated to the interpolated angular velocity.
    This assumes that B's angular velocity is strictly monotonic for roughly 5
    data points to either side.
  
"""

%feature("docstring") GWFrames::Waveform::BoostWeight """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::Eigensystem """


  Parameters
  ----------
    Matrix& M
  
  Returns
  -------
    vector<double>
  
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

%feature("docstring") GWFrames::Waveform::AlignTimeAndFrame """
Align time and frame over extended region.
==========================================
  Parameters
  ----------
    const Waveform& A
      Fixed Waveform in corotating frame to which this Waveform is aligned
    const double t1
      Initial time of region over which differences are minimized
    const double t2
      Final time of region over which differences are minimized
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    Note that this function operates in place; the Waveform to which it is
    applied will change. However, the modes are not altered; only the t and
    frame data are.
    
    As noted above, it is implicitly assumed that both Waveforms are in their
    corotating frames, with the modes appropriately aligned to the frames at
    t_fid. The assumption is that the frames actually represent something
    physically meaningful, so that it is meaningful to insist that they be the
    same.
    
    Then, this function adjust the time and orientation of this Waveform, so
    that the difference between the two frames is minimized. That difference is
    measured by finding the rotor R_Delta required to rotate one frame into the
    other, taking the angle of that rotor, and integrating over the region [t1,
    t2].
    
    Relative to the inertial basis, the physical measurables (angular-velocity
    vector and dominant eigenvector of $<LL>$) of this Waveform are rotated.
    
    AlignDecompositionFrameToModes
  
"""

%feature("docstring") MatrixC::~MatrixC """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    ~MatrixC
  
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

%feature("docstring") WignerCoefficientCalculator """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    We need (2*ell+1)*(2*ell+1) coefficients for each value of ell from 0 (for
    completenes) up to ellMax_Utilities (hard coded in the header file). That's
    a total of from sympy import summation, symbols, simplify from
    sympy.polys.polyfuncs import horner ell, ellMax_Utilities, m, mp =
    symbols('ell ellMax_Utilities m mp', integer=True)
    horner(simplify(summation((2*ell+1)**2, (ell, 0, ellMax_Utilities))))
    
    ellMax_Utilities*(ellMax_Utilities*(4*ellMax_Utilities/3 + 4) + 11/3) + 1
    With a similar calculation, we can see that the associated access operator
    needs element horner(summation((2*ell+1)**2, (ell, 0, ell-1)) +
    (2*ell+1)*(ell+mp) + ell + m)
    
    ell*(ell*(4*ell/3 + 2) + 5/3) + mp*(2*ell + 1) + m of the array.
  
"""

%feature("docstring") GWFrames::Waveform::operator/ """


  Parameters
  ----------
    const Waveform& B
  
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

%feature("docstring") GWFrames::FactorialFunctor::FactorialFunctor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    FactorialFunctor
  
"""

%feature("docstring") MatrixC::resize """


  Parameters
  ----------
    int newn
    int newm
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::SetHistory """


  Parameters
  ----------
    const string& Hist
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::SetData """


  Parameters
  ----------
    const vector<vector<complex<double>>>& a
  
  Returns
  -------
    void
  



  Parameters
  ----------
    const unsigned int i_Mode
    const unsigned int i_Time
    const complex<double>& a
  
  Returns
  -------
    void
  
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

%feature("docstring") SchmidtEtAl_f """


  Parameters
  ----------
    const gsl_vector * v
    void * params
  
  Returns
  -------
    double
  
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
    vector<Quaternion>
  
  Description
  -----------
    This function combines the steps required to obtain the corotating frame.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
  
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
    void
  
"""

%feature("docstring") WignerDMatrix::WignerDMatrix """
Construct the D matrix object given the (optional) rotor.
=========================================================
  Parameters
  ----------
    const Quaternion& iR = Quaternion(1, 0, 0, 0)
  
  Returns
  -------
    WignerDMatrix
  
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

%feature("docstring") SchmidtEtAl_df """


  Parameters
  ----------
    const gsl_vector * v
    void * params
    gsl_vector * df
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::~Waveform """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    ~Waveform
  
"""

%feature("docstring") GWFrames::Matrix::gslobj """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const gsl_matrix *
  
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

%feature("docstring") GWFrames::Quaternion::sqrtOfRotor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  
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

%feature("docstring") GWFrames::operator+ """


  Parameters
  ----------
    const double a
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const double a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const Quaternion& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const Quaternion& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const double a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& Q
    const vector<double>& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const vector<double>& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& P
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const vector<double>& b
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const vector<double>& a
    const double b
  
  Returns
  -------
    vector<double>
  
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

%feature("docstring") GWFrames::LadderOperatorFactorFunctor """
class GWFrames::LadderOperatorFactorFunctor
===========================================
  Object for pre-computing and retrieving values of the ladder operators.
  
  Member variables
  ----------------
    const vector<double> FactorTable
  
"""

%feature("docstring") GWFrames::SWSH::SetAngles """


  Parameters
  ----------
    const double vartheta
    const double varphi
  
  Returns
  -------
    SWSH&
  
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

%feature("docstring") GWFrames::angle """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    double
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::operator* """


  Parameters
  ----------
    const double a
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const double a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const Quaternion& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const Quaternion& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const double a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& Q
    const vector<double>& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const vector<double>& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& P
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const Matrix& b
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const Quaternion& a
    const Matrix& b
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Waveform::operator+ """


  Parameters
  ----------
    const Waveform& B
  
  Returns
  -------
    Waveform
  
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

%feature("docstring") GWFrames::Quaternion::normalized """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Waveform::HistoryStream """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    stringstream&
  
"""

%feature("docstring") GWFrames::FactorialFunctor::operator() """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Union """
Return the union of two time sequences.
=======================================
  Parameters
  ----------
    const vector<double>& t1
    const vector<double>& t2
    const double MinStep = 0.005
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    On the overlap between the two sequences, the time is built up by taking
    the smaller time step in either of the two sequences, or MinStep if that
    step is smaller.
    
    The input to this function is assumed to be strictly monotonic.
  
"""

%feature("docstring") GWFrames::MatrixC::operator[] """


  Parameters
  ----------
    const int i
  
  Returns
  -------
    complex<double> *
  
"""

%feature("docstring") GWFrames::MatrixC::operator[] """


  Parameters
  ----------
    const int i
  
  Returns
  -------
    const complex<double> *
  
"""

%feature("docstring") GWFrames::vec """


  Parameters
  ----------
    const Quaternion& Q
  
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

%feature("docstring") GWFrames::Waveform::GetAlignmentOfTime """
Find the time offset aligning this waveform to the other at the fiducial time.
==============================================================================
  Parameters
  ----------
    const Waveform& A
      Fixed Waveform in inertial frame to which this Waveform is aligned
    const double t_fid
    double& deltat
      The value to be returned
  
  Returns
  -------
    void
  
  Description
  -----------
    This function simply finds the appropriate time offset, rather than
    applying it. This is called by AlignTime and probably does not need to be
    called directly; see that function's documentation for more details.
    
    AlignTime
  
"""

%feature("docstring") GWFrames::conjugate """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  
"""

%feature("docstring") GWFrames::Squad """
Squad interpolation of Quaternion time series.
==============================================
  Parameters
  ----------
    const vector<Quaternion>& RIn
      Vector of rotors
    const vector<double>& tIn
      Vector of corresponding times
    const vector<double>& tOut
      Vector of times to which RIn will be interpolated
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    This function implements a version of cubic-spline interpolation designed
    for unit quaternions, which delivers more accurate, smooth, and physical
    rotations than other forms of interpolation.
  
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

%feature("docstring") GWFrames::Quaternion::inverse """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Component3 """


  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::normalized """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  
"""

%feature("docstring") GWFrames::Component1 """


  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::Interpolate """
Interpolate the Waveform to a new set of time instants.
=======================================================
  Parameters
  ----------
    const vector<double>& NewTime
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") GWFrames::Waveform::SetMIsScaledOut """


  Parameters
  ----------
    const bool Scaled
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::HistoryStr """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    string
  
"""

%feature("docstring") GWFrames::Waveform::AlignFrame """
Change this Waveform by aligning the frame to the other's at the given time.
============================================================================
  Parameters
  ----------
    const Waveform& A
      Fixed Waveform in corotating frame to which this Waveform is aligned
    const double t_fid
      Fiducial time at which to equate frames
  
  Returns
  -------
    Waveform&
  
  Description
  -----------
    Note that this function operates in place; the Waveform to which it is
    applied will change. However, the modes are not altered; only the frame
    data is.
    
    As noted above, it is implicitly assumed that both Waveforms are in their
    corotating frames, with the modes appropriately aligned to the frames at
    t_fid. The assumption is that the frames actually represent something
    physically meaningful, so that it is meaningful to insist that they be the
    same.
    
    Then, this function aligns the frames at t_fid by multiplying this->frame
    on the left by a constant rotor such that this->frame at t_fid is exactly
    A.frame at t_fid. The resulting frame is now corotating with an
    angular-velocity vector that has been rotated by that constant rotor,
    relative to the inertial basis.
    
    AlignDecompositionFrameToModes
  
"""

%feature("docstring") GWFrames::Waveform::NPEdthBar """
Newman-Penrose edth operator conjugate.
=======================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This operator is the one defined by Newman and Penrose (1966) and further
    described by Goldberg et al. (1967). It lowers the spin weight of any field
    on the sphere by 1. Note that this operator does not preserve boost weights
    in any nice way  except in special cases. The GHP version does.
    
    Note that the boost weight is set to the value of WeightError, which is
    just meant to be large enough that it will give improbable values if used.
    This is not fool-proof.
    
    NPEdth
    
    GHPEdth
    
    GHPEdthBar
  
"""

%feature("docstring") GWFrames::Quaternion::vec """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::CumulativeVectorIntegral """
Integrate vector function by simple trapezoidal rule.
=====================================================
  Parameters
  ----------
    const vector<vector<double>>& fdot
      Vector of vectors (first index time).
    const vector<double>& t
      Vector of corresponding time steps.
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::dot """


  Parameters
  ----------
    const Quaternion& Q
    const Quaternion& P
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::WignerCoefficientFunctor::operator() """


  Parameters
  ----------
    const int ell
    const int mp
    const int m
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Quaternion::normsquared """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::SWSH """
class GWFrames::SWSH
====================
  Object for computing values of the spin-weighted spherical harmonics.
  
  Member variables
  ----------------
    WignerDMatrix D
    int spin
    double sign
     __pad0__
  
"""

%feature("docstring") GWFrames::Quaternion::cross """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::LadderOperatorFactorFunctor::LadderOperatorFactorFunctor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    LadderOperatorFactorFunctor
  
"""

%feature("docstring") GWFrames::Waveform::ResizeData """


  Parameters
  ----------
    const unsigned int NModes
    const unsigned int NTimes
  
  Returns
  -------
    void
  
"""

%feature("docstring") WignerDMatrix::operator() """
Evaluate the D matrix element for the given (ell, mp, m) indices.
=================================================================
  Parameters
  ----------
    const int ell
    const int mp
    const int m
  
  Returns
  -------
    complex<double>
  
"""

%feature("docstring") GWFrames::Waveform::PNEquivalentPrecessionalAV """
Deduce PN-equivalent precessional angular velocity from Waveform.
=================================================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate
  
  Returns
  -------
    vector<vector<double>>
  
  Description
  -----------
    This function subtracts the PN-equivalent orbital angular velocity (given
    by PNEquivalentOrbitalAV) from the field's angular velocity. This should be
    equivalent to the precessional angular velocity of the PN system. Note that
    the returned vector is relative to the inertial frame.
    
    This may be essentially numerical noise if there is no precession, or if
    precession has oscillated to zero.
    
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.
    
    PNEquivalentOrbitalAV
  
"""

%feature("docstring") SchmidtEtAl_fdf """


  Parameters
  ----------
    const gsl_vector * v
    void * params
    double * f
    gsl_vector * df
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::ScalarDerivative """
Three-point finite-differencing of vector of scalars.
=====================================================
  Parameters
  ----------
    const vector<double>& f
      Vector of scalars.
    const vector<double>& t
      Vector of corresponding time steps.
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    Sundquist and Veronis, Tellus XXII (1970), 1
  
"""

%feature("docstring") GWFrames::Waveform::GHPEdth """
Geroch-Held-Penrose edth operator.
==================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This operator is the one defined by Geroch et al. (1973). It raises the
    spin weight of any field on the sphere by 1, while leaving the boost weight
    unchanged.
    
    This operator is very similar to the basic Newman-Penrose edth operator,
    except that it preserves boost weights. In certain cases, it reduces to the
    NP edth (up to a factor of $\\sqrt{2}$), but not generally. In those
    certain cases, we have NPEdth() = sqrt(2)*GHPEdth().
    
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
    
    NPEdth
    
    NPEdthBar
    
    GHPEdthBar
  
"""

%feature("docstring") GWFrames::DifferentiateRotorByLogarithm """
Calculate the derivative of a rotor by the logarithm formula.
=============================================================
  Parameters
  ----------
    const vector<Quaternion>& RIn
    const vector<double>& tIn
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    This is a much more complicated way of evaluating the derivative of a
    quaternion function of time, as compared to finite differencing by
    'QuaternionDerivative'. However, there may be accuracy advantages when the
    logarithm is smooth, and  at the least  this can serve as a good test of
    the correctness of the logarithm formula.
  
"""

%feature("docstring") GWFrames::Matrix::set """


  Parameters
  ----------
    const unsigned int r
    const unsigned int c
    const double v
  
  Returns
  -------
    Matrix&
  
"""

%feature("docstring") MatrixC::swap """


  Parameters
  ----------
    MatrixC& b
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::cross """


  Parameters
  ----------
    const Quaternion& Q
    const Quaternion& P
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") LadderOperatorFactorCalculator """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    We need (2*ell+1) coefficients for each value of ell from 0 (for
    completeness) up to ellMax_Utilities (hard coded in the header file).
    That's a total of from sympy import summation, symbols ell,
    ellMax_Utilities, m, mp = symbols('ell ellMax_Utilities m mp',
    integer=True) summation(2*ell+1, (ell, 0, ellMax_Utilities))
    
    ellMax_Utilities**2 + 2*ellMax_Utilities + 1 With a similar calculation, we
    can see that the associated access operator needs element
    summation(2*ell+1, (ell, 0, ell-1)) + ell + m
    
    ell**2 + ell + m
  
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

%feature("docstring") GWFrames::DominantPrincipalAxis """


  Parameters
  ----------
    Matrix& M
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::FrameFromZ """
Construct minimal-rotation frame from Z basis vector of that frame.
===================================================================
  Parameters
  ----------
    const vector<Quaternion>& Z
      Vector of Quaternions
    const vector<double>& T
      Vector of corresponding times
    const unsigned int NIterations = 5
      Number of refinements [default: 5]
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    The input vector of Quaternions, assumed to be pure unit vectors, represent
    the Z basis vectors of the frame at each instant of time. The returned
    vector of rotors will rotate the stationary frame's (x,y,z) vectors into
    the new frame's (X,Y,Z) vectors. The X and Y vectors are deduced by
    imposing the minimal-rotation condition. Note that this leaves an unfixed
    initial rotation about z.
  
"""

%feature("docstring") GWFrames::Determinant """


  Parameters
  ----------
    Matrix& M
  
  Returns
  -------
    double
  
"""

%feature("docstring") Matrix::resize """


  Parameters
  ----------
    unsigned int newNRows
    unsigned int newNCols
    const double  = 0.0
  
  Returns
  -------
    void
  
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

%feature("docstring") GWFrames::Waveform::DataType """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::SWSH::operator() """


  Parameters
  ----------
    const int ell
    const int m
  
  Returns
  -------
    complex<double>
  
"""

%feature("docstring") GWFrames::Waveform::RIsScaledOut """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    bool
  
"""

%feature("docstring") GWFrames::SWSH::sign """


  Parameters
  ----------
    s% 2 = =0?1.0:-1.0
  
  Returns
  -------
    sign
  
"""

%feature("docstring") GWFrames::Quaternion::angle """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Waveform::RotateDecompositionBasis """
Rotate the basis in which this Waveform is measured by a constant rotor.
========================================================================
  Parameters
  ----------
    const Quaternion& R_frame
  
  Returns
  -------
    Waveform&
  

Rotate the basis in which this Waveform is measured.
====================================================
  Parameters
  ----------
    const vector<Quaternion>& R_frame
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

%feature("docstring") GWFrames::Slerp """


  Parameters
  ----------
    const double tau
    const Quaternion& Qa
    const Quaternion& Qb
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Waveform::SetBoostWeight """


  Parameters
  ----------
    const int NewBoostWeight
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::RotateDecompositionBasisOfUncertainties """
Rotate the basis in which this Waveform's uncertainties are measured.
=====================================================================
  Parameters
  ----------
    const vector<Quaternion>& R_frame
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

%feature("docstring") GWFrames::Waveform::SetFrame """


  Parameters
  ----------
    const vector<Quaternion>& a
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Quaternion::operator== """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    bool
  
"""

%feature("docstring") GWFrames::Waveform::FindModeIndexWithoutError """


  Parameters
  ----------
    const int L
    const int M
  
  Returns
  -------
    unsigned int
  
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

%feature("docstring") GWFrames::LadderOperatorFactorFunctor::operator() """


  Parameters
  ----------
    const int ell
    const int m
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::operator/ """


  Parameters
  ----------
    const double a
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const double a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const Quaternion& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const Quaternion& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const double a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& Q
    const vector<double>& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const vector<double>& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& P
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const double b
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const vector<vector<double>>& a
    const vector<double>& b
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::PrescribedRotation """
Input frame with prescribed rate of rotation about Z axis.
==========================================================
  Parameters
  ----------
    const vector<double>& RotationRateAboutZ
      Vector of rotation rates about the new frame's Z axis.
    const vector<Quaternion>& R
      Vector of rotors.
    const vector<double>& T
      Vector of corresponding time steps.
    const unsigned int NIterations = 5
      Number of refinements [default: 5]
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    This function returns a copy of the input R, which takes the z axis to the
    same point as R, but adjusts the rotation about that new point by imposing
    the minimal-rotation condition, and then including an additional rotation
    about the new Z axis to agree with the given rotation rate.
  
"""

%feature("docstring") GWFrames::operator- """


  Parameters
  ----------
    const double a
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const double a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const Quaternion& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& a
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const Quaternion& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const double a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const Quaternion& Q
    const vector<double>& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& Q
    const vector<double>& a
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<Quaternion>& P
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  



  Parameters
  ----------
    const vector<double>& a
    const vector<double>& b
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const vector<double>& a
    const double b
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::commutator """


  Parameters
  ----------
    const Quaternion& Q
    const Quaternion& P
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::sqrt """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  
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

%feature("docstring") GWFrames::Matrix::ncols """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") GWFrames::Waveform::FrameType """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::RDelta """
Difference between frame rotors.
================================
  Parameters
  ----------
    const vector<Quaternion>& R1
      Vector of rotors
    const vector<Quaternion>& R2
      Vector of rotors
    const unsigned int IndexOfFiducialTime = 0
      Integer index of time at which difference is set to zero [default: 0]
  
  Returns
  -------
    vector<Quaternion>
  
"""

%feature("docstring") GWFrames::MatrixC """
class GWFrames::MatrixC
=======================
  Rectangular array of complex data; probably not needed directly.
  
  Member variables
  ----------------
    int nn
    int mm
    complex<double> ** v
  
"""

%feature("docstring") Matrix::clear """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::FactorialFunctor """
class GWFrames::FactorialFunctor
================================
  Object for pre-computing and retrieving factorials.
  
  Member variables
  ----------------
    const vector<double> FactorialTable
  
"""

%feature("docstring") GWFrames::Unwrap """
Unwrap phase so that it is (roughly) continuous.
================================================
  Parameters
  ----------
    const vector<double>& In
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Component2 """


  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::sqrtOfRotor """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  
"""

%feature("docstring") NumberToString """


  Parameters
  ----------
    typename T 
    T Number
  
  Returns
  -------
    typename T
  
"""

%feature("docstring") GWFrames::operator<< """
Print the quaternion nicely to stream.
======================================
  Parameters
  ----------
    ostream& out
    const Quaternion& q
  
  Returns
  -------
    ostream&
  
"""

%feature("docstring") GWFrames::Quaternion::operator[] """
Get component of Quaternion.
============================
  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    double
  
  Description
  -----------
    The 0 component is the scalar part, and the 13 components are the vector
    components.
  

Get reference to component of Quaternion.
=========================================
  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    double&
  
  Description
  -----------
    Note: This is unavailable from python.
  
"""

%feature("docstring") WaveformModes """


  Parameters
  ----------
    const double delta
    const double v
    const double chisl
    const double chial
    vector<complex<double>>& modes
  
  Returns
  -------
    void
  
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

%feature("docstring") GWFrames::UnflipRotors """
Remove sign-ambiguity of rotors.
================================
  Parameters
  ----------
    const vector<Quaternion>& R
      Vector of rotors
    const double discont = 1.4142135623730951
      Acceptable discontinuity [default: sqrt(2)]
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    Because of the two-sided nature of quaternion rotations, the sign of a
    rotor may be undetermined in many cases. Discontinuous flips in that sign
    for rotor-valued functions of time can cause significant problems. This
    function removes those flips by ensuring that the output rotors at
    successive instants are within 'discont' of each other.
  
"""

%feature("docstring") GWFrames::Waveforms """
class GWFrames::Waveforms
=========================
  Object storing a collection of Waveform objects to be operated on uniformly.
  
  Member variables
  ----------------
    vector<Waveform> Ws
    bool CommonTimeSet
  
"""

%feature("docstring") GWFrames::Waveform::DescriptorString """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    string
  
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

%feature("docstring") GWFrames::Matrix """
class GWFrames::Matrix
======================
  3x3 object wrapping GSL matrix; probably not needed directly
  
  Member variables
  ----------------
    gsl_matrix * m
  
"""

%feature("docstring") BinomialCoefficientCalculator """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    We need (n+1) coefficients for each value of n from 0 (for completeness) up
    to 2*ellMax_Utilities (hard coded in the header file). That's a total of
    from sympy import summation, symbols ellMax_Utilities, n, k =
    symbols('ellMax_Utilities n k', integer=True) summation(n+1, (n, 0,
    2*ellMax_Utilities))
    
    2*ellMax_Utilities**2 + 3*ellMax_Utilities + 1 With a similar calculation,
    we can see that the associated access operator needs element (n*(n+1)/2 +
    k) of the array.
  
"""

%feature("docstring") GWFrames::Waveforms::size """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") GWFrames::abs """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    double
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const vector<double>& v
  
  Returns
  -------
    double
  



  Parameters
  ----------
    const vector<vector<double>>& v
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Quaternion::exp """
Return exponent of Quaternion.
==============================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Waveform::OShaughnessyEtAlVector """
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

%feature("docstring") tolower """


  Parameters
  ----------
    const string& A
  
  Returns
  -------
    string
  
"""

%feature("docstring") GWFrames::Matrix::nrows """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
"""

%feature("docstring") GWFrames::Quaternion::dot """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames """
namespace GWFrames
==================
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

%feature("docstring") GWFrames::Quaternion::sqrt """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  
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

%feature("docstring") GWFrames::Matrix::operator() """


  Parameters
  ----------
    const unsigned int row
    const unsigned int col
  
  Returns
  -------
    double
  



  Parameters
  ----------
    const unsigned int row
    const unsigned int col
  
  Returns
  -------
    double&
  
"""

%feature("docstring") GWFrames::Component0 """


  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::SetSpinWeight """


  Parameters
  ----------
    const int NewSpinWeight
  
  Returns
  -------
    void
  
"""

%feature("docstring") GWFrames::Waveform::TransformUncertaintiesToCorotatingFrame """
Transform Waveform uncertainties to corotating frame.
=====================================================
  Parameters
  ----------
    const vector<Quaternion>& R_frame
      Vector of rotors giving corotating frame of the data.
  
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
    void
  
"""

%feature("docstring") GWFrames::VectorIntegral """
Integrate vector function by simple trapezoidal rule.
=====================================================
  Parameters
  ----------
    const vector<vector<double>>& fdot
      Vector of vectors (first index time).
    const vector<double>& t
      Vector of corresponding time steps.
  
  Returns
  -------
    vector<vector<double>>
  
"""

%feature("docstring") GWFrames::log """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  



  Parameters
  ----------
    const vector<Quaternion>& Q
  
  Returns
  -------
    vector<Quaternion>
  
"""

%feature("docstring") GWFrames::Quaternion::conjugate """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Waveform::MaxNormTime """


  Parameters
  ----------
    const unsigned int SkipFraction = 4
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Waveform::GetAlignmentOfDecompositionFrameToModes """
Find the appropriate rotation to fix the orientation of the corotating frame.
=============================================================================
  Parameters
  ----------
    const double t_fid
      Fiducial time at which the alignment should happen
    Quaternion& R_delta
      Returned rotor
    const vector<int>& Lmodes = vector<int>(0)
      Lmodes to use in computing $<LL>$
  
  Returns
  -------
    void
  
  Description
  -----------
    This function simply finds the rotation necessary to align the corotating
    frame to the waveform at the fiducial time, rather than applying it. This
    is called by AlignDecompositionFrameToModes and probably does not need to
    be called directly; see that function's documentation for more details.
    
    AlignDecompositionFrameToModes
  
"""

%feature("docstring") Matrix::operator* """


  Parameters
  ----------
    const vector<double>& b
  
  Returns
  -------
    vector<double>
  



  Parameters
  ----------
    const Quaternion& b
  
  Returns
  -------
    Quaternion
  
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
  
"""

%feature("docstring") Matrix::operator- """


  Parameters
  ----------
    const Matrix& rhs
  
  Returns
  -------
    Matrix
  
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

%feature("docstring") GWFrames::PNWaveform::chi1Mag """


  Parameters
  ----------
    const unsigned int iTime
  
  Returns
  -------
    double
  
"""

%feature("docstring") GWFrames::Quaternion::log """
Return logarithm of Quaternion.
===============================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::MatrixC::ncols """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    int
  
"""

%feature("docstring") GWFrames::MinimalRotation """
Minimal-rotation version of the input frame.
============================================
  Parameters
  ----------
    const vector<Quaternion>& R
      Vector of rotors.
    const vector<double>& T
      Vector of corresponding time steps.
    const unsigned int NIterations = 5
      Number of refinements [default: 5]
  
  Returns
  -------
    vector<Quaternion>
  
  Description
  -----------
    This function returns a copy of the input R, which takes the z axis to the
    same point as R, but adjusts the rotation about that new point by imposing
    the minimal-rotation condition.
  
"""

%feature("docstring") GWFrames::Quaternion::commutator """


  Parameters
  ----------
    const Quaternion& Q
  
  Returns
  -------
    Quaternion
  
"""

%feature("docstring") GWFrames::Waveform::Segment """
Extract a segment of a Waveform.
================================
  Parameters
  ----------
    const unsigned int i1
      Index of initial time
    const unsigned int i2
      Index just beyond final time
  
  Returns
  -------
    Waveform
  
"""

%feature("docstring") GWFrames::Waveforms::Waveforms """
Empty constructor of N empty objects.
=====================================
  Parameters
  ----------
    const int N = 0
  
  Returns
  -------
    Waveforms
  
  Description
  -----------
    Waveforms (plural!) //
  

Basic copy constructor.
=======================
  Parameters
  ----------
    const Waveforms& In
  
  Returns
  -------
    Waveforms
  

Basic copy constructor.
=======================
  Parameters
  ----------
    const vector<Waveform>& In
  
  Returns
  -------
    Waveforms
  
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

%feature("docstring") GWFrames::ScalarIntegral """
Integrate scalar function by simple trapezoidal rule.
=====================================================
  Parameters
  ----------
    const vector<double>& fdot
      Vector of scalars.
    const vector<double>& t
      Vector of corresponding time steps.
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") Matrix::operator= """


  Parameters
  ----------
    const Matrix& rhs
  
  Returns
  -------
    Matrix&
  



  Parameters
  ----------
    const vector<vector<double>>& newData
  
  Returns
  -------
    Matrix&
  
"""

