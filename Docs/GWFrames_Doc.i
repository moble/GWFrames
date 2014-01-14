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
  



  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
"""

%feature("docstring") GWFrames::Waveform::SetRIsScaledOut """


  Parameters
  ----------
    const bool Scaled
  
  Returns
  -------
    void
  
"""

%feature("docstring") SQR """


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
    
    NPEdth
    
    NPEdthBar
    
    GHPEdth
    
    IntegrateNPEdth
    
    IntegrateNPEdthBar
    
    IntegrateGHPEdth
    
    IntegrateGHPEdthBar
  
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
    in any nice way  except in special cases. The GHP version does. Note that,
    in this implementation, the only difference between the NP and GHP versions
    is the factor of $\\sqrt{2}$. The additional GHP term that keeps the boost
    weight meaningful is zero in any given frame  though it transforms
    nontrivially.
    
    Note that the boost weight is set to the value of WeightError, which is
    just meant to be large enough that it will give improbable values if used.
    This is not fool-proof.
    
    NPEdthBar
    
    GHPEdth
    
    GHPEdthBar
    
    IntegrateNPEdth
    
    IntegrateNPEdthBar
    
    IntegrateGHPEdth
    
    IntegrateGHPEdthBar
  
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

%feature("docstring") GWFrames::Waveform::IntegrateGHPEdthBar """
Integrate the Geroch-Held-Penrose edth operator conjugate.
==========================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This operator inverts the action of the GHP edth operator. This is not a
    perfect inverse, because the l=s-1 term is set to zero. To be precise, if
    Waveform A has spins weight $s$, then A.GHPEdth().IntegrateGHPEdth() has
    the effect of setting the $\\ell=s$ term in A to zero.
    
    NPEdth
    
    NPEdthBar
    
    GHPEdth
    
    GHPEdthBar
    
    IntegrateNPEdth
    
    IntegrateNPEdthBar
    
    IntegrateGHPEdth
  
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
    void
  
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
    void
  
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

%feature("docstring") InverseConformalFactorFunctor::operator() """


  Parameters
  ----------
    const Quaternions::Quaternion& R
  
  Returns
  -------
    double
  
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

%feature("docstring") GWFrames::Waveforms::clear """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    void
  
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

%feature("docstring") StringForm """


  Parameters
  ----------
    const vector<int>& Lmodes
  
  Returns
  -------
    string
  
"""

%feature("docstring") Boost """
Return a rotor taking n into its boosted version.
=================================================
  Parameters
  ----------
    ThreeVector v
      Three-vector velocity of the new frame WRT this frame
    ThreeVector n
      Three-vector direction to be boosted by the rotor
  
  Returns
  -------
    Quaternion
  
  Description
  -----------
    This function returns a rotor $R_b$ that takes the input vector $\\hat{n}$
    (which will be normalized) on the future null sphere into its boosted
    version. Note that this rotor is a function of both the vector being
    boosted and the boost itself.
  
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
    Quaternions::Quaternion& R_delta
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

%feature("docstring") GWFrames::PNWaveform::PNWaveform """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    PNWaveform
  



  Parameters
  ----------
    const PNWaveform& W
  
  Returns
  -------
    PNWaveform
  



  Parameters
  ----------
    const string& Approximant
    const double delta
    const vector<double>& chi1_i
    const vector<double>& chi2_i
    const double Omega_orb_i
    const Quaternions::Quaternion& R_frame_i = Quaternions::Quaternion(1, 0, 0, 0)
    const double PNOrder = 4.0
    double v_0 = -1.0
  
  Returns
  -------
    PNWaveform
  
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
    Quaternions::Quaternion& R_delta
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

%feature("docstring") GWFrames::Waveform::Boost """
Apply a boost to a boost-weighted function.
===========================================
  Parameters
  ----------
    const vector<double>& v
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This function does three things. First, it evaluates the Waveform on what
    will become an equi-angular grid after transformation by the boost. Second,
    it multiplies each of those points by the appropriate conformal factor
    $K^b(\\vartheta, \\varphi)$, where $b$ is the boost weight stored with the
    Waveform. Finally, it transforms back to Fourier space using that new
    equi-angular grid.
  
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

%feature("docstring") GWFrames::Waveform::IntegrateGHPEdth """
Integrate the Geroch-Held-Penrose edth operator.
================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This operator inverts the action of the GHP edth operator. This is not a
    perfect inverse, because the l=s-1 term is set to zero. To be precise, if
    Waveform A has spins weight $s$, then A.GHPEdth().IntegrateGHPEdth() has
    the effect of setting the $\\ell=s$ term in A to zero.
    
    NPEdth
    
    NPEdthBar
    
    GHPEdth
    
    GHPEdthBar
    
    IntegrateNPEdth
    
    IntegrateNPEdthBar
    
    IntegrateGHPEdthBar
  
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

%feature("docstring") GWFrames::Waveform::SetHistory """


  Parameters
  ----------
    const string& Hist
  
  Returns
  -------
    void
  
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

%feature("docstring") SchmidtEtAl_f """


  Parameters
  ----------
    const gsl_vector * v
    void * params
  
  Returns
  -------
    double
  
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

%feature("docstring") GWFrames::Waveform::IntegrateNPEdthBar """
Integrate the Newman-Penrose edth operator conjugate.
=====================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This operator inverts the action of the conjugated Newman-Penrose edth
    operator. This is not a perfect inverse, because the l=s-1 term is set to
    zero. To be precise, if Waveform A has spin weight $s$, then
    A.NPEdthBar().IntegrateNPEdthBar() has the effect of setting the $\\ell=s$
    term in A to zero.
    
    Note that the N-P edth operator does not preserve boost weights, so the
    boost weight is set to the value of WeightError, which is just meant to be
    large enough that it will give improbable values if used. This is not
    fool-proof. See the GHP edth operator for a weight-preserving version.
    
    NPEdth
    
    NPEdthBar
    
    GHPEdth
    
    GHPEdthBar
    
    IntegrateNPEdthBar
    
    IntegrateGHPEdth
    
    IntegrateGHPEdthBar
  
"""

%feature("docstring") GWFrames::PNWaveform::Omega_tot """


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
    void
  
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

%feature("docstring") DataGrid::operator/ """


  Parameters
  ----------
    const DataGrid& 
  
  Returns
  -------
    DataGrid
  
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
  
  Returns
  -------
    Waveform
  
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
    in any nice way  except in special cases. The GHP version does. Note that,
    in this implementation, the only difference between the NP and GHP versions
    is the factor of $\\sqrt{2}$. The additional GHP term that keeps the boost
    weight meaningful is zero in any given frame  though it transforms
    nontrivially.
    
    Note that the boost weight is set to the value of WeightError, which is
    just meant to be large enough that it will give improbable values if used.
    This is not fool-proof.
    
    NPEdth
    
    GHPEdth
    
    GHPEdthBar
    
    IntegrateNPEdth
    
    IntegrateNPEdthBar
    
    IntegrateGHPEdth
    
    IntegrateGHPEdthBar
  
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

%feature("docstring") GWFrames::ScriFunctor """
class GWFrames::ScriFunctor
===========================
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
    
    NPEdth
    
    NPEdthBar
    
    GHPEdthBar
    
    IntegrateNPEdth
    
    IntegrateNPEdthBar
    
    IntegrateGHPEdth
    
    IntegrateGHPEdthBar
  
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

%feature("docstring") std """
namespace std
=============
  STL namespace.
  
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
  



  Parameters
  ----------
    const double b
    const Waveform& A
  
  Returns
  -------
    Waveform
  
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
  
"""

%feature("docstring") GWFrames::Waveform::IntegrateNPEdth """
Integrate the Newman-Penrose edth operator.
===========================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This operator inverts the action of the Newman-Penrose edth operator. This
    is not a perfect inverse, because the l=s-1 term is set to zero. To be
    precise, if Waveform A has spin weight $s$, then
    A.NPEdth().IntegrateNPEdth() has the effect of setting the $\\ell=s$ term
    in A to zero.
    
    Note that the N-P edth operator does not preserve boost weights, so the
    boost weight is set to the value of WeightError, which is just meant to be
    large enough that it will give improbable values if used. This is not
    fool-proof. See the GHP edth operator for a weight-preserving version.
    
    NPEdth
    
    NPEdthBar
    
    GHPEdth
    
    GHPEdthBar
    
    IntegrateNPEdthBar
    
    IntegrateGHPEdth
    
    IntegrateGHPEdthBar
  
"""

%feature("docstring") GWFrames::PNWaveform::Omega_totMag """


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

%feature("docstring") GWFrames::DataGrid::SetNTheta """


  Parameters
  ----------
    const int N_theta
  
  Returns
  -------
    DataGrid&
  
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

%feature("docstring") GWFrames::Waveforms::~Waveforms """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    ~Waveforms
  
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

%feature("docstring") GWFrames::Waveforms """
class GWFrames::Waveforms
=========================
  Object storing a collection of Waveform objects to be operated on uniformly.
  
  Member variables
  ----------------
    vector<Waveform> Ws
    bool CommonTimeSet
  
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

%feature("docstring") GWFrames::Waveforms::size """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    unsigned int
  
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

%feature("docstring") GWFrames::Waveform::ApplySupertranslation """
Re-interpolate data to new time slices given by this supertranslation.
======================================================================
  Parameters
  ----------
    vector<complex<double>>& gamma
  
  Returns
  -------
    Waveform
  
  Description
  -----------
    This function takes the current data decomposed as spherical harmonics on a
    given slicing, transforms to physical space, re-interpolates the data at
    each point to a new set of time slices, and transforms back to
    spherical-harmonic coefficients.
    
    The supertranslation data input gamma is a vector of complex numbers
    representing the (scalar) spherical-harmonic components of the
    supertranslation, stored in the order (0,0), (1,-1), (1,0), (1,1), (2,-2),
    ... The overall time translation is given by the first component; the
    spatial translation is given by the second through fourth componentes;
    higher components give the proper supertranslations. In particular, a
    proper supertranslation will have its first four coefficients equal to 0.0.
    
    Note that, for general spin-weighted spherical-harmonic components
    ${}_{s}a_{l,m}$, a real function results when ${}_{-s}a_{l,-m} =
    {}_{s}a_{l,m}^\\ast$. In particular, the input gamma data are assumed to
    satisfy this formula with $s=0$.
  
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

%feature("docstring") GWFrames::Waveform::operator() """


  Parameters
  ----------
    const unsigned int Mode
    const unsigned int TimeIndex
  
  Returns
  -------
    complex<double>
  
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
    const vector<Quaternions::Quaternion>& R_frame
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

%feature("docstring") GWFrames::Modes::Data """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<complex<double>>
  
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
    Quaternions::Quaternion& R_delta
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

%feature("docstring") GWFrames::Waveform::SetFrame """


  Parameters
  ----------
    const vector<Quaternions::Quaternion>& a
  
  Returns
  -------
    void
  
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
    void
  
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
Waveforms (plural!) //.
=======================
  Parameters
  ----------
    const int N = 0
  
  Returns
  -------
    Waveforms
  
  Description
  -----------
    Empty constructor of N empty objects.
  

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

%feature("docstring") Rapidity """
Returns the rapidity of a Lorentz boost with velocity three-vector v.
=====================================================================
  Parameters
  ----------
    const vector<double>& v
  
  Returns
  -------
    double
  
  Description
  -----------
    The vector v is expected to be the velocity three-vector of the new frame
    relative to the current frame, in units where c=1.
  
"""

