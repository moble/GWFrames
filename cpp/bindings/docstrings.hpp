// Docstrings for the Python bindings, converted once from the
// Doxygen-generated SWIG files `Docs/GWFrames_Doc.i` and
// `Quaternions/docs/Quaternions_Doc.i`.  Keys omit the leading
// `GWFrames::` or `Quaternions::` namespace.

#ifndef GWFRAMES_BINDINGS_DOCSTRINGS_HPP
#define GWFRAMES_BINDINGS_DOCSTRINGS_HPP

#include <string>
#include <unordered_map>

namespace GWFramesBindings {

inline const std::unordered_map<std::string, const char*>& Docstrings() {
  static const std::unordered_map<std::string, const char*> docstrings = {
    {"AdvLIGO_NSNSOptimal", R"doc(Parameters
----------
  const vector<double>& F
  const bool Invert = false
  const double NoiseFloor = 0.0

Returns
-------
  vector<double>)doc"},
    {"AdvLIGO_ZeroDet_HighP", R"doc(Parameters
----------
  const vector<double>& F
  const bool Invert = false
  const double NoiseFloor = 0.0

Returns
-------
  vector<double>)doc"},
    {"AdvLIGO_ZeroDet_LowP", R"doc(Parameters
----------
  const vector<double>& F
  const bool Invert = false
  const double NoiseFloor = 0.0

Returns
-------
  vector<double>)doc"},
    {"AlignWaveforms", R"doc(Do everything necessary to align two waveform objects.
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
    const bool Debug = false

  Returns
  -------
    void

  Description
  -----------
    This function aligns the frame to the waveform modes for both input
    Waveform objects at time t_mid = (t_1+t_2)/2. It also optimizes the
    alignment of W_B by adjusting its time and overall attitudes to align with
    W_A as well as possible. While doing so, it re-adjusts the frame alignment
    to the modes for W_B to account for the changing meaning of t_mid.

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
    remotely well sampled, and assuming the code's logic does not fail due to
    otherwise weird or unexpected inputs, I (Mike Boyle) do hereby guarantee
    that this algorithm will find the optimal alignment in both time and
    attitude. Or your money back.)doc"},
    {"ApproximateMeanRotor", R"doc('Mean' rotor, minimizing 'distance' between R(t) and R_mean
===========================================================
  Parameters
  ----------
    const vector<Quaternion>& R
      Vector of rotors to be averaged.

  Returns
  -------
    Quaternion

  Description
  -----------
    The true mean of a rotor function of time can only be calculated by an
    optimization algorithm because the metric is too complicated to solve
    analytically. Instead, this algorithm uses a simple approximation that is
    true in the limit of R(t) being constant. In particular, the ideal mean
    would be calculated using the angle between R(t) and R_mean; this function
    minimizes the squared sin of half that angle. Another way of thinking of
    this is that this function isn't approximate; it just uses a different
    metric, which is approximately equal to what we would like to use.

    This algorithm calculates the averaging formula found in 'Means
and
    averaging in the group of rotations'
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.16.5040.


'Mean' rotor, approximately minimizing integrated 'distance' between R(t) and R_mean
====================================================================================
  Parameters
  ----------
    const vector<Quaternion>& R
      Vector of rotors to be averaged.
    const vector<double>& t
      Vector of corresponding time steps.
    const double t1 = -1e300
      Optional early bound of time to average over
    const double t2 = 1e300
      Optional late bound of time to average over

  Returns
  -------
    Quaternion

  Description
  -----------
    The true mean of a rotor function of time can only be calculated by an
    optimization algorithm because the metric is too complicated to solve
    analytically. Instead, this algorithm uses a simple approximation that is
    true in the limit of R(t) being constant. In particular, the ideal mean
    would be calculated using the angle between R(t) and R_mean; this function
    minimizes the squared sin of half that angle. Another way of thinking of
    this is that this function isn't approximate; it just uses a different
    metric, which is approximately equal to what we would like to use.

    This algorithm calculates the continuous limit of the averaging formula
    found in 'Means and averaging in the group of rotations'
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.16.5040.)doc"},
    {"ApproximateOptimalAlignment", R"doc(Time and rotor offsets such that Rdelta*Rb(t+deltat) is as close to Ra(t) as possible.
======================================================================================
  Parameters
  ----------
    const double t1
      Early bound of time to average over
    const double t2
      Late bound of time to average over
    const vector<Quaternion>& Ra
      Vector of fixed rotors.
    const vector<double>& ta
      Vector of corresponding time steps.
    const vector<Quaternion>& Rb
      Vector of rotors to be aligned to Ra.
    const vector<double>& tb
      Vector of corresponding time steps.
    double& deltat
      Output optimal time offset
    Quaternions::Quaternion& R_delta
      Output optimal rotor

  Returns
  -------
    void

  Description
  -----------
    This routine actively optimizes over the time offset of Rb and
    approximately optimizes over the rotational degrees of freedom in Rb such
    that Rdelta*Rb(t+deltat) is as close to Ra(t) as possible. The rotational
    optimization for any given time offset is performed analytically by the
    routine ApproximateOptimalAlignmentRotor.

    The times t1 and t2 are measured relative to ta, and all are left fixed; Rb
    and tb are shifted to achieve alignment.)doc"},
    {"ApproximateOptimalAlignmentRotor", R"doc(Approximate Rdelta such that Rdelta*Rb is as nearly equal to Ra as possible.
============================================================================
  Parameters
  ----------
    const vector<Quaternion>& Ra
      Vector of fixed rotors.
    const vector<Quaternion>& Rb
      Vector of rotors to be aligned to Ra.
    const vector<double>& t
      Vector of corresponding time steps.
    const double t1 = -1e300
      Optional early bound of time to average over
    const double t2 = 1e300
      Optional late bound of time to average over

  Returns
  -------
    Quaternion

  Description
  -----------
    This routine essentially just returns the MeanRotor of the difference
    between Ra and Rb. The result is an Rdelta which can be applied to Rb, by
    multiplication on the left, to give a new time-series for Rb that is
    approximately as close to Ra as possible.)doc"},
    {"BinomialCoefficientFunctor", R"doc(class GWFrames::BinomialCoefficientFunctor
==========================================
  Storage for pre-calculated binomial values.

  The object is a functor for the binomial coefficients. E.g.,
  BinomialCoefficientFunctor Binom;
const double Binom_4_2 = Binom(4,2);

  Note that values up to Binom(2*ellMax_GWFrames, ...) are calculated, where
  ellMax_GWFrames is set in Utilities.hpp. Requests for greater values are
  answered with an exception.)doc"},
    {"Boost", R"doc(Parameters
----------
  ThreeVector v
  ThreeVector n

Returns
-------
  Quaternions::Quaternion)doc"},
    {"BoostRotor", R"doc(Return a rotor taking the given vector into its boosted version.
================================================================
  Parameters
  ----------
    vector<double> ThreeVelocity
      Three-vector velocity of the new frame WRT this frame
    vector<double> ThreeVectorToBeBoosted
      Three-vector direction to be boosted by the rotor

  Returns
  -------
    Quaternion

  Description
  -----------
    This function returns a rotor $R_b$ that takes the input vector $\hat{w}$
    (which will be normalized)  considered as a direction on the future null
    sphere  into its boosted version. Note that this rotor is a function of
    both the vector being boosted and the boost itself.)doc"},
    {"ComplexI", R"doc(Parameters
----------
  0. 0
  1. 0

Returns
-------
  const complex<double>)doc"},
    {"Component", R"doc(Parameters
----------
  const vector<Quaternions::Quaternion>& Q
  const unsigned int i

Returns
-------
  vector<double>)doc"},
    {"Component0", R"doc(Parameters
----------
  const vector<Quaternions::Quaternion>& Q

Returns
-------
  vector<double>)doc"},
    {"Component1", R"doc(Parameters
----------
  const vector<Quaternions::Quaternion>& Q

Returns
-------
  vector<double>)doc"},
    {"Component2", R"doc(Parameters
----------
  const vector<Quaternions::Quaternion>& Q

Returns
-------
  vector<double>)doc"},
    {"Component3", R"doc(Parameters
----------
  const vector<Quaternions::Quaternion>& Q

Returns
-------
  vector<double>)doc"},
    {"ConformalFactorGrid", R"doc(Construct a grid with the conformal factor at each point.
=========================================================
  Parameters
  ----------
    const ThreeVector& v
    const int n_theta
    const int n_phi

  Returns
  -------
    DataGrid)doc"},
    {"DataGrid", R"doc(class GWFrames::DataGrid
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
    vector<complex<double>> data)doc"},
    {"DataGrid::Data", R"doc(Parameters
----------
  (none)

Returns
-------
  vector<complex<double>>)doc"},
    {"DataGrid::DataGrid", R"doc(  Parameters
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
    alignment of the point to be evaluated. In particular, the rotor takes the
    $\hat{z}$ vector into the point at which the field is to be measured, and
    takes $\hat{x} + i \hat{y}$ into the $m$ vector (within normalization)
    needed for spin-weighted fields.

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
  DataGrid)doc"},
    {"DataGrid::N_phi", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"DataGrid::N_theta", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"DataGrid::SetNPhi", R"doc(Parameters
----------
  const int N_phi

Returns
-------
  DataGrid&)doc"},
    {"DataGrid::SetNTheta", R"doc(Parameters
----------
  const int N_theta

Returns
-------
  DataGrid&)doc"},
    {"DataGrid::SetSpin", R"doc(Parameters
----------
  const int ess

Returns
-------
  DataGrid&)doc"},
    {"DataGrid::Spin", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"DataGrid::operator*", R"doc(Parameters
----------
  const DataGrid& 

Returns
-------
  DataGrid)doc"},
    {"DataGrid::operator+", R"doc(Parameters
----------
  const DataGrid& 

Returns
-------
  DataGrid)doc"},
    {"DataGrid::operator-", R"doc(Parameters
----------
  const DataGrid& 

Returns
-------
  DataGrid)doc"},
    {"DataGrid::operator/", R"doc(Parameters
----------
  const DataGrid& 

Returns
-------
  DataGrid)doc"},
    {"DataGrid::operator[]", R"doc(Parameters
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
  complex<double>&)doc"},
    {"DataGrid::pow", R"doc(Parameters
----------
  const int p

Returns
-------
  DataGrid)doc"},
    {"DataGrid::size", R"doc(Parameters
----------
  (none)

Returns
-------
  unsigned int)doc"},
    {"DifferentiateRotorByLogarithm", R"doc(Calculate the derivative of a rotor by the logarithm formula.
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
    the correctness of the logarithm formula.)doc"},
    {"FactorialFunctor", R"doc(class GWFrames::FactorialFunctor
================================
  Storage for pre-calculated factorials.

  Rather than re-calculating factorial values each time they are needed, we
  simply cache the values, which are pre-calculated when this object is
  constructed. They may then be accessed by querying this object as a functor.
  E.g., FactorialFunctor Factorial;
const double Fac17 = Factorial(17);

  Note that only values up to Factorial(171) are calculated, because that is
  the largest value that can be stored as an int in standard C++.)doc"},
    {"FrameAngularVelocity", R"doc(Parameters
----------
  const vector<Quaternion>& f
  const vector<double>& t

Returns
-------
  vector<Quaternion>)doc"},
    {"FrameFromAngularVelocity", R"doc(Find the frame with the given angular velocity data.
====================================================
  Parameters
  ----------
    const vector<Quaternion>& Omega
      Vector of Quaternions.
    const vector<double>& T
      Vector of corresponding times.

  Returns
  -------
    vector<Quaternions::Quaternion>

  Description
  -----------
    Note that each element of Omega should be a pure-vector Quaternion,
    corresponding to the angular-velocity vector at the instant of time.

    This function is not as accurate as might be hoped, because the time step
    of the input data is not adjustable, so you can't do anything like
    'adaptive' integration. If you have a function that returns Omega given
    time, you should use the other version of this (overloaded) function that
    uses such a function. Alternatively, if you are doing an integration as
    part of a larger system, you can use the code for the other version of this
    function as a guide on how to do so.

    FrameFromAngularVelocity(std::vector<double> (* Omega)(const double t),
    const double t0, const double t1, std::vector<Quaternion>& Qs,
    std::vector<double>& Ts)


Find the frame with the given angular velocity function.
========================================================
  Parameters
  ----------
    vector<double>(*)(const double t) Omega
      Function pointer returning angular velocity
    const double t0
      Initial time
    const double t1
      Final time
    vector<Quaternion>& Qs
      Output frame rotors
    vector<double>& Ts
      Output time steps

  Returns
  -------
    void

  Description
  -----------
    This function takes a function pointer Omega (which returns a 3-vector,
    given the time) and integrates to find the frame with that angular velocity.

    This function may not be very useful in general, because the angular
    velocity may not be known as a function of time. However, there are
    situations where the angular velocity is known at an instant of time, given
    other information. The code for this function should serve as a useful
    guide when implementing such integrations.

    In particular, the key piece in this integration is to reset the value of
    the quaternion logarithm (denoted below as r) between integration steps
    when the magnitude of r is too large. It gets reset to a value that is
    identical in terms of the resulting rotation, but has a smaller magnitude,
    so that the final result doesn't wander too much. This is equivalent to
    changing branches of a complex logarithm.

    There are two important things to note about this resetting procedure.
    First, the time stepper may want to take a very small step immediately
    after the reset, and should not be cause for alarm. Below, this is dealt
    with by also restting nSteps, and making sure that we take at least 10 more
    steps after that to let the time stepper adjust its step sizes accordingly.

    Second, the resulting rotor (which is the exponential of the logarithm)
    will flip signs when the logarithm is reset. This will have no effect on
    the physical frame deduced from the rotor, but could be bad news for
    interpolations. So we simply 'unflip' the signs when returning.)doc"},
    {"FrameFromAngularVelocity_2D", R"doc(Find the frame with the given angular velocity function.
========================================================
  Parameters
  ----------
    vector<double>(*)(const double t) Omega
      Function pointer returning angular velocity
    const double t0
      Initial time
    const double t1
      Final time
    vector<Quaternion>& Qs
      Output frame rotors
    vector<double>& Ts
      Output time steps

  Returns
  -------
    void

  Description
  -----------
    This function takes a function pointer Omega (which returns a 3-vector,
    given the time) and integrates to find the frame with that angular velocity.

    This function may not be very useful in general, because the angular
    velocity may not be known as a function of time. However, there are
    situations where the angular velocity is known at an instant of time, given
    other information. The code for this function should serve as a useful
    guide when implementing such integrations.

    In particular, the key piece in this integration is to reset the value of
    the quaternion logarithm (denoted below as r) between integration steps
    when the magnitude of r is too large. It gets reset to a value that is
    identical in terms of the resulting rotation, but has a smaller magnitude,
    so that the final result doesn't wander too much. This is equivalent to
    changing branches of a complex logarithm.

    There are two important things to note about this resetting procedure.
    First, the time stepper may want to take a very small step immediately
    after the reset, and should not be cause for alarm. Below, this is dealt
    with by also restting nSteps, and making sure that we take at least 10 more
    steps after that to let the time stepper adjust its step sizes accordingly.

    Second, the resulting rotor (which is the exponential of the logarithm)
    will flip signs when the logarithm is reset. This will have no effect on
    the physical frame deduced from the rotor, but could be bad news for
    interpolations. So we simply 'unflip' the signs when returning.

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
  corresponding to the angular-velocity vector at the instant of time.)doc"},
    {"FrameFromAngularVelocity_2D_Integrand", R"doc(Time-derivative of 2-D quaternion logarithm for vector with given angular velocity.
===================================================================================
  Parameters
  ----------
    const double r_x
    const double r_y
    vector<double> Omega
    double& rdot_x
    double& rdot_y

  Returns
  -------
    void)doc"},
    {"FrameFromAngularVelocity_Integrand", R"doc(Time-derivative of the quaternion logarithm for frame with given angular velocity.
==================================================================================
  Parameters
  ----------
    const vector<double>& r
      Quaternion logarithm corresponding to the current frame
    const vector<double>& Omega
      Quaternion vector giving instantaneous orbital velocity of frame

  Returns
  -------
    vector<double>

  Description
  -----------
    This function returns the time-derivative of the quaternion logarithm of
    the frame with angular velocity Omega. This can be integrated to give the
    quaternion logarithm of the frame, which can then be exponentiated to give
    the frame rotor.)doc"},
    {"FrameFromPrescribedRotation", R"doc(Construct minimal-rotation frame from Z basis vector of that frame.
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
    in the XY plane.)doc"},
    {"FrameFromXY", R"doc(Construct frame given the X and Y basis vectors of that frame.
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
    (x,y,z) vectors into the new frame's (X,Y,Z) vectors.)doc"},
    {"FrameFromZ", R"doc(Construct minimal-rotation frame from Z basis vector of that frame.
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
    initial rotation about z.)doc"},
    {"GWFrames", R"doc(namespace GWFrames
==================
  Collection of objects and functions for manipulating gravitational waves.)doc"},
    {"ImaginaryI", R"doc(Parameters
----------
  0. 0
  1. 0

Returns
-------
  const complex<double>)doc"},
    {"IniLIGO_Approx", R"doc(Parameters
----------
  const vector<double>& F
  const bool Invert = false
  const double NoiseFloor = 0.0

Returns
-------
  vector<double>)doc"},
    {"Interpolator::Interpolator", R"doc(Parameters
----------
  const vector<double>& x
  const vector<double>& y
  int m

Returns
-------
  Interpolator)doc"},
    {"Interpolator::hunt", R"doc(Parameters
----------
  const double x

Returns
-------
  int)doc"},
    {"Interpolator::locate", R"doc(Parameters
----------
  const double x

Returns
-------
  int)doc"},
    {"InverseConformalFactorBoostedGrid", R"doc(Construct a boosted grid with the conformal factor at each point.
=================================================================
  Parameters
  ----------
    const ThreeVector& v
    const int n_theta
    const int n_phi

  Returns
  -------
    DataGrid)doc"},
    {"InverseConformalFactorFunctor", R"doc(class InverseConformalFactorFunctor
===================================
  Member variables
  ----------------
    double gamma
    Quaternion v)doc"},
    {"InverseConformalFactorFunctor::InverseConformalFactorFunctor", R"doc(Parameters
----------
  const ThreeVector& vi

Returns
-------
  InverseConformalFactorFunctor)doc"},
    {"InverseConformalFactorFunctor::operator()", R"doc(Parameters
----------
  const Quaternions::Quaternion& R

Returns
-------
  double)doc"},
    {"InverseConformalFactorGrid", R"doc(Construct a grid with the conformal factor at each point.
=========================================================
  Parameters
  ----------
    const ThreeVector& v
    const int n_theta
    const int n_phi

  Returns
  -------
    DataGrid)doc"},
    {"LadderOperatorFactorFunctor", R"doc(class GWFrames::LadderOperatorFactorFunctor
===========================================
  Storage for pre-calculated values of the ladder operator factor.

  The object is a functor for the coefficients of the standard angular-momentum
  ladder operator. E.g., LadderOperatorFactorFunctor LadderCoeff;
const double
  Ladder_4_3 = LadderCoeff(4,3);

  Note that values up to LadderCoeff(ellMax_GWFrames, ...) are calculated,
  where ellMax_GWFrames is set in Utilities.hpp. Requests for greater values
  are answered with an exception.)doc"},
    {"MeanRotor", R"doc(Mean rotor, minimizing integrated distance between R(t_i) and R_mean.
=====================================================================
  Parameters
  ----------
    const vector<Quaternion>& R
      Vector of rotors to be averaged.

  Returns
  -------
    Quaternion

  Description
  -----------
    This algorithm minimizes the integrated distance between the input data's
    R(t) and R_mean by varying R_mean and finding the minimum using the
    Nedler-Mead simplex algorithm(2) of GSL. This is the correct generalization
    of the similar ApproximateMeanRotor function, and is essentially described
    in my paper http://arxiv.org/abs/1302.2919/.


Mean rotor, minimizing integrated distance between R(t) and R_mean.
===================================================================
  Parameters
  ----------
    const vector<Quaternion>& R
      Vector of rotors to be averaged.
    const vector<double>& t
      Vector of corresponding time steps.
    double t1 = -1e300
      Optional early bound of time to average over
    double t2 = 1e300
      Optional late bound of time to average over

  Returns
  -------
    Quaternion

  Description
  -----------
    This algorithm minimizes the integrated distance between the input data's
    R(t) and R_mean by varying R_mean and finding the minimum using the
    Nedler-Mead simplex algorithm(2) of GSL. This is the correct generalization
    of the similar ApproximateMeanRotor function, and is essentially described
    in my paper http://arxiv.org/abs/1302.2919/.)doc"},
    {"MinimalRotation", R"doc(Minimal-rotation version of the input frame.
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
    the minimal-rotation condition.)doc"},
    {"Modes", R"doc(class GWFrames::Modes
=====================
  Member variables
  ----------------
    int s

            This object holds complex spin-weighted data on the sphere in a
      spin-weighted spherical-harmonic mode representation. All modes are
      present, even for $\ell<|s|$. The modes are also assumed to be stored in
      order, as $(\ell,m) = (0,0), (1,-1), (1,0), (1,1), (2,-2), \ldots$.    int ellMax
    vector<complex<double>> data)doc"},
    {"Modes::Data", R"doc(Parameters
----------
  (none)

Returns
-------
  vector<complex<double>>)doc"},
    {"Modes::EllMax", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"Modes::EvaluateAtPoint", R"doc(Evaluate Waveform at a particular sky location.
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
      Quaternion giving point by rotation of $\hat{z}$

  Returns
  -------
    complex<double>

  Description
  -----------
    Note that the argument R might typically be thought of as the rotor taking
    the unit $z$ vector into a point $(\vartheta, \varphi)$. However, more
    general arguments are possible; this feature is used to greatly simplify
    the process of constructing a DataGrid object from a Modes object with a
    boost.)doc"},
    {"Modes::Modes", R"doc(Parameters
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
  Modes)doc"},
    {"Modes::SetEllMax", R"doc(Parameters
----------
  const int ell

Returns
-------
  Modes&)doc"},
    {"Modes::SetSpin", R"doc(Parameters
----------
  const int ess

Returns
-------
  Modes&)doc"},
    {"Modes::Spin", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"Modes::bar", R"doc(Parameters
----------
  (none)

Returns
-------
  Modes)doc"},
    {"Modes::edth", R"doc(Parameters
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
  is identical (up to a factor of $\sqrt{2}$) to the NP edth. There is an
  additional term in the definition of the GHP operator, but its value is
  zero. (It transforms nontrivially, though.) In this context, we have
  NPEdth() = sqrt(2)*GHPEdth().

  The complex shear $\sigma$ has spin weight +2 and boost weight +1. The
  radial coordinate $r$ has boost weight -1, and the derivative with respect
  to time $d/du$ has boost weight -1. The asymptotic metric shear $r\, h$
  has spin weight -2 and boost weight -1. In particular, it seems that $r\,
  h = r^2\, \bar{\sigma}$.

  The Newman-Penrose scalars $\Psi_i$ have spin weight and boost weight
  equal to $2-i$. (E.g., $\Psi_4$ has $s = b = -2$.) However, when these are
  multiplied by the appropriate factors of $r$ to find the leading-order
  terms, they acquire boost weights. In particular, we need to multiply
  $\Psi_i$ by $r^{5-i}$ to get nonzero values at scri, which adds $i-5$ to
  the boost weight, so that the asymptotic NP scalars all have boost weight
  -3.)doc"},
    {"Modes::edth2edthbar2", R"doc(The operator edth^2 bar{edth}^2.
================================
  Parameters
  ----------
    (none)

  Returns
  -------
    Modes)doc"},
    {"Modes::edthbar", R"doc(Parameters
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
  is identical (up to a factor of $\sqrt{2}$) to the NP edth. There is an
  additional term in the definition of the GHP operator, but its value is
  zero. (It transforms nontrivially, though.) In this context, we have
  NPEdthBar() = sqrt(2)*GHPEdthBar().

  The complex shear $\sigma$ has spin weight +2 and boost weight +1. The
  radial coordinate $r$ has boost weight -1, and the derivative with respect
  to time $d/du$ has boost weight -1. The asymptotic metric shear $r\, h$
  has spin weight -2 and boost weight -1. In particular, it seems that $r\,
  h = r^2\, \bar{\sigma}$.

  The Newman-Penrose scalars $\Psi_i$ have spin weight and boost weight
  equal to $2-i$. (E.g., $\Psi_4$ has $s = b = -2$.) However, when these are
  multiplied by the appropriate factors of $r$ to find the leading-order
  terms, they acquire boost weights. In particular, we need to multiply
  $\Psi_i$ by $r^{5-i}$ to get nonzero values at scri, which adds $i-5$ to
  the boost weight, so that the asymptotic NP scalars all have boost weight
  -3.)doc"},
    {"Modes::operator*", R"doc(Parameters
----------
  const Modes& M

Returns
-------
  Modes)doc"},
    {"Modes::operator+", R"doc(Parameters
----------
  const Modes& M

Returns
-------
  Modes)doc"},
    {"Modes::operator-", R"doc(Parameters
----------
  const Modes& M

Returns
-------
  Modes)doc"},
    {"Modes::operator/", R"doc(Parameters
----------
  const Modes& M

Returns
-------
  Modes)doc"},
    {"Modes::operator=", R"doc(Parameters
----------
  const Modes& B

Returns
-------
  Modes&)doc"},
    {"Modes::operator[]", R"doc(Parameters
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
  complex<double>&)doc"},
    {"Modes::pow", R"doc(Parameters
----------
  const int p

Returns
-------
  Modes)doc"},
    {"Modes::size", R"doc(Parameters
----------
  (none)

Returns
-------
  unsigned int)doc"},
    {"OptimalAlignment", R"doc(Time and rotor offsets such that Rdelta*Rb(t+deltat) is as close to Ra(t) as possible.
======================================================================================
  Parameters
  ----------
    const double t1
      Initial bound of alignment range
    const double t2
      Final bound of alignment range
    const vector<Quaternion>& Ra
      Vector of fixed rotors.
    const vector<double>& ta
      Vector of corresponding time steps.
    const vector<Quaternion>& Rb
      Vector of rotors to be aligned to Ra.
    const vector<double>& tb
      Vector of corresponding time steps.
    double& deltat
      Output optimal time offset
    Quaternions::Quaternion& R_delta
      Output optimal rotor

  Returns
  -------
    void

  Description
  -----------
    This routine actively optimizes the time offset and all three rotational
    degrees of freedom for Rb, so that Rdelta*Rb(t+deltat) is as close to Ra(t)
    as possible over the range (t1, t2). The error measure is just the usual
    angle(UnflipRotors(RDelta(Ra,Rb))) integrated over (t1,t2).

    The times t1 and t2 are measured relative to ta, and all are left fixed; Rb
    and tb are shifted to achieve alignment.)doc"},
    {"PNWaveform", R"doc(class GWFrames::PNWaveform
==========================
  Object for calculating a post-Newtonian Waveform with (optional) precession.

  Fundamental object creating precessing PN waveform.

  This object is a subclass of the GWFrames::Waveform object. In addition to
  the data stored in a Waveform, this stores the two spins chi1 and chi2; the
  orbital angular-velocity vectors Omega_orb, Omega_prec, and Omega_tot; the PN
  angular momentum; and the orbital phase Phi_orb. Various methods also exist
  for retrieving the vectors, their magnitudes, and their normalized versions.

  history

  String recording history of data

  t

  Vector of time indices

  frame

  Vector of quaternions. This represents the rotation applied to the physical
  system to go from the standard inertial frame to the frame in which the
  Waveform is decomposed. If this vector has length 0, the frame is assumed to
  be constant, and the standard frame. If this vector has length 1, the frame
  is assumed to be constant, but given by that value. Otherwise, this vector is
  assumed to have length NTimes(), corresponding to the frame at each moment.

  lm

  Vector of vectors of int's, giving (ell,m) data

  data

  Vector of complex numbers representing h or Psi4.

  chi1

  Vector of vectors of doubles representing the spin of BH 1.

  chi2

  Vector of vectors of doubles representing the spin of BH 2.

  Omega_orb

  Vector of vectors of doubles representing the orbital angular-velocity vector.

  Omega_prec

  Vector of vectors of doubles representing the precessional angular-velocity
  vector.

  Omega_tot

  Vector of vectors of doubles representing the total angular-velocity vector.

  L

  Vector of vectors of doubles representing the angular-momentum vector.

  Phi_orb

  Vector of doubles representing the orbital phase.

  Member variables
  ----------------
    vector<vector<double>> mchi1
    vector<vector<double>> mchi2
    vector<vector<double>> mOmega_orb
    vector<vector<double>> mOmega_prec
    vector<vector<double>> mL
    vector<double> mPhi_orb)doc"},
    {"PNWaveform::L", R"doc(Parameters
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
  const vector<vector<double>>&)doc"},
    {"PNWaveform::LHat", R"doc(Parameters
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
  vector<vector<double>>)doc"},
    {"PNWaveform::LMag", R"doc(  Parameters
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
    vector<double>)doc"},
    {"PNWaveform::OmegaHat_orb", R"doc(Parameters
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
  vector<vector<double>>)doc"},
    {"PNWaveform::OmegaHat_prec", R"doc(Parameters
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
  vector<vector<double>>)doc"},
    {"PNWaveform::OmegaHat_tot", R"doc(Parameters
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
  vector<vector<double>>)doc"},
    {"PNWaveform::Omega_orb", R"doc(Parameters
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
  const vector<vector<double>>&)doc"},
    {"PNWaveform::Omega_orbMag", R"doc(  Parameters
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
    vector<double>)doc"},
    {"PNWaveform::Omega_prec", R"doc(Parameters
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
  const vector<vector<double>>&)doc"},
    {"PNWaveform::Omega_precMag", R"doc(  Parameters
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
    vector<double>)doc"},
    {"PNWaveform::Omega_tot", R"doc(Total angular velocity of PN binary at an instant of time.
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
    vector<vector<double>>)doc"},
    {"PNWaveform::Omega_totMag", R"doc(  Parameters
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
    vector<double>)doc"},
    {"PNWaveform::PNWaveform", R"doc(Default constructor for an empty object.
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
    const double delta
    const vector<double>& chi1_i
    const vector<double>& chi2_i
    const double Omega_orb_i
    double Omega_orb_0 = -1.0
    const Quaternions::Quaternion& R_frame_i = Quaternions::Quaternion(1, 0, 0, 0)
    const unsigned int MinStepsPerOrbit = 32
    const double PNWaveformModeOrder = 3.5
    const double PNOrbitalEvolutionOrder = 4.0

  Returns
  -------
    PNWaveform

  Description
  -----------
    See GWFrames/Code/SWIG/Extensions.py for the docstring for this object)doc"},
    {"PNWaveform::Phi_orb", R"doc(Parameters
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
  const vector<double>&)doc"},
    {"PNWaveform::chi1", R"doc(Parameters
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
  const vector<vector<double>>&)doc"},
    {"PNWaveform::chi1Mag", R"doc(Parameters
----------
  const unsigned int iTime

Returns
-------
  double)doc"},
    {"PNWaveform::chi2", R"doc(Parameters
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
  const vector<vector<double>>&)doc"},
    {"PNWaveform::chi2Mag", R"doc(Parameters
----------
  const unsigned int iTime

Returns
-------
  double)doc"},
    {"PNWaveform::chiHat1", R"doc(Parameters
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
  vector<vector<double>>)doc"},
    {"PNWaveform::chiHat2", R"doc(Parameters
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
  vector<vector<double>>)doc"},
    {"PNWaveform::~PNWaveform", R"doc(Parameters
----------
  (none)

Returns
-------
  ~PNWaveform)doc"},
    {"PolynomialInterpolator::rawinterp", R"doc(Parameters
----------
  int jl
  double x

Returns
-------
  double)doc"},
    {"PrescribedRotation", R"doc(Input frame with prescribed rate of rotation about Z axis.
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
    about the new Z axis to agree with the given rotation rate.)doc"},
    {"Quaternion", R"doc(class Quaternions::Quaternion
=============================
  Object representing an individual quaternion.

  Member variables
  ----------------
    double w
    double x
    double y
    double z)doc"},
    {"Quaternion::Quaternion", R"doc(Empty constructor  initialized to 0s.
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
    rotation about the given axis through the given angle.)doc"},
    {"Quaternion::abs", R"doc(Parameters
----------
  (none)

Returns
-------
  double)doc"},
    {"Quaternion::angle", R"doc(Parameters
----------
  (none)

Returns
-------
  double)doc"},
    {"Quaternion::commutator", R"doc(Parameters
----------
  const Quaternion& Q

Returns
-------
  Quaternion)doc"},
    {"Quaternion::conjugate", R"doc(Parameters
----------
  (none)

Returns
-------
  Quaternion)doc"},
    {"Quaternion::cross", R"doc(Parameters
----------
  const Quaternion& Q

Returns
-------
  Quaternion)doc"},
    {"Quaternion::dot", R"doc(Parameters
----------
  const Quaternion& Q

Returns
-------
  double)doc"},
    {"Quaternion::exp", R"doc(Return exponent of Quaternion.
==============================
  Parameters
  ----------
    (none)

  Returns
  -------
    Quaternion)doc"},
    {"Quaternion::inverse", R"doc(Parameters
----------
  (none)

Returns
-------
  Quaternion)doc"},
    {"Quaternion::log", R"doc(Return logarithm of Quaternion.
===============================
  Parameters
  ----------
    (none)

  Returns
  -------
    Quaternion)doc"},
    {"Quaternion::logRotor", R"doc(Return logarithm of a rotor.
============================
  Parameters
  ----------
    (none)

  Returns
  -------
    Quaternion

  Description
  -----------
    This function is just like the standard log function, except that a
    negative scalar does not raise an exception; instead, the scalar pi is
    returned.)doc"},
    {"Quaternion::normalized", R"doc(Parameters
----------
  (none)

Returns
-------
  Quaternion)doc"},
    {"Quaternion::normsquared", R"doc(Parameters
----------
  (none)

Returns
-------
  double)doc"},
    {"Quaternion::operator!=", R"doc(Parameters
----------
  const Quaternion& Q

Returns
-------
  bool)doc"},
    {"Quaternion::operator*", R"doc(  Parameters
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
    Quaternion)doc"},
    {"Quaternion::operator+", R"doc(Parameters
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
  Quaternion)doc"},
    {"Quaternion::operator+=", R"doc(Parameters
----------
  const Quaternion& Q

Returns
-------
  Quaternion&)doc"},
    {"Quaternion::operator-", R"doc(Parameters
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
  Quaternion)doc"},
    {"Quaternion::operator/", R"doc(Parameters
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
  Quaternion)doc"},
    {"Quaternion::operator=", R"doc(Parameters
----------
  const Quaternion& Q

Returns
-------
  Quaternion&)doc"},
    {"Quaternion::operator==", R"doc(Parameters
----------
  const Quaternion& Q

Returns
-------
  bool)doc"},
    {"Quaternion::operator[]", R"doc(Get component of Quaternion.
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
    Note: This is unavailable from python.)doc"},
    {"Quaternion::pow", R"doc(Parameters
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
  Quaternion)doc"},
    {"Quaternion::sqrt", R"doc(Parameters
----------
  (none)

Returns
-------
  Quaternion)doc"},
    {"Quaternion::sqrtOfRotor", R"doc(Parameters
----------
  (none)

Returns
-------
  Quaternion)doc"},
    {"Quaternion::vec", R"doc(Parameters
----------
  (none)

Returns
-------
  vector<double>)doc"},
    {"QuaternionArray", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"QuaternionDerivative", R"doc(Five-point finite-differencing of vector of Quaternions.
========================================================
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
    The formula for this finite difference comes from Eq. (A 5b) of 'Derivative
    formulas and errors for non-uniformly spaced points' by M. K. Bowen and
    Ronald Smith. As explained in their Eqs. (B 9b) and (B 10b), this is a
    fourth-order formula.

    If there are fewer than five points, the function reverts to simpler
    formulas. If there are fewer than two points, or there are different
    numbers of points in the two input vectors, an exception is thrown.

    Sundqvist and Veronis, Tellus XXII (1970), 1)doc"},
    {"Quaternions", R"doc(namespace Quaternions
=====================)doc"},
    {"RDelta", R"doc(Difference between frame rotors R1*R2.inverse()
===============================================
  Parameters
  ----------
    const vector<Quaternion>& R1
      Vector of rotors
    const vector<Quaternion>& R2
      Vector of rotors
    const int IndexOfFiducialTime = -1
      Integer index of time at which difference is set to zero [default: -1]

  Returns
  -------
    vector<Quaternion>

  Description
  -----------
    If the optional IndexOfFiducialTime is given, the returned quantity is
    R1*Offset*R2.inverse(), where Offset is R1[I].inverse()*R2[I], thus
    ensuring that RDelta[I]=1.)doc"},
    {"SQR", R"doc(Parameters
----------
  const double a

Returns
-------
  double)doc"},
    {"SWSH", R"doc(class GWFrames::SWSH
====================
  Calculate SWSHs directly from Wigner D matrices.

  This object calculates values of SWSHs directly from the Wigner D matrices,
  in terms of rotors, with the conventions described in
  http://arxiv.org/abs/0709.0093. Here, the value at some point (theta, phi) is
  given by evaluating this object for a rotor that takes the north pole (0,0,1)
  into the point (theta, phi).)doc"},
    {"ScalarIntegral", R"doc(Integrate scalar function by simple trapezoidal rule.
=====================================================
  Parameters
  ----------
    const vector<double>& fdot
      Vector of scalars.
    const vector<double>& t
      Vector of corresponding time steps.

  Returns
  -------
    vector<double>)doc"},
    {"Scri", R"doc(class GWFrames::Scri
====================
  Member variables
  ----------------
    vector<double> t

            A Scri object contains all the information needed to express the
      asymptotic geometry of null infinity, and symmetry transformations
      thereof. This geometry is encoded in the NewmanPenrose curvature scalars
      $\psi_0, \ldots, \psi_4$ and the complex shear $\sigma$ of outgoing
      null rays (strain in gravitational-wave detectors). The general symmetry
      transformation is an element of the BondiMetznerSachs (BMS) group, which
      transforms the data contained by Scri among itself.    vector<SliceModes> slices)doc"},
    {"Scri::BMSTransformation", R"doc(Apply a (constant) BMS transformation to data on null infinity.
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
    K(u-\delta)$, where $K$ is the conformal factor (which is a function of
    angle). So we need to interpolate at each point to the original time
    coordinate $u = u'/K + \delta$, which again depends on angle. That is, we
    have to interpolate to a different time for each grid point. In this case,
    we arbitrarily set $u' = 0$, because any other choice can be absorbed into
    a time- and space-translation. This does not matter, of course, because
    that choice is not stored in any way.)doc"},
    {"Scri::NTimes", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"Scri::Scri", R"doc(Parameters
----------
  const Scri& S

Returns
-------
  Scri

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
  Scri)doc"},
    {"Scri::T", R"doc(Parameters
----------
  (none)

Returns
-------
  const vector<double>)doc"},
    {"Scri::operator[]", R"doc(Parameters
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
  SliceModes&)doc"},
    {"ScriFunctor", R"doc(class GWFrames::ScriFunctor
===========================)doc"},
    {"ScriFunctor::operator()", R"doc(Parameters
----------
  const Quaternions::Quaternion& 

Returns
-------
  double)doc"},
    {"Slerp", R"doc(Parameters
----------
  const double tau
  const Quaternion& Qa
  const Quaternion& Qb

Returns
-------
  Quaternion)doc"},
    {"SliceModes", R"doc(class GWFrames::SliceModes
==========================)doc"},
    {"SliceModes::BMSTransformationOnSlice", R"doc(Execute a BMS transformation except for the supertranslation of points.
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
    to interpolate to the supertranslated time.)doc"},
    {"SliceModes::EllMax", R"doc(Find largest ell value in the data on this slice.
=================================================
  Parameters
  ----------
    (none)

  Returns
  -------
    int)doc"},
    {"SliceModes::FourMomentum", R"doc(Calculate the four-momentum of the system from the supermomentum.
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
    supermomentum.)doc"},
    {"SliceModes::Mass", R"doc(Calculate the mass of the system from the four-momentum.
========================================================
  Parameters
  ----------
    (none)

  Returns
  -------
    double)doc"},
    {"SliceModes::MoreschiIteration", R"doc(Find the next iteration of the BMS transformation via Moreschi's algorithm.
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
    transformed by the BMS transformation represented by $K_i$ and $\delta_i$.
    This then takes the data on that slice and computes the values of $K_{i+1}$
    and $\delta_{i+1}$, returning them by reference.)doc"},
    {"SliceModes::SliceModes", R"doc(Parameters
----------
  const SliceModes& S

Returns
-------
  SliceModes

Constructor from ellMax.
========================
  Parameters
  ----------
    const int ellMax = 0

  Returns
  -------
    SliceModes)doc"},
    {"SliceModes::SuperMomentum", R"doc(Find the Moreschi supermomentum.
================================
  Parameters
  ----------
    (none)

  Returns
  -------
    Modes

  Description
  -----------
    $\psi = \psi_2 + \sigma \dot{\bar{\sigma}} + \eth^2 \bar{\sigma}$)doc"},
    {"SliceOfScri", R"doc(class GWFrames::SliceOfScri
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
    D sigmadot)doc"},
    {"SliceOfScri::SliceOfScri", R"doc(Empty constructor with reserved storage.
========================================
  Parameters
  ----------
    const int size = 0

  Returns
  -------
    SliceOfScri)doc"},
    {"SliceOfScri<D>::SliceOfScri", R"doc(Parameters
----------
  const SliceOfScri& S

Returns
-------
  SliceOfScri)doc"},
    {"SliceOfScri<D>::operator[]", R"doc(Parameters
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
  D&)doc"},
    {"SplineIntegrator::CumulativeIntegral", R"doc(Parameters
----------
  (none)

Returns
-------
  double)doc"},
    {"SplineIntegrator::operator()", R"doc(Parameters
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
  vector<double>)doc"},
    {"SplineInterpolator::rawinterp", R"doc(Parameters
----------
  int jl
  double xv

Returns
-------
  double)doc"},
    {"SplineInterpolator::sety2", R"doc(Parameters
----------
  const vector<double>& xv
  const vector<double>& yv
  double yp1
  double ypn

Returns
-------
  void)doc"},
    {"Squad", R"doc(Squad interpolation of Quaternion time series.
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
    rotations than other forms of interpolation.)doc"},
    {"StringForm", R"doc(Parameters
----------
  const vector<int>& Lmodes

Returns
-------
  string)doc"},
    {"SuperMomenta", R"doc(class GWFrames::SuperMomenta
============================
  Member variables
  ----------------
    vector<double> t
    vector<Modes> Psi)doc"},
    {"SuperMomenta::BMSTransform", R"doc(Return value of Psi on u'=const slice centered at delta[0].
===========================================================
  Parameters
  ----------
    const Modes& OneOverK
    const Modes& delta

  Returns
  -------
    Modes)doc"},
    {"SuperMomenta::MoreschiIteration", R"doc(Transform to given slice with given BMS transformation, and return next step in Moreschi algorithm.
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
    algorithm.)doc"},
    {"SuperMomenta::NTimes", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"SuperMomenta::SuperMomenta", R"doc(Parameters
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
  SuperMomenta)doc"},
    {"SuperMomenta::T", R"doc(Parameters
----------
  (none)

Returns
-------
  const vector<double>)doc"},
    {"SuperMomenta::operator[]", R"doc(Parameters
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
  Modes&)doc"},
    {"TransitionFunction_Smooth", R"doc(Local utility function.
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
    but non-analytic function.)doc"},
    {"UnflipRotors", R"doc(Remove sign-ambiguity of rotors.
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
    successive instants are within 'discont' of each other.)doc"},
    {"VectorStringForm", R"doc(Parameters
----------
  const vector<double>& V

Returns
-------
  string)doc"},
    {"Waveform", R"doc(class GWFrames::Waveform
========================
  Object storing data and other information for a single waveform.

  Fundamental object encapsulating waveform data, such as time, (l,m)
  information, and complex data.

  This object provides the main user interface for this collection of code. The
  various methods for this class are intended to provide all manipulations that
  might be necessary in the course of waveform analysis.

  history

  String recording history of data

  t

  Vector of time indices

  frame

  Vector of unit quaternions (rotors). This rotation, applied to the standard
  basis (x,y,z) results in the basis (X,Y,Z), with respect to which the
  waveform is decomposed. If this vector has length 0, the frame is assumed to
  be constant. If this vector has length 1, the frame is assumed to be
  constant, but given by that value. Otherwise, this vector is assumed to have
  length NTimes(), corresponding to the frame at each moment.

  FrameType

  Enum type defined as {UnknownFrameType, Inertial, Aligned, Coorbital,
  Corotating}. This is used to check that appropriate operations are performed,
  when such operations are sensitive to the frame.

  DataType

  Enum type defined as {UnknownDataType, h, hdot, Psi4}. This is used to check
  that appropriate operations are performed, when such operations are sensitive
  to the data type (finding the flux, for example).

  RScaledOut

  Boolean stating whether or not radial dependence has been scaled out of the
  data. If true, the data represents r*h or r*Psi4, for example.

  MScaledOut

  Boolean stating whether or not the mass dependence has been scaled out of the
  data. If true, the data represents h/M or M*Psi4, for example, and the time
  represents t/M.

  lm

  Vector of vectors of int's, giving (ell,m) data

  data

  Vector of complex numbers representing h or Psi4. Note that this data is
  stored as real and imaginary parts. In the corotating frame, this is probably
  optimal, as decomposing into modulus and phase is no longer useful, and can
  actually lead to discontinuities in m=0 modes.

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
    complex<double> XParityConjugate
    complex<double> YParityConjugate
    complex<double> ZParityConjugate
    complex<double> ParityConjugate
    complex<double> ConjugateAntipodalEvaluation
    vector<double> InvolutionViolationSquared
    Waveform Involution
    Waveform InvolutionSymmetricPart
    Waveform InvolutionAntisymmetricPart
    Waveform& TransformModesToRotatedFrame
      Rotate modes of the Waveform object.)doc"},
    {"Waveform::Abs", R"doc(  Parameters
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
    vector<vector<double>>)doc"},
    {"Waveform::AlignDecompositionFrameToModes", R"doc(Fix the attitude of the corotating frame.
=========================================
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
    or [2,3,4], for example, restricts the range of the sum.)doc"},
    {"Waveform::AngularVelocityVector", R"doc(Calculate the angular velocity of the Waveform.
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
    'Angular velocity of gravitational radiation and the corotating frame'. The
    vector is given in the (possibly rotating) mode frame (X,Y,Z), rather than
    the inertial frame (x,y,z).

    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.)doc"},
    {"Waveform::AngularVelocityVectorRelativeToInertial", R"doc(Calculate the angular velocity of the Waveform.
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
    'Angular velocity of gravitational radiation and the corotating frame'. The
    vector is given in the inertial frame (x,y,z), rather than the (possibly
    rotating) mode frame (X,Y,Z), which is unlike most other methods for
    Waveform.

    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.)doc"},
    {"Waveform::AppendHistory", R"doc(Parameters
----------
  const string& Hist

Returns
-------
  Waveform&)doc"},
    {"Waveform::Arg", R"doc(  Parameters
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
    vector<vector<double>>)doc"},
    {"Waveform::ArgUnwrapped", R"doc(  Parameters
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
    vector<vector<double>>)doc"},
    {"Waveform::BinaryOp", R"doc(Pointwise multiply this object by another Waveform object.
==========================================================
  Parameters
  ----------
    typename Op 
    const Waveform& B

  Returns
  -------
    typename Op

Parameters
----------
  typename Op 
  const Waveform& b

Returns
-------
  typename Op)doc"},
    {"Waveform::BoostHFaked", R"doc(Apply a boost to Psi4 data.
===========================
  Parameters
  ----------
    const vector<vector<double>>& v

  Returns
  -------
    Waveform&

  Description
  -----------
    This function does three things. First, it evaluates the Waveform on what
    will become an equi-angular grid after transformation by the boost. Second,
    at each point of that grid, it takes the appropriate combinations of the
    present value of h and its conjugate to give the value of h as observed in
    the boosted frame. Finally, it transforms back to Fourier space using that
    new equi-angular grid.

    The input three-velocities are assumed to give the velocities of the
    boosted frame relative to the present frame.

    Note that this function simply uses the correct transformation of Psi_4,
    then multiplies by the appropriate power of gamma (-2) at each point. This,
    of course, assumes that $\ddot{h} = \Psi_4$ in both frames. That need not
    be the case, which is why 'Faked' is in the name of this function.)doc"},
    {"Waveform::BoostPsi4", R"doc(Apply a boost to Psi4 data.
===========================
  Parameters
  ----------
    const vector<vector<double>>& v

  Returns
  -------
    Waveform&

  Description
  -----------
    This function does three things. First, it evaluates the Waveform on what
    will become an equi-angular grid after transformation by the boost. Second,
    at each point of that grid, it takes the appropriate combinations of the
    present value of Psi_4 and its conjugate to give the value of Psi_4 as
    observed in the boosted frame. Finally, it transforms back to Fourier space
    using that new equi-angular grid.

    The input three-velocities are assumed to give the velocities of the
    boosted frame relative to the present frame.)doc"},
    {"Waveform::BoostWeight", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"Waveform::Compare", R"doc(Return a Waveform with differences between the two inputs.
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
    modulus and phase.)doc"},
    {"Waveform::ConjugateAntipodalAntisymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::ConjugateAntipodalEvaluation", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::ConjugateAntipodalSymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::CopyWithoutData", R"doc(Copy the Waveform, except for the data (t, frame, lm, data)
===========================================================
  Parameters
  ----------
    (none)

  Returns
  -------
    Waveform)doc"},
    {"Waveform::CorotatingFrame", R"doc(Frame in which the rotation is minimal.
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
    The result is given relative to the (possibly rotating) mode frame (X,Y,Z),
    rather than the inertial frame (x,y,z).

    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.)doc"},
    {"Waveform::Data", R"doc(  Parameters
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
    vector<vector<complex<double>>>)doc"},
    {"Waveform::DataDot", R"doc(Return time derivative of data.
===============================
  Parameters
  ----------
    const unsigned int Mode

  Returns
  -------
    vector<complex<double>>)doc"},
    {"Waveform::DataType", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"Waveform::DataTypeLaTeXString", R"doc(Parameters
----------
  (none)

Returns
-------
  string)doc"},
    {"Waveform::DataTypeString", R"doc(Parameters
----------
  (none)

Returns
-------
  string)doc"},
    {"Waveform::DescriptorString", R"doc(Parameters
----------
  (none)

Returns
-------
  string)doc"},
    {"Waveform::Differentiate", R"doc(Differentiate the waveform as a function of time.
=================================================
  Parameters
  ----------
    (none)

  Returns
  -------
    Waveform&)doc"},
    {"Waveform::DipoleMoment", R"doc(Evaluate the dipole moment of the waveform.
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
    defined as $\vec{d} = \int \hat{n} \lvert f \rvert^2 d\Omega$. Up to
    a geometric factor, this function applied to $\dot{h}$ is the rate of
    emission of momentum in gravitational waves.)doc"},
    {"Waveform::DropEllModes", R"doc(Remove data relating to the given ell modes.
============================================
  Parameters
  ----------
    const vector<unsigned int>& EllModesToDrop

  Returns
  -------
    Waveform&)doc"},
    {"Waveform::DropTimesOutside", R"doc(Remove all data relating to times outside of the given range.
=============================================================
  Parameters
  ----------
    const double ta
    const double tb

  Returns
  -------
    Waveform&)doc"},
    {"Waveform::EllMax", R"doc(Return greatest ell value present in the data.
==============================================
  Parameters
  ----------
    (none)

  Returns
  -------
    int)doc"},
    {"Waveform::EvaluateAtPoint", R"doc(Evaluate Waveform at a particular sky location.
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
    UnknownFrameType nor Inertial.)doc"},
    {"Waveform::FindModeIndex", R"doc(Find index of mode with given (l,m) data.
=========================================
  Parameters
  ----------
    const int L
    const int M

  Returns
  -------
    unsigned int)doc"},
    {"Waveform::FindModeIndexWithoutError", R"doc(Find index of mode with given (l,m) data without the chance of throwing an exception.
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
    end of the mode vector.)doc"},
    {"Waveform::Frame", R"doc(Parameters
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
  const vector<Quaternions::Quaternion>&)doc"},
    {"Waveform::FrameType", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"Waveform::FrameTypeString", R"doc(Parameters
----------
  (none)

Returns
-------
  string)doc"},
    {"Waveform::GetAlignmentOfDecompositionFrameToModes", R"doc(Find the appropriate rotation to fix the attitude of the corotating frame.
==========================================================================
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

    AlignDecompositionFrameToModes)doc"},
    {"Waveform::GetAlignmentsOfDecompositionFrameToModes", R"doc(Find the appropriate rotations to fix the attitude of the corotating frame.
===========================================================================
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
    $R_{\varepsilon}$ so that the decomposition frame is aligned to the
    waveform. This particular version finds the appropriate $R_{\varepsilon}$
    at each time in the input Waveform. This is useful in cases where we need
    to try many such alignments, because the setup for interpolation is very
    slow.

    Note that this function has no option to choose the direction of X based on
    some nHat vector, as other similar functions have. That issue is assumed to
    be handled elsewhere.)doc"},
    {"Waveform::HistoryStr", R"doc(Parameters
----------
  (none)

Returns
-------
  string)doc"},
    {"Waveform::HistoryStream", R"doc(Parameters
----------
  (none)

Returns
-------
  stringstream&)doc"},
    {"Waveform::Hybridize", R"doc(Hybridize this Waveform with another.
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

    The transition function is a $C^\infty$ function, meaning that the output
    data has exactly this Waveform's data before t1, exactly Waveform B's data
    after t2, and a smooth blend in between.

    Note that this function does NOT operate in place; a new Waveform object is
    constructed and returned.)doc"},
    {"Waveform::Im", R"doc(  Parameters
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
    vector<vector<double>>)doc"},
    {"Waveform::Interpolate", R"doc(Interpolate the Waveform to a new set of time instants.
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
    such times are requested, an error will be thrown.)doc"},
    {"Waveform::InterpolateInPlace", R"doc(Interpolate the Waveform to a new set of time instants.
=======================================================
  Parameters
  ----------
    const vector<double>& NewTime

  Returns
  -------
    Waveform&)doc"},
    {"Waveform::InterpolateToPoint", R"doc(Evaluate Waveform at a particular sky location and an instant of time.
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
    UnknownFrameType nor Inertial.)doc"},
    {"Waveform::KeepOnlyEll2", R"doc(Remove data relating to all but the ell=2 modes.
================================================
  Parameters
  ----------
    (none)

  Returns
  -------
    Waveform&)doc"},
    {"Waveform::KeepOnlyEllModes", R"doc(Remove data relating to all but the given ell modes.
====================================================
  Parameters
  ----------
    const vector<unsigned int>& EllModesToKeep

  Returns
  -------
    Waveform&)doc"},
    {"Waveform::LLDominantEigenvector", R"doc(Calculate the principal axis of the LL matrix, as prescribed by O'Shaughnessy et al.
====================================================================================
  Parameters
  ----------
    const vector<int>& Lmodes = vector<int>(0)
      L modes to evaluate (optional)
    const Quaternions::Quaternion& RoughInitialEllDirection = Quaternions::zHat
      Vague guess about the preferred initial (optional)

  Returns
  -------
    vector<vector<double>>

  Description
  -----------
    If Lmodes is empty (default), all L modes are used. Setting Lmodes to [2]
    or [2,3,4], for example, restricts the range of the sum.

    Ell is the direction of the angular velocity for a PN system, so some rough
    guess about that direction allows us to choose the direction of the
    eigenvectors output by this function to be more parallel than anti-parallel
    to that direction. The default is to simply choose the z axis, since this
    is most often the correct choice anyway.

    The vector is given in the (possibly rotating) mode frame (X,Y,Z), rather
    than the inertial frame (x,y,z).)doc"},
    {"Waveform::LLMatrix", R"doc(Calculate the $<LL>$ quantity defined in the paper.
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
    or [2,3,4], for example, restricts the range of the sum. The matrix is
    given in the (possibly rotating) mode frame (X,Y,Z), rather than the
    inertial frame (x,y,z).

    $<LL>^{ab} = \sum_{\ell,m,m'} [\bar{f}^{\ell,m'} <\ell,m' | L_a L_b |
    \ell,m> f^{\ell,m} ]$)doc"},
    {"Waveform::LM", R"doc(Parameters
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
  const vector<vector<int>>&)doc"},
    {"Waveform::LdtVector", R"doc(Calculate the $<L \partial_t>$ quantity defined in the paper.
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
    or [2,3,4], for example, restricts the range of the sum. The vector is
    given with respect to the (possibly rotating) mode frame (X,Y,Z), rather
    than the inertial frame (x,y,z).

    $<L \partial_t>^a = \sum_{\ell,m,m'} \Im [ \bar{f}^{\ell,m'}
    <\ell,m' | L_a | \ell,m> \dot{f}^{\ell,m} ]$)doc"},
    {"Waveform::MIsScaledOut", R"doc(Parameters
----------
  (none)

Returns
-------
  bool)doc"},
    {"Waveform::MaxNormIndex", R"doc(Return the data index corresponding to the time of the largest norm.
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

    MaxNormTime())doc"},
    {"Waveform::MaxNormTime", R"doc(Parameters
----------
  const unsigned int SkipFraction = 4

Returns
-------
  double)doc"},
    {"Waveform::MinimalParityViolation", R"doc(Measure the relative magnitude of the violation of parity in the z direction.
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
    and the parity violation is minimized.)doc"},
    {"Waveform::NModes", R"doc(Parameters
----------
  (none)

Returns
-------
  unsigned int)doc"},
    {"Waveform::NTimes", R"doc(Parameters
----------
  (none)

Returns
-------
  unsigned int)doc"},
    {"Waveform::Norm", R"doc(Return the norm (sum of squares of modes) of the waveform.
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

    MaxNormTime)doc"},
    {"Waveform::NormalizedAntisymmetry", R"doc(Return the normalized asymmetry as a function of time.
======================================================
  Parameters
  ----------
    vector<int> LModesForAsymmetry = vector<int>(0)
      $\ell$ modes to use when calculating numerator

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
    will always be used to calculate the normalization.)doc"},
    {"Waveform::Output", R"doc(Output Waveform object to data file.
====================================
  Parameters
  ----------
    const string& FileName
    const unsigned int precision = 14

  Returns
  -------
    const Waveform&)doc"},
    {"Waveform::ParityAntisymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::ParityInvolution", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::ParitySymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::ParityViolationNormalized", R"doc(Parameters
----------
  vector<int> Lmodes = vector<int>(0)

Returns
-------
  vector<double>)doc"},
    {"Waveform::ParityViolationSquared", R"doc(Parameters
----------
  vector<int> Lmodes = vector<int>(0)

Returns
-------
  vector<double>)doc"},
    {"Waveform::RIsScaledOut", R"doc(Parameters
----------
  (none)

Returns
-------
  bool)doc"},
    {"Waveform::Re", R"doc(  Parameters
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
    vector<vector<double>>)doc"},
    {"Waveform::ResizeData", R"doc(Parameters
----------
  const unsigned int NModes
  const unsigned int NTimes

Returns
-------
  Waveform&)doc"},
    {"Waveform::RotateDecompositionBasis", R"doc(Rotate the basis in which this Waveform is measured by a constant rotor.
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
    calling function.)doc"},
    {"Waveform::RotatePhysicalSystem", R"doc(Rotate the physical content of the Waveform by a constant rotor.
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
    original frame data times $\bar{R}_{phys}$.

    Note that this function does not change the frameType; that is left to the
    calling function.)doc"},
    {"Waveform::SetBoostWeight", R"doc(Parameters
----------
  const int NewBoostWeight

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetData", R"doc(Parameters
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
  Waveform&)doc"},
    {"Waveform::SetDataType", R"doc(Parameters
----------
  const WaveformDataType Type

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetFrame", R"doc(Parameters
----------
  const vector<Quaternions::Quaternion>& a

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetFrameType", R"doc(Parameters
----------
  const WaveformFrameType Type

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetHistory", R"doc(Parameters
----------
  const string& Hist

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetLM", R"doc(Parameters
----------
  const vector<vector<int>>& a

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetMIsScaledOut", R"doc(Parameters
----------
  const bool Scaled

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetRIsScaledOut", R"doc(Parameters
----------
  const bool Scaled

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetSpinWeight", R"doc(Parameters
----------
  const int NewSpinWeight

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetT", R"doc(Parameters
----------
  const vector<double>& a

Returns
-------
  Waveform&)doc"},
    {"Waveform::SetTime", R"doc(Parameters
----------
  const vector<double>& a

Returns
-------
  Waveform&)doc"},
    {"Waveform::SliceOfTimeIndices", R"doc(Copy the Waveform between indices i_t_a and i_t_b.
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
    necessary))doc"},
    {"Waveform::SliceOfTimeIndicesWithEll2", R"doc(Copy of the Waveform between indices i_t_a and i_t_b, only ell=2 modes.
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
    necessary))doc"},
    {"Waveform::SliceOfTimeIndicesWithoutModes", R"doc(Copy of the Waveform between indices i_t_a and i_t_b without mode data.
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
    necessary))doc"},
    {"Waveform::SliceOfTimes", R"doc(Copy the Waveform between t_a and t_b.
======================================
  Parameters
  ----------
    const double t_a = -1e300
    const double t_b = 1e300

  Returns
  -------
    Waveform)doc"},
    {"Waveform::SliceOfTimesWithEll2", R"doc(Copy of the Waveform between t_a and t_b, only ell=2 modes.
===========================================================
  Parameters
  ----------
    const double t_a = -1e300
    const double t_b = 1e300

  Returns
  -------
    Waveform)doc"},
    {"Waveform::SliceOfTimesWithoutModes", R"doc(Copy of the Waveform between t_a and t_b without mode data.
===========================================================
  Parameters
  ----------
    const double t_a = -1e300
    const double t_b = 1e300

  Returns
  -------
    Waveform)doc"},
    {"Waveform::SpinWeight", R"doc(Parameters
----------
  (none)

Returns
-------
  int)doc"},
    {"Waveform::T", R"doc(Parameters
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
  const vector<double>&)doc"},
    {"Waveform::TransformToAngularVelocityFrame", R"doc(Transform Waveform to frame aligned with angular-velocity vector.
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
    or [2,3,4], for example, restricts the range of the sum.)doc"},
    {"Waveform::TransformToCoprecessingFrame", R"doc(Transform Waveform to co-precessing frame.
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
    or [2,3,4], for example, restricts the range of the sum.)doc"},
    {"Waveform::TransformToCorotatingFrame", R"doc(Transform Waveform to corotating frame.
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
    or [2,3,4], for example, restricts the range of the sum.)doc"},
    {"Waveform::TransformToInertialFrame", R"doc(Transform Waveform to an inertial frame.
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
    observations should be made.)doc"},
    {"Waveform::Translate", R"doc(Translate the waveform data by some series of spatial translations.
===================================================================
  Parameters
  ----------
    const vector<vector<double>>& deltax

  Returns
  -------
    Waveform)doc"},
    {"Waveform::Waveform", R"doc(Default constructor for an empty object.
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
    Arguments are T, LM, Data, which consist of the explicit data.)doc"},
    {"Waveform::XParityAntisymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::XParityInvolution", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::XParitySymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::XParityViolationNormalized", R"doc(Parameters
----------
  vector<int> Lmodes = vector<int>(0)

Returns
-------
  vector<double>)doc"},
    {"Waveform::XParityViolationSquared", R"doc(Parameters
----------
  vector<int> Lmodes = vector<int>(0)

Returns
-------
  vector<double>)doc"},
    {"Waveform::YParityAntisymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::YParityInvolution", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::YParitySymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::YParityViolationNormalized", R"doc(Parameters
----------
  vector<int> Lmodes = vector<int>(0)

Returns
-------
  vector<double>)doc"},
    {"Waveform::YParityViolationSquared", R"doc(Parameters
----------
  vector<int> Lmodes = vector<int>(0)

Returns
-------
  vector<double>)doc"},
    {"Waveform::ZParityAntisymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::ZParityInvolution", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::ZParitySymmetricPart", R"doc(Parameters
----------
  (none)

Returns
-------
  Waveform)doc"},
    {"Waveform::ZParityViolationNormalized", R"doc(Parameters
----------
  vector<int> Lmodes = vector<int>(0)

Returns
-------
  vector<double>)doc"},
    {"Waveform::ZParityViolationSquared", R"doc(Parameters
----------
  vector<int> Lmodes = vector<int>(0)

Returns
-------
  vector<double>)doc"},
    {"Waveform::operator()", R"doc(Parameters
----------
  const unsigned int Mode

Returns
-------
  const complex<double> *

Parameters
----------
  const unsigned int Mode
  const unsigned int TimeIndex

Returns
-------
  complex<double>)doc"},
    {"Waveform::operator*", R"doc(Parameters
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
  Waveform)doc"},
    {"Waveform::operator+", R"doc(Parameters
----------
  const Waveform& B

Returns
-------
  Waveform)doc"},
    {"Waveform::operator-", R"doc(Parameters
----------
  const Waveform& B

Returns
-------
  Waveform)doc"},
    {"Waveform::operator/", R"doc(Parameters
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
  Waveform)doc"},
    {"Waveform::operator=", R"doc(Assignment operator.
====================
  Parameters
  ----------
    const Waveform& 

  Returns
  -------
    Waveform&)doc"},
    {"Waveform::swap", R"doc(Efficiently swap data between two Waveform objects.
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
    pointers to data, for efficiency.)doc"},
    {"Waveform::~Waveform", R"doc(Parameters
----------
  (none)

Returns
-------
  ~Waveform)doc"},
    {"WaveformAtAPointFT", R"doc(class GWFrames::WaveformAtAPointFT
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
    bool mNormalized)doc"},
    {"WaveformAtAPointFT::F", R"doc(Returns the physical frequencies in Hz.
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
    const double&)doc"},
    {"WaveformAtAPointFT::Im", R"doc(Parameters
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
  const double&)doc"},
    {"WaveformAtAPointFT::InnerProduct", R"doc(Parameters
----------
  const WaveformAtAPointFT& B
  const vector<double>& InversePSD

Returns
-------
  double)doc"},
    {"WaveformAtAPointFT::InversePSD", R"doc(Parameters
----------
  const string& Detector = 'AdvLIGO_ZeroDet_HighP'

Returns
-------
  vector<double>)doc"},
    {"WaveformAtAPointFT::IsNormalized", R"doc(Parameters
----------
  (none)

Returns
-------
  bool)doc"},
    {"WaveformAtAPointFT::Match", R"doc(Compute the match between two WaveformAtAPointFT.
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
    void)doc"},
    {"WaveformAtAPointFT::NFreq", R"doc(Parameters
----------
  (none)

Returns
-------
  unsigned int)doc"},
    {"WaveformAtAPointFT::Normalize", R"doc(Parameters
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
  WaveformAtAPointFT&)doc"},
    {"WaveformAtAPointFT::Re", R"doc(Parameters
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
  const double&)doc"},
    {"WaveformAtAPointFT::SNR", R"doc(Parameters
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
  double)doc"},
    {"WaveformAtAPointFT::WaveformAtAPointFT", R"doc(Parameters
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
  WaveformAtAPointFT)doc"},
    {"WaveformAtAPointFT::ZeroAbove", R"doc(Parameters
----------
  const double Frequency

Returns
-------
  WaveformAtAPointFT&)doc"},
    {"WaveformUtilities", R"doc(namespace WaveformUtilities
===========================)doc"},
    {"WaveformUtilities::Interpolate", R"doc(Parameters
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
  vector<double>)doc"},
    {"WaveformUtilities::Interpolator", R"doc(class WaveformUtilities::Interpolator
=====================================
  Member variables
  ----------------
    int n
    int mm
    int jsav
    int cor
    int dj
    const vector<double>& xx
    const vector<double>& yy)doc"},
    {"WaveformUtilities::Interpolator::interp", R"doc(Parameters
----------
  double x

Returns
-------
  double)doc"},
    {"WaveformUtilities::Interpolator::rawinterp", R"doc(Parameters
----------
  int jlo
  double x

Returns
-------
  double)doc"},
    {"WaveformUtilities::InverseNoiseCurve", R"doc(Parameters
----------
  const vector<double>& F
  const string& Detector = 'AdvLIGO_ZeroDet_HighP'
  const double NoiseFloor = 0.0

Returns
-------
  vector<double>)doc"},
    {"WaveformUtilities::NoiseCurve", R"doc(Parameters
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
  AdvLIGO_NSNSOptimal  Collin Capano's fit for the NS-NS optimized noise curve)doc"},
    {"WaveformUtilities::PolynomialInterpolator", R"doc(class WaveformUtilities::PolynomialInterpolator
===============================================
  Member variables
  ----------------
    double dy)doc"},
    {"WaveformUtilities::PolynomialInterpolator::PolynomialInterpolator", R"doc(Parameters
----------
  const vector<double>& xv
  const vector<double>& yv
  int m

Returns
-------
  PolynomialInterpolator)doc"},
    {"WaveformUtilities::SplineCumulativeIntegral", R"doc(Parameters
----------
  const vector<double>& X1
  const vector<double>& Y1

Returns
-------
  double)doc"},
    {"WaveformUtilities::SplineIntegral", R"doc(Parameters
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
  vector<double>)doc"},
    {"WaveformUtilities::SplineIntegrator", R"doc(class WaveformUtilities::SplineIntegrator
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
    void SetUpIntegrationConstants)doc"},
    {"WaveformUtilities::SplineIntegrator::SplineIntegrator", R"doc(Parameters
----------
  const vector<double>& xv
  const vector<double>& yv
  double yp1 = 1.e99
  double ypn = 1.e99

Returns
-------
  SplineIntegrator)doc"},
    {"WaveformUtilities::SplineInterpolator", R"doc(class WaveformUtilities::SplineInterpolator
===========================================
  Member variables
  ----------------
    vector<double> y2)doc"},
    {"WaveformUtilities::SplineInterpolator::SplineInterpolator", R"doc(Parameters
----------
  const vector<double>& xv
  const vector<double>& yv
  double yp1 = 1.e99
  double ypn = 1.e99

Returns
-------
  SplineInterpolator)doc"},
    {"WaveformUtilities::TimeToFrequency", R"doc(This function creates a frequency vector in (0 -> positives -> negatives -> 0) order.
=====================================================================================
  Parameters
  ----------
    const vector<double>& Time

  Returns
  -------
    vector<double>)doc"},
    {"WaveformUtilities::TimeToPositiveFrequencies", R"doc(This function returns the positive half of the frequencies, so returned size is 1/2 input size + 1.
===================================================================================================
  Parameters
  ----------
    const vector<double>& Time

  Returns
  -------
    vector<double>)doc"},
    {"WaveformUtilities::WrapVecDoub", R"doc(class WaveformUtilities::WrapVecDoub
====================================
  This helper class is useful for untangling NR's fft routines.

  Member variables
  ----------------
    vector<double> vvec
    vector<double>& v
    int n
    int mask)doc"},
    {"WaveformUtilities::WrapVecDoub::WrapVecDoub", R"doc(Parameters
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
  WrapVecDoub)doc"},
    {"WaveformUtilities::WrapVecDoub::imag", R"doc(Parameters
----------
  int i

Returns
-------
  double&)doc"},
    {"WaveformUtilities::WrapVecDoub::operator[]", R"doc(Parameters
----------
  int i

Returns
-------
  complex<double>&)doc"},
    {"WaveformUtilities::WrapVecDoub::real", R"doc(Parameters
----------
  int i

Returns
-------
  double&)doc"},
    {"WaveformUtilities::WrapVecDoub::validate", R"doc(Parameters
----------
  (none)

Returns
-------
  void)doc"},
    {"WaveformUtilities::convlv", R"doc(Parameters
----------
  const vector<double>& data
  const vector<double>& respns
  const int isign

Returns
-------
  vector<double>)doc"},
    {"WaveformUtilities::dft", R"doc(Parameters
----------
  vector<double>& data

Returns
-------
  void

Description
-----------
  The following call the Numerical Recipes fft routines (from fourier.h) Note
  that the returned quantities represent the bare fft sum, with no
  normalization constants)doc"},
    {"WaveformUtilities::idft", R"doc(Parameters
----------
  vector<double>& data

Returns
-------
  void)doc"},
    {"WaveformUtilities::realdft", R"doc(Parameters
----------
  vector<double>& data

Returns
-------
  void)doc"},
    {"WignerCoefficientFunctor", R"doc(class GWFrames::WignerCoefficientFunctor
========================================
  Storage for pre-calculated values of coefficients used in calculating Wigner
  D matrices.

  The object is a functor for certain coefficients needed when computing the
  Wigner D matrices. E.g., WignerCoefficientFunctor WignerCoeff;
const double
  Coeff_4_3_m2 = WignerCoeff(4,3,-2);

  Note that values up to WignerCoeff(ellMax_GWFrames, ...) are calculated,
  where ellMax_GWFrames is set in Utilities.hpp. Requests for greater values
  are answered with an exception.)doc"},
    {"WignerD", R"doc(class GWFrames::WignerD
=======================
  Calculate values of Wigner's D matrices.

  This object is a functor for calculating the values of Wigner's D matrices
  for many values of (ell,m',m), given some rotation. E.g., WignerD D(R);
const
  double D_5_2_m4 = D(5,2,-4);
const double D_5_2_m3 = D(5,2,-3);

  The idea is that the rotation can be set to some rotor R, appropriate for
  some instant in time (possibly using the SetRotation member function), which
  caches certain values that need to be re-computed for any R. Work at that
  instant of time will usually require computing many values, making this an
  efficient method of doing so.)doc"},
    {"abs", R"doc(Parameters
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
  vector<double>)doc"},
    {"angle", R"doc(Parameters
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
  vector<double>)doc"},
    {"commutator", R"doc(Parameters
----------
  const Quaternion& Q
  const Quaternion& P

Returns
-------
  Quaternion)doc"},
    {"complexi", R"doc(Parameters
----------
  0. 0
  1. 0

Returns
-------
  const complex<double>)doc"},
    {"conjugate", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"cross", R"doc(Parameters
----------
  const Quaternion& Q
  const Quaternion& P

Returns
-------
  Quaternion)doc"},
    {"dot", R"doc(Parameters
----------
  const Quaternion& Q
  const Quaternion& P

Returns
-------
  double)doc"},
    {"dotproduct", R"doc(Parameters
----------
  const double * a
  const double * b

Returns
-------
  double)doc"},
    {"ellmax", R"doc(Parameters
----------
  const int maxindex

Returns
-------
  int)doc"},
    {"exp", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"hunt", R"doc(Find the largest index i with x[i]<=x_i.
========================================
  Parameters
  ----------
    const vector<double>& x
    const double x_i
    unsigned int i = 0

  Returns
  -------
    unsigned int

  Description
  -----------
    Based on the Numerical Recipes routine of the same name)doc"},
    {"inverse", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"log", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"logRotor", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"myGSLErrorHandler", R"doc(Parameters
----------
  const char * reason
  const char * file
  int line
  int gsl_errno

Returns
-------
  void)doc"},
    {"normalized", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"normsquared", R"doc(Parameters
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
  vector<double>)doc"},
    {"operator*", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"operator+", R"doc(Parameters
----------
  const double& a
  const DataGrid& b

Returns
-------
  DataGrid

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
  vector<double>)doc"},
    {"operator-", R"doc(Parameters
----------
  const double& a
  const DataGrid& b

Returns
-------
  DataGrid

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
  vector<Quaternion>)doc"},
    {"operator/", R"doc(Parameters
----------
  const double& a
  const DataGrid& b

Returns
-------
  DataGrid

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
  vector<Quaternion>)doc"},
    {"operator<<", R"doc(Print the quaternion nicely to stream.
======================================
  Parameters
  ----------
    ostream& out
    const Quaternions::Quaternion& q

  Returns
  -------
    ostream&)doc"},
    {"pow", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"sqrt", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"sqrtOfRotor", R"doc(Parameters
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
  vector<Quaternion>)doc"},
    {"std", R"doc(namespace std
=============
  STL namespace.)doc"},
    {"tolower", R"doc(Parameters
----------
  const string& A

Returns
-------
  string)doc"},
    {"vFromOneOverK", R"doc(Derive three-velocity from the inverse conformal metric.
========================================================
  Parameters
  ----------
    const Modes& OneOverK

  Returns
  -------
    ThreeVector)doc"},
    {"vec", R"doc(Parameters
----------
  const Quaternion& Q

Returns
-------
  vector<double>




Parameters
----------
  const vector<Quaternions::Quaternion>& Q

Returns
-------
  vector<vector<double>>)doc"},
    {"virtualWaveformUtilities::Interpolator::~Interpolator", R"doc(Parameters
----------
  (none)

Returns
-------
  ~Interpolator)doc"},
    {"zHat", R"doc(Parameters
----------
  0 
  0 
  0 
  1 

Returns
-------
  const Quaternion)doc"},
    {"zero", R"doc(Parameters
----------
  0. 0
  0. 0

Returns
-------
  const complex<double>)doc"},
  };
  return docstrings;
}

}  // namespace GWFramesBindings

#endif  // GWFRAMES_BINDINGS_DOCSTRINGS_HPP
