// Docstrings for the Python bindings, in the NumPy style that Sphinx renders
// on the documentation website and that `help()` shows in a terminal.  Each
// key is the C++ name of the bound function, without the leading `GWFrames::`
// or `Quaternions::` namespace, or the Python class name for a class.  The
// `Doc` function in `bindings.hpp` looks them up; for overloaded functions,
// the docstring belongs on the last overload bound.

#ifndef GWFRAMES_BINDINGS_DOCSTRINGS_HPP
#define GWFRAMES_BINDINGS_DOCSTRINGS_HPP

#include <string>
#include <unordered_map>

namespace GWFramesBindings {

inline const std::unordered_map<std::string, const char*>& Docstrings() {
  static const std::unordered_map<std::string, const char*> docstrings = {
    {"Quaternion", R"doc(A quaternion, used mainly to represent rotations.

A quaternion is stored as its four components ``[w, x, y, z]``, where ``w`` is
the scalar part and ``(x, y, z)`` is the vector part.  The components are
read and written by index, from ``Q[0]`` (the scalar part) to ``Q[3]``, and
negative indices count from the end, as for a list; ``len(Q)`` is 4.
Quaternions may be copied, deep-copied, and pickled.

The arithmetic operators ``+``, ``-``, ``*``, and ``/`` act between two
quaternions, and between a quaternion and a float, which is treated as a
quaternion with only a scalar part.  The product ``*`` is the quaternion
product, and ``P / Q`` is ``P * Q.inverse()``.  The power ``Q ** t`` is
`Quaternions.Quaternion.pow`.  The comparisons ``==`` and ``!=`` compare all
four components exactly.

A unit quaternion (a rotor) ``R`` represents the rotation that takes a vector
``v``, written as a quaternion with zero scalar part, to
``R * v * R.inverse()``.  A sequence of three or four floats is accepted
wherever a `Quaternion` is expected, and is converted as in the constructor
``Quaternion(q)`` below.  The module also defines the constants
``Quaternions.Zero``, ``Quaternions.One``, ``Quaternions.xHat``,
``Quaternions.yHat``, and ``Quaternions.zHat``.

A `Quaternion` may be constructed in any of these ways:

``Quaternion()``
    The zero quaternion.
``Quaternion(Q)``
    A copy of the quaternion ``Q``.
``Quaternion(vartheta, varphi)``
    The rotor :math:`e^{\varphi \hat{z}/2}\, e^{\vartheta \hat{y}/2}`, which
    rotates the z axis onto the point with polar angle ``vartheta`` and
    azimuthal angle ``varphi``.
``Quaternion(alpha, beta, gamma)``
    The rotor :math:`e^{\alpha \hat{z}/2}\, e^{\beta \hat{y}/2}\,
    e^{\gamma \hat{z}/2}` for the Euler angles ``alpha``, ``beta``, and
    ``gamma`` in the z-y-z convention: a rotation by ``gamma`` about the z
    axis, followed by a rotation by ``beta`` about the fixed y axis, and
    finally a rotation by ``alpha`` about the fixed z axis.
``Quaternion(w, x, y, z)``
    The quaternion with the given components.
``Quaternion(q)``
    The quaternion with components taken from the sequence ``q`` of floats.
    Four elements are the components ``[w, x, y, z]``; three elements are the
    vector part ``[x, y, z]`` of a quaternion with zero scalar part.  Any
    other length raises an error.
``Quaternion(angle, axis)``
    The rotor :math:`\cos(\theta/2) + \sin(\theta/2)\, \hat{n}` for a rotation
    through ``angle`` :math:`\theta` about the three-vector ``axis``
    :math:`\hat{n}`, which is assumed to be normalized.)doc"},
    {"Quaternion::pow", R"doc(Return this quaternion raised to a real or quaternion power.

The power is computed as ``exp(log(self) * t)`` or ``exp(log(self) * Q)``,
with the logarithm of `Quaternions.Quaternion.log`; it therefore raises an
error when this quaternion is a negative real number.  The same operation is
written ``self ** t`` or ``self ** Q``.  `Quaternions.pow` applies it to
arrays.

Parameters
----------
t : float
    The real exponent, in the form ``pow(t)``.
Q : Quaternion
    The quaternion exponent, in the form ``pow(Q)``.  It multiplies the
    logarithm from the right.

Returns
-------
Quaternion
    The power.)doc"},
    {"Quaternion::cross", R"doc(Return the cross product of the vector parts of two quaternions.

The scalar parts are ignored, and the result is a quaternion with zero scalar
part whose vector part is the cross product of the vector part of this
quaternion with that of ``Q``.

Parameters
----------
Q : Quaternion
    The second factor.

Returns
-------
Quaternion
    The cross product, as a quaternion with zero scalar part.)doc"},
    {"Quaternion::commutator", R"doc(Return the commutator of this quaternion with another.

The commutator is ``self * Q - Q * self``, which equals twice
`Quaternions.Quaternion.cross` and has zero scalar part.

Parameters
----------
Q : Quaternion
    The second quaternion.

Returns
-------
Quaternion
    The commutator.)doc"},
    {"Quaternion::dot", R"doc(Return the dot product of the vector parts of two quaternions.

The scalar parts are ignored, so this is
``self[1]*Q[1] + self[2]*Q[2] + self[3]*Q[3]``.

Parameters
----------
Q : Quaternion
    The second factor.

Returns
-------
float
    The dot product of the vector parts.)doc"},
    {"Quaternion::abs", R"doc(Return the magnitude of this quaternion.

The magnitude is the square root of the sum of the squares of all four
components; see also `Quaternions.Quaternion.normsquared`.

Returns
-------
float
    The magnitude.)doc"},
    {"Quaternion::log", R"doc(Return the natural logarithm of this quaternion.

For a quaternion with magnitude :math:`|Q|`, vector part :math:`\vec{v}`, and
scalar part :math:`w`, the logarithm has scalar part :math:`\ln |Q|` and
vector part :math:`\operatorname{atan2}(|\vec{v}|, w)\, \vec{v}/|\vec{v}|`.
For a rotor :math:`e^{\theta \hat{n}/2}` with :math:`0 \le \theta \le 2\pi`,
it is therefore :math:`\theta \hat{n}/2`.

The logarithm of a negative real number has infinitely many values, so this
function raises an error for one; `Quaternions.Quaternion.logRotor` returns a
value instead.

Returns
-------
Quaternion
    The logarithm.)doc"},
    {"Quaternion::logRotor", R"doc(Return the natural logarithm, choosing a value for negative real numbers.

This is the same as `Quaternions.Quaternion.log`, except for negative real
numbers, for which the logarithm has infinitely many values.  Instead of
raising an error, this function prints a warning and returns
:math:`\ln|w| + \pi \hat{x}`; for the rotor -1, which represents the identity
rotation, that is :math:`\pi \hat{x}`.

Returns
-------
Quaternion
    The logarithm.)doc"},
    {"Quaternion::exp", R"doc(Return the exponential of this quaternion.

For scalar part :math:`w` and vector part :math:`\vec{v}`, the exponential is
:math:`e^{w} \left[\cos|\vec{v}| + \sin|\vec{v}|\,
\vec{v}/|\vec{v}|\right]`.  In particular, the exponential of
:math:`\theta \hat{n}/2` is the rotor for a rotation through :math:`\theta`
about :math:`\hat{n}`.

Returns
-------
Quaternion
    The exponential.)doc"},
    {"Quaternion::sqrt", R"doc(Return the principal square root of this quaternion.

The result is ``(self/abs(self) + 1) * sqrt(abs(self) / (2 + 2*w/abs(self)))``,
where ``w`` is the scalar part.  It is undefined (NaN) for a negative real
number.  For a unit quaternion, `Quaternions.Quaternion.sqrtOfRotor` gives
the same result more cheaply, and handles the rotor -1.

Returns
-------
Quaternion
    The square root.)doc"},
    {"Quaternion::angle", R"doc(Return the full angle of the rotation that this rotor represents.

The result is ``2 * abs(self.logRotor())``.  For a unit quaternion
:math:`e^{\theta \hat{n}/2}` this is the full rotation angle :math:`\theta`,
not half of it, and lies between 0 and :math:`2\pi`; a rotor and its negative
give angles that add to :math:`2\pi`.  For a quaternion that is not
normalized, the logarithm of its magnitude also contributes.

Returns
-------
float
    The rotation angle.)doc"},
    {"Quaternion::inverse", R"doc(Return the multiplicative inverse of this quaternion.

The inverse is the conjugate divided by `Quaternions.Quaternion.normsquared`;
for a unit quaternion, it equals the conjugate.

Returns
-------
Quaternion
    The inverse.)doc"},
    {"Quaternion::conjugate", R"doc(Return the conjugate of this quaternion.

The conjugate has the same scalar part and the negated vector part.

Returns
-------
Quaternion
    The conjugate.)doc"},
    {"Quaternion::normalized", R"doc(Return this quaternion divided by its magnitude.

Returns
-------
Quaternion
    The unit quaternion in the same direction.)doc"},
    {"Quaternion::normsquared", R"doc(Return the squared magnitude of this quaternion.

This is the sum of the squares of all four components.

Returns
-------
float
    The squared magnitude.)doc"},
    {"Quaternion::sqrtOfRotor", R"doc(Return the square root of this quaternion, assuming that it is a rotor.

The result is ``(1 + self) / sqrt(2 + 2*w)``, where ``w`` is the scalar part,
which is the principal square root only when this quaternion has unit
magnitude.  If ``w`` is within 1e-14 of -1, the result is ``xHat``, one of
the infinitely many square roots of -1.

Returns
-------
Quaternion
    The square root.)doc"},
    {"Quaternion::vec", R"doc(Return the vector part of this quaternion.

Returns
-------
numpy.ndarray of float, shape (3,)
    The components ``[x, y, z]``.)doc"},
    {"Quaternion::XParityConjugateSpinor", R"doc(Return the image of this rotor under the reflection that reverses the x axis.

Reflecting the space through the y-z plane takes the rotation through an angle
about the axis :math:`(n_x, n_y, n_z)` to the rotation through the same angle
about :math:`(n_x, -n_y, -n_z)`, because a rotation axis is a pseudovector.
The result is therefore ``Quaternion(w, x, -y, -z)``.

Returns
-------
Quaternion
    The reflected rotor.)doc"},
    {"Quaternion::YParityConjugateSpinor", R"doc(Return the image of this rotor under the reflection that reverses the y axis.

Reflecting the space through the x-z plane takes the rotation through an angle
about the axis :math:`(n_x, n_y, n_z)` to the rotation through the same angle
about :math:`(-n_x, n_y, -n_z)`, because a rotation axis is a pseudovector.
The result is therefore ``Quaternion(w, -x, y, -z)``.

Returns
-------
Quaternion
    The reflected rotor.)doc"},
    {"Quaternion::ZParityConjugateSpinor", R"doc(Return the image of this rotor under the reflection that reverses the z axis.

Reflecting the space through the x-y plane takes the rotation through an angle
about the axis :math:`(n_x, n_y, n_z)` to the rotation through the same angle
about :math:`(-n_x, -n_y, n_z)`, because a rotation axis is a pseudovector.
The result is therefore ``Quaternion(w, -x, -y, z)``.

Returns
-------
Quaternion
    The reflected rotor.)doc"},
    {"Quaternion::ParityConjugateSpinor", R"doc(Return the image of this rotor under the parity transformation.

The parity transformation reverses all three axes.  A rotation axis is a
pseudovector, which parity leaves unchanged, so the result is a copy of this
quaternion.

Returns
-------
Quaternion
    A copy of this rotor.)doc"},
    {"Quaternion::IntrinsicDistance", R"doc(Return the intrinsic (geodesic) distance between this rotor and another.

The distance is ``2 * abs((self * Q.inverse()).logRotor())``, which for unit
quaternions is the angle of the rotation that takes ``Q`` to this rotor, as
given by `Quaternions.Quaternion.angle`.  It lies between 0 and
:math:`2\pi`.

Parameters
----------
Q : Quaternion
    The other rotor.

Returns
-------
float
    The intrinsic distance.)doc"},
    {"Quaternion::ChordalDistance", R"doc(Return the chordal distance between this quaternion and another.

The distance is the magnitude of the difference, ``abs(self - Q)``.

Parameters
----------
Q : Quaternion
    The other quaternion.

Returns
-------
float
    The chordal distance.)doc"},
    {"Quaternion::str", R"doc(Return the four components separated by spaces.

The components appear in the order ``w x y z``, with six significant digits.

Returns
-------
str
    The components.)doc"},
    {"pow", R"doc(Raise quaternions to real or quaternion powers.

Each power is computed as ``exp(log(Q) * x)`` or ``exp(log(Q) * P)``, as in
`Quaternions.Quaternion.pow`, so it raises an error when a base is a negative
real number.  Either argument may be a single value or an array.  A single
value is combined with every element of an array, and two arrays, which must
have the same length, are combined element by element.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The base or bases.
x : float or array_like of float
    The real exponent or exponents.  Pass either ``x`` or ``P``.
P : Quaternion or array_like of Quaternion
    The quaternion exponent or exponents, which multiply the logarithm from
    the right.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The power, or an array of powers if either argument is an array.)doc"},
    {"cross", R"doc(Return the cross product of the vector parts of two quaternions.

The scalar parts are ignored, and the result has zero scalar part; see
`Quaternions.Quaternion.cross`.

Parameters
----------
Q : Quaternion
    The first factor.
P : Quaternion
    The second factor.

Returns
-------
Quaternion
    The cross product, as a quaternion with zero scalar part.)doc"},
    {"commutator", R"doc(Return the commutator of two quaternions.

The commutator is ``Q * P - P * Q``, which equals twice `Quaternions.cross`.

Parameters
----------
Q : Quaternion
    The first quaternion.
P : Quaternion
    The second quaternion.

Returns
-------
Quaternion
    The commutator.)doc"},
    {"dot", R"doc(Return the dot product of the vector parts of two quaternions.

The scalar parts are ignored; see `Quaternions.Quaternion.dot`.

Parameters
----------
Q : Quaternion
    The first factor.
P : Quaternion
    The second factor.

Returns
-------
float
    The dot product of the vector parts.)doc"},
    {"abs", R"doc(Return the magnitude of a quaternion or of each quaternion in an array.

See `Quaternions.Quaternion.abs`.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
float or numpy.ndarray of float
    The magnitude or magnitudes.)doc"},
    {"log", R"doc(Return the natural logarithm of a quaternion or of each in an array.

See `Quaternions.Quaternion.log`; an error is raised for a negative real
number, for which `Quaternions.logRotor` returns a value instead.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The logarithm or logarithms.)doc"},
    {"logRotor", R"doc(Return the logarithm of quaternions, choosing a value for negative reals.

This is the same as `Quaternions.log`, except that a negative real number
gives :math:`\ln|w| + \pi \hat{x}`, with a warning, instead of an error; see
`Quaternions.Quaternion.logRotor`.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The logarithm or logarithms.)doc"},
    {"exp", R"doc(Return the exponential of a quaternion or of each quaternion in an array.

See `Quaternions.Quaternion.exp`.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The exponential or exponentials.)doc"},
    {"sqrt", R"doc(Return the principal square root of a quaternion or of each in an array.

See `Quaternions.Quaternion.sqrt`.  For unit quaternions,
`Quaternions.sqrtOfRotor` gives the same result more cheaply.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The square root or square roots.)doc"},
    {"angle", R"doc(Return the full rotation angle of a rotor or of each rotor in an array.

The angle is ``2 * abs(logRotor(Q))``, which for a unit quaternion
:math:`e^{\theta \hat{n}/2}` is the full rotation angle :math:`\theta`, not
half of it, between 0 and :math:`2\pi`; see `Quaternions.Quaternion.angle`.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The rotor or rotors.

Returns
-------
float or numpy.ndarray of float
    The rotation angle or angles.)doc"},
    {"inverse", R"doc(Return the inverse of a quaternion or of each quaternion in an array.

See `Quaternions.Quaternion.inverse`.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The inverse or inverses.)doc"},
    {"conjugate", R"doc(Return the conjugate of a quaternion or of each quaternion in an array.

The conjugate has the same scalar part and the negated vector part.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The conjugate or conjugates.)doc"},
    {"normalized", R"doc(Return a quaternion, or each quaternion in an array, divided by its magnitude.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The unit quaternion or quaternions.)doc"},
    {"normsquared", R"doc(Return the squared magnitude of a quaternion or of each in an array.

The squared magnitude is the sum of the squares of all four components.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
float or numpy.ndarray of float
    The squared magnitude or magnitudes.)doc"},
    {"sqrtOfRotor", R"doc(Return the square root of a rotor or of each rotor in an array.

The result is the principal square root only for unit quaternions; see
`Quaternions.Quaternion.sqrtOfRotor`.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The rotor or rotors.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The square root or square roots.)doc"},
    {"vec", R"doc(Return the vector part of a quaternion or of each quaternion in an array.

Parameters
----------
Q : Quaternion or array_like of Quaternion
    The quaternion or quaternions.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The components ``[x, y, z]``, with one row per quaternion if ``Q`` is an
    array of ``N`` quaternions.)doc"},
    {"XParityConjugateSpinor", R"doc(Return the image of a rotor under the reflection that reverses the x axis.

The result is ``Quaternion(w, x, -y, -z)``; see
`Quaternions.Quaternion.XParityConjugateSpinor`.

Parameters
----------
Q : Quaternion
    The rotor.

Returns
-------
Quaternion
    The reflected rotor.)doc"},
    {"YParityConjugateSpinor", R"doc(Return the image of a rotor under the reflection that reverses the y axis.

The result is ``Quaternion(w, -x, y, -z)``; see
`Quaternions.Quaternion.YParityConjugateSpinor`.

Parameters
----------
Q : Quaternion
    The rotor.

Returns
-------
Quaternion
    The reflected rotor.)doc"},
    {"ZParityConjugateSpinor", R"doc(Return the image of a rotor under the reflection that reverses the z axis.

The result is ``Quaternion(w, -x, -y, z)``; see
`Quaternions.Quaternion.ZParityConjugateSpinor`.

Parameters
----------
Q : Quaternion
    The rotor.

Returns
-------
Quaternion
    The reflected rotor.)doc"},
    {"ParityConjugateSpinor", R"doc(Return the image of a rotor under the parity transformation.

Parity leaves a rotation axis unchanged, so the result is a copy of ``Q``; see
`Quaternions.Quaternion.ParityConjugateSpinor`.

Parameters
----------
Q : Quaternion
    The rotor.

Returns
-------
Quaternion
    A copy of the rotor.)doc"},
    {"IntrinsicDistance", R"doc(Return the intrinsic (geodesic) distance between two rotors.

The distance is ``2 * abs(logRotor(R_A * R_B.inverse()))``, which for unit
quaternions is the angle of the rotation that takes ``R_B`` to ``R_A``.  It
lies between 0 and :math:`2\pi`.

Parameters
----------
R_A : Quaternion
    The first rotor.
R_B : Quaternion
    The second rotor.

Returns
-------
float
    The intrinsic distance.)doc"},
    {"ChordalDistance", R"doc(Return the chordal distance between two quaternions.

The distance is the magnitude of the difference, ``abs(R_A - R_B)``.

Parameters
----------
R_A : Quaternion
    The first quaternion.
R_B : Quaternion
    The second quaternion.

Returns
-------
float
    The chordal distance.)doc"},
    {"Slerp", R"doc(Interpolate spherically and linearly between two rotors.

The result is ``pow(Qb / Qa, tau) * Qa``, which moves from ``Qa`` at
``tau = 0`` to ``Qb`` at ``tau = 1`` along a geodesic.  Values of ``tau``
outside that interval extrapolate along the same geodesic.

Parameters
----------
tau : float
    The interpolation parameter.
Qa : Quaternion
    The rotor at ``tau = 0``.
Qb : Quaternion
    The rotor at ``tau = 1``.

Returns
-------
Quaternion
    The interpolated rotor.)doc"},
    {"QuaternionArray", R"doc(Construct an array of quaternions from arrays of numbers.

The forms of this function match the constructors of `Quaternions.Quaternion`,
applied element by element:

``QuaternionArray(vartheta, varphi)``
    Rotors that rotate the z axis onto the points with polar angles
    ``vartheta`` and azimuthal angles ``varphi``.
``QuaternionArray(alpha, beta, gamma)``
    Rotors for the z-y-z Euler angles ``alpha``, ``beta``, and ``gamma``.
``QuaternionArray(w, x, y, z)``
    Quaternions with the given components.
``QuaternionArray(q)``
    Quaternions whose components are the rows of ``q``; each row has three
    elements (the vector part) or four (all components).
``QuaternionArray(angle, axis)``
    Rotors for rotations through each ``angle`` about the corresponding row of
    ``axis``, which is assumed to be normalized.

The length of the result is that of the first argument, and the other
arguments must be at least as long.

Parameters
----------
vartheta, varphi : array_like of float
    Polar and azimuthal angles.
alpha, beta, gamma : array_like of float
    Euler angles.
w, x, y, z : array_like of float
    Components.
q : array_like of float, shape (N, 3) or (N, 4)
    Components, one quaternion per row.
angle : array_like of float
    Rotation angles.
axis : array_like of float, shape (N, 3)
    Unit rotation axes, one per row.

Returns
-------
numpy.ndarray of Quaternion
    The quaternions.)doc"},
    {"Component", R"doc(Return one component of each quaternion in an array.

Parameters
----------
Q : array_like of Quaternion
    The quaternions.
i : int
    Index of the component: 0 for the scalar part, and 1, 2, or 3 for the x,
    y, or z component of the vector part.  Any other value raises an error.

Returns
-------
numpy.ndarray of float
    The ``i``-th component of each quaternion.)doc"},
    {"Component0", R"doc(Return the scalar component of each quaternion in an array.

This is ``Component(Q, 0)``.

Parameters
----------
Q : array_like of Quaternion
    The quaternions.

Returns
-------
numpy.ndarray of float
    The scalar component of each quaternion.)doc"},
    {"Component1", R"doc(Return the x component of each quaternion in an array.

This is ``Component(Q, 1)``.

Parameters
----------
Q : array_like of Quaternion
    The quaternions.

Returns
-------
numpy.ndarray of float
    The x component of each quaternion.)doc"},
    {"Component2", R"doc(Return the y component of each quaternion in an array.

This is ``Component(Q, 2)``.

Parameters
----------
Q : array_like of Quaternion
    The quaternions.

Returns
-------
numpy.ndarray of float
    The y component of each quaternion.)doc"},
    {"Component3", R"doc(Return the z component of each quaternion in an array.

This is ``Component(Q, 3)``.

Parameters
----------
Q : array_like of Quaternion
    The quaternions.

Returns
-------
numpy.ndarray of float
    The z component of each quaternion.)doc"},
    {"IndefiniteIntegral", R"doc(Integrate a quaternion function of time, returning the integral at each time.

Each component is integrated separately by the trapezoidal rule, starting from
zero at the first time.

Parameters
----------
RIn : array_like of Quaternion
    Values of the function.
tIn : array_like of float
    The corresponding times, of the same length as ``RIn``.

Returns
-------
numpy.ndarray of Quaternion
    The integral from ``tIn[0]`` to each time in ``tIn``.)doc"},
    {"DefiniteIntegral", R"doc(Integrate a quaternion function of time over all of the given times.

Each component is integrated separately by the trapezoidal rule.  The result
equals the last element of `Quaternions.IndefiniteIntegral`.

Parameters
----------
RIn : array_like of Quaternion
    Values of the function.
tIn : array_like of float
    The corresponding times, of the same length as ``RIn``.

Returns
-------
Quaternion
    The integral from ``tIn[0]`` to ``tIn[-1]``.)doc"},
    {"DifferentiateRotorByLogarithm", R"doc(Differentiate a rotor function of time by way of its logarithm.

The logarithm of each rotor is computed with `Quaternions.logRotor` and
differentiated with `Quaternions.QuaternionDerivative`, and the derivative of
the rotor is then found from the analytic formula for the derivative of the
exponential.  This is more complicated than differentiating the rotors
directly with `Quaternions.QuaternionDerivative`, but it may be more accurate
when the logarithm is smooth, and it serves as a test of the logarithm
formula.

Parameters
----------
RIn : array_like of Quaternion
    The rotors.
tIn : array_like of float
    The corresponding times.

Returns
-------
numpy.ndarray of Quaternion
    The time derivative of the rotors.)doc"},
    {"QuaternionDerivative", R"doc(Differentiate a quaternion function of time by finite differences.

The derivative is computed with the five-point formula of Eq. (A 5b) of
Bowen and Smith, "Derivative formulas and errors for non-uniformly spaced
points", which is fourth-order accurate, as they explain in Eqs. (B 9b) and
(B 10b).  The times need not be uniformly spaced.  With three or four points,
the derivative uses three-point formulas, at interior points from Sundqvist
and Veronis, Tellus XXII (1970), 1; with two points, it is the slope between
them.  Fewer than two points, or arrays of different lengths, raise an error.

Parameters
----------
f : array_like of Quaternion
    Values of the function.
t : array_like of float
    The corresponding times.

Returns
-------
numpy.ndarray of Quaternion
    The time derivative at each time.)doc"},
    {"MinimalRotation", R"doc(Return the minimal-rotation version of a frame.

The returned rotors take the z axis to the same direction as the rotors
``R`` do, but the rotation about that direction is adjusted to impose the
minimal-rotation condition: the angular velocity of the frame has no
component along its rotated z axis.  The adjustment is computed by
integrating the unwanted rotation rate over time and removing it, and the
process is repeated ``NIterations`` times to refine the result.

Parameters
----------
R : array_like of Quaternion
    The rotors of the input frame.
T : array_like of float
    The corresponding times.
NIterations : int, default 5
    Number of refinements.

Returns
-------
numpy.ndarray of Quaternion
    The rotors of the minimal-rotation frame.)doc"},
    {"PrescribedRotation", R"doc(Return a frame with the same z axis as the input and a prescribed spin rate.

The returned rotors take the z axis to the same direction as the rotors
``R`` do, but the rotation about that direction is adjusted, as in
`Quaternions.MinimalRotation`, so that the frame rotates about its rotated z
axis at the rate given by ``RotationRateAboutZ``.  The process is repeated
``NIterations`` times to refine the result, and progress messages are
printed to standard output.

Parameters
----------
RotationRateAboutZ : array_like of float
    Rate of rotation of the new frame about its rotated z axis at each time.
R : array_like of Quaternion
    The rotors of the input frame.
T : array_like of float
    The corresponding times.
NIterations : int, default 5
    Number of refinements.

Returns
-------
numpy.ndarray of Quaternion
    The rotors of the new frame.)doc"},
    {"FrameFromXY", R"doc(Construct the frame with the given X and Y basis vectors.

At each time, the returned rotor rotates the fixed basis vectors
:math:`(\hat{x}, \hat{y}, \hat{z})` onto the frame's basis vectors
:math:`(\hat{X}, \hat{Y}, \hat{Z})`.

Parameters
----------
X : array_like of Quaternion
    The X basis vector at each time, as a unit quaternion with zero scalar
    part.
Y : array_like of Quaternion
    The Y basis vector at each time, in the same form.

Returns
-------
numpy.ndarray of Quaternion
    The rotors of the frame.)doc"},
    {"FrameFromZ", R"doc(Construct the minimal-rotation frame with the given Z basis vector.

At each time, the returned rotor rotates the fixed basis vectors
:math:`(\hat{x}, \hat{y}, \hat{z})` onto the frame's basis vectors
:math:`(\hat{X}, \hat{Y}, \hat{Z})`.  The X and Y basis vectors are determined
by the minimal-rotation condition of `Quaternions.MinimalRotation`, which
leaves the initial rotation about the Z axis arbitrary.

Parameters
----------
Z : array_like of Quaternion
    The Z basis vector at each time, as a unit quaternion with zero scalar
    part.
T : array_like of float
    The corresponding times.
NIterations : int, default 5
    Number of refinements in `Quaternions.MinimalRotation`.

Returns
-------
numpy.ndarray of Quaternion
    The rotors of the frame.)doc"},
    {"FrameFromPrescribedRotation", R"doc(Construct the frame whose Z axis and spin rate follow an angular velocity.

At each time, the returned rotor rotates the fixed basis vectors
:math:`(\hat{x}, \hat{y}, \hat{z})` onto the frame's basis vectors
:math:`(\hat{X}, \hat{Y}, \hat{Z})`.  The Z basis vector is parallel to
``omega``, and the X and Y basis vectors are determined, with
`Quaternions.PrescribedRotation`, by the condition that the frame rotates
about its Z axis at the rate ``abs(omega)``.  This leaves the initial
rotation in the X-Y plane arbitrary.  Progress messages are printed to
standard output.

Parameters
----------
omega : array_like of Quaternion
    The angular-velocity vector at each time, as a quaternion with zero
    scalar part.
T : array_like of float
    The corresponding times.
NIterations : int, default 5
    Number of refinements in `Quaternions.PrescribedRotation`.

Returns
-------
numpy.ndarray of Quaternion
    The rotors of the frame.)doc"},
    {"FrameFromAngularVelocity", R"doc(Find the frame with the given angular velocity.

The frame starts as the identity at ``T[0]``.  The logarithm of its rotor is
integrated in time, using `Quaternions.FrameFromAngularVelocity_Integrand`
with ``Omega`` interpolated by cubic splines, and the result is exponentiated
and passed through `Quaternions.UnflipRotors`.  The logarithm is replaced by
an equivalent one of smaller magnitude whenever its magnitude exceeds
:math:`\pi/2`, much as a complex logarithm changes branch.

The accuracy is limited, because the integration cannot adapt its steps to
anything finer than the interpolated data.

Parameters
----------
Omega : array_like of Quaternion
    The angular-velocity vector at each time, as a quaternion with zero
    scalar part.
T : array_like of float
    The corresponding times.

Returns
-------
numpy.ndarray of Quaternion
    The rotors of the frame at the times ``T``.)doc"},
    {"FrameFromAngularVelocity_2D", R"doc(Find a frame from an angular velocity, keeping its logarithm in the x-y plane.

This is like `Quaternions.FrameFromAngularVelocity`, except that the
logarithm of each rotor is constrained to have no z component.  The time
derivative of the logarithm is given by
`Quaternions.FrameFromAngularVelocity_2D_Integrand`, which adjusts the
angular velocity to maintain that constraint.  The frame starts as the
identity at ``T[0]``.

Parameters
----------
Omega : array_like of Quaternion
    The angular-velocity vector at each time, as a quaternion with zero
    scalar part.
T : array_like of float
    The corresponding times.

Returns
-------
numpy.ndarray of Quaternion
    The rotors of the frame at the times ``T``.)doc"},
    {"FrameFromAngularVelocity_Integrand", R"doc(Return the time derivative of the logarithm of a rotating frame's rotor.

If the frame's rotor is :math:`R = e^{\mathfrak{r}}`, with angular velocity
:math:`\vec{\Omega}` defined by :math:`\dot{R} = \vec{\Omega} R/2`, this
function returns :math:`\dot{\mathfrak{r}}`.  Integrating it gives the
logarithm of the frame, which can then be exponentiated to give the rotor.

Parameters
----------
r : array_like of float, shape (3,)
    The vector part of the logarithm of the current rotor.
Omega : array_like of float, shape (3,)
    The instantaneous angular-velocity vector of the frame.

Returns
-------
numpy.ndarray of float, shape (3,)
    The time derivative of ``r``.)doc"},
    {"UnflipRotors", R"doc(Remove sign flips from a time series of rotors.

A rotor and its negative represent the same rotation, so the sign of a rotor
is often arbitrary, and discontinuous flips in that sign cause problems for
rotor-valued functions of time.  This function negates each rotor whose
distance from the previous output rotor exceeds ``discont``, so that
successive rotors are within ``discont`` of each other.

Parameters
----------
R : array_like of Quaternion
    The rotors.
discont : float, default 1.4142135623730951
    The largest acceptable chordal distance between successive rotors; the
    default is :math:`\sqrt{2}`.

Returns
-------
numpy.ndarray of Quaternion
    The rotors, with signs made continuous.)doc"},
    {"RDelta", R"doc(Return the rotors ``R1 * R2.inverse()`` relating two frames.

If ``IndexOfFiducialTime`` is not negative, the result is instead
``R1 * Offset * R2.inverse()``, where ``Offset`` is
``R1[I].inverse() * R2[I]`` for ``I = IndexOfFiducialTime``, which ensures
that the result is 1 at that index.

Parameters
----------
R1 : array_like of Quaternion
    The first rotors.
R2 : array_like of Quaternion
    The second rotors, of the same length as ``R1``.
IndexOfFiducialTime : int, default -1
    Index of the time at which the difference is set to the identity, or a
    negative number for no offset.

Returns
-------
numpy.ndarray of Quaternion
    The difference rotors.)doc"},
    {"Squad", R"doc(Interpolate a time series of rotors with Squad.

Squad is a version of cubic-spline interpolation designed for unit
quaternions, which gives more accurate, smoother, and more physical rotations
than other forms of interpolation.  The input rotors are assumed to have unit
magnitude, and the output rotors are passed through `Quaternions.UnflipRotors`.
Times in ``tOut`` outside the range of ``tIn`` raise an error.

Parameters
----------
RIn : array_like of Quaternion
    The rotors.
tIn : array_like of float
    The corresponding times.
tOut : array_like of float
    The times to which the rotors are interpolated.

Returns
-------
numpy.ndarray of Quaternion
    The interpolated rotors.)doc"},
    {"FrameAngularVelocity", R"doc(Return the angular velocity of a frame given by rotors.

The angular velocity is :math:`2 \dot{R} \bar{R}`, with the derivative
computed by `Quaternions.QuaternionDerivative`.  For unit rotors, each result
has zero scalar part, and its vector part is the angular-velocity vector in
the fixed frame.

Parameters
----------
f : array_like of Quaternion
    The rotors of the frame.
t : array_like of float
    The corresponding times.

Returns
-------
numpy.ndarray of Quaternion
    The angular velocity at each time.)doc"},
    {"BoostRotor", R"doc(Return a rotor taking a direction to its boosted version.

This function returns a rotor :math:`R_b` that takes the input direction
:math:`\hat{w}`, considered as a direction on the future null sphere, to the
direction that it has after the boost.  The rotor therefore depends on both
the boost and the direction.  It is the identity when the velocity is
negligible or parallel to the direction.

Parameters
----------
ThreeVelocity : array_like of float, shape (3,)
    Velocity of the new frame relative to the current frame, in units where
    the speed of light is 1.
ThreeVectorToBeBoosted : array_like of float, shape (3,)
    The direction :math:`\hat{w}` to be boosted, which is normalized first.
    A zero vector raises an error.

Returns
-------
Quaternion
    The rotor :math:`R_b`.)doc"},
    {"ApproximateMeanRotor", R"doc(Return an approximate mean of a set of rotors.

The mean is the normalized sum of the rotors ``R`` or, if the times ``t`` are
given, the normalized trapezoidal integral of ``R`` over the times from ``t1``
to ``t2``.  This minimizes the sum, or integral, of the squared chordal
distance to the mean, rather than of the squared intrinsic distance that
`Quaternions.MeanRotor` minimizes.  The two agree in the limit that the rotors
are all equal, and this function requires no numerical optimization.  The
rotors should be continuous in sign, as `Quaternions.UnflipRotors` makes them.

Parameters
----------
R : array_like of Quaternion
    The rotors to be averaged.
t : array_like of float, optional
    The corresponding times.  If omitted, the rotors are weighted equally.
t1 : float, default -1e300
    Earliest time included in the average.
t2 : float, default 1e300
    Latest time included in the average.

Returns
-------
Quaternion
    The mean rotor.

References
----------
The discrete formula is from Moakher, "Means and averaging in the group of
rotations", SIAM J. Matrix Anal. Appl. 24, 1 (2002); the version with times
is its continuous limit.)doc"},
    {"ApproximateOptimalAlignmentRotor", R"doc(Return the rotor ``R_delta`` for which ``R_delta * Rb`` best matches ``Ra``.

The result is the `Quaternions.ApproximateMeanRotor` of
``RDelta(Ra, Rb)``, averaged over the times from ``t1`` to ``t2``.
Multiplying ``Rb`` on the left by it gives a new time series that is
approximately as close to ``Ra`` as possible.

Parameters
----------
Ra : array_like of Quaternion
    The fixed rotors.
Rb : array_like of Quaternion
    The rotors to be aligned to ``Ra``.
t : array_like of float
    The times of both ``Ra`` and ``Rb``.
t1 : float, default -1e300
    Earliest time included in the average.
t2 : float, default 1e300
    Latest time included in the average.

Returns
-------
Quaternion
    The rotor ``R_delta``.)doc"},
    {"ApproximateOptimalAlignment", R"doc(Find a time shift and rotor that approximately align two rotor time series.

The time offset ``deltat`` is optimized numerically, and for each trial offset
the rotor ``R_delta`` is found analytically by
`Quaternions.ApproximateOptimalAlignmentRotor`.  Together they make
``R_delta * Rb``, with its times shifted to ``tb + deltat``, as close as
possible to ``Ra`` over the interval from ``t1`` to ``t2``, as measured by the
integral of the squared intrinsic distance.  The times ``t1`` and ``t2`` are
fixed with respect to ``ta``, and must lie within the ranges of both ``ta``
and ``tb``.  `Quaternions.OptimalAlignment` refines the result.

Parameters
----------
t1 : float
    Beginning of the alignment interval.
t2 : float
    End of the alignment interval.
Ra : array_like of Quaternion
    The fixed rotors.
ta : array_like of float
    The times of ``Ra``.
Rb : array_like of Quaternion
    The rotors to be aligned to ``Ra``.
tb : array_like of float
    The times of ``Rb``.

Returns
-------
deltat : float
    The optimal time offset.
R_delta : Quaternion
    The optimal rotor.)doc"},
    {"OptimalAlignment", R"doc(Find the time shift and rotor that best align two rotor time series.

Starting from the result of `Quaternions.ApproximateOptimalAlignment`, the
time offset ``deltat`` and all three rotational degrees of freedom of the
rotor ``R_delta`` are optimized together, with the Nelder-Mead simplex
algorithm of GSL.  The quantity minimized is the integral, over the interval
from ``t1`` to ``t2``, of the squared intrinsic distance between ``Ra`` and
``R_delta * Rb``, with the times of ``Rb`` shifted to ``tb + deltat``.  The
times ``t1`` and ``t2`` are fixed with respect to ``ta``, and must lie within
the ranges of both ``ta`` and ``tb``.

Parameters
----------
t1 : float
    Beginning of the alignment interval.
t2 : float
    End of the alignment interval.
Ra : array_like of Quaternion
    The fixed rotors.
ta : array_like of float
    The times of ``Ra``.
Rb : array_like of Quaternion
    The rotors to be aligned to ``Ra``.
tb : array_like of float
    The times of ``Rb``.

Returns
-------
deltat : float
    The optimal time offset.
R_delta : Quaternion
    The optimal rotor.)doc"},
    {"MeanRotor", R"doc(Return the mean of a set of rotors.

The mean minimizes the sum of the squared intrinsic distances between the
rotors ``R`` and the mean or, if the times ``t`` are given, the integral of
the squared intrinsic distance over the times from ``t1`` to ``t2``.  The
minimum is found with the Nelder-Mead simplex algorithm of GSL, starting from
`Quaternions.ApproximateMeanRotor`.  This mean is described in
<https://arxiv.org/abs/1302.2919>.

Parameters
----------
R : array_like of Quaternion
    The rotors to be averaged.
t : array_like of float, optional
    The corresponding times.  If omitted, the rotors are weighted equally.
t1 : float, default -1e300
    Earliest time included in the average; the default means ``t[0]``.
t2 : float, default 1e300
    Latest time included in the average; the default means ``t[-1]``.

Returns
-------
Quaternion
    The mean rotor.)doc"},
    {"hunt", R"doc(Return the largest index ``j`` with ``x[j] <= x_i`` in a sorted array.

The search starts from the guess ``i`` and expands outward, as in the
Numerical Recipes routine of the same name, so it is fast when the guess is
close.  If ``x_i`` is less than ``x[0]``, the result is 0, and if ``x_i`` is at
least ``x[-1]``, the result is ``len(x) - 1``.  An empty array raises an
error.

Parameters
----------
x : array_like of float
    The sorted array.
x_i : float
    The value to locate.
i : int, default 0
    Initial guess for the index.

Returns
-------
int
    The index.)doc"},
    {"huntRight", R"doc(Return one more than the index found by `Quaternions.hunt`.

For ``x[0] <= x_i < x[-1]``, this is the smallest index ``j`` with
``x[j] > x_i``.  If ``x_i`` is less than ``x[0]``, the result is 1, and if
``x_i`` is at least ``x[-1]``, the result is ``len(x)``.

Parameters
----------
x : array_like of float
    The sorted array.
x_i : float
    The value to locate.
i : int, default 0
    Initial guess for the index.

Returns
-------
int
    The index.)doc"},
    {"Matrix", R"doc(Real matrix stored in a GSL matrix.

This class is a thin wrapper around a GSL matrix of floats.  It is used as the
input to the eigen-system helpers, such as `DominantPrincipalAxis` and
`Eigensystem`, most of which require a symmetric :math:`3\times3` matrix.  It
is rarely needed directly.  The element in row ``row`` and column ``col`` is
returned by ``M(row, col)``, and the method `Matrix.numpy` returns a copy of
all the elements as a numpy array.  The difference ``M1 - M2`` of two matrices
of the same shape is taken element by element, and the product ``M * b`` with
a sequence ``b`` of ``M.ncols()`` floats is the matrix-vector product, returned
as a numpy array.

A matrix with zero rows or zero columns is empty: it has no storage, and it
reports zero rows and zero columns.  An empty matrix may be copied, printed,
converted to a numpy array of shape ``(0, 0)``, resized, and cleared.

A `Matrix` may be constructed in any of these ways:

``Matrix()``
    An empty matrix, which may be given a size with `Matrix.resize`.
``Matrix(rows, cols)``
    A matrix of zeros with ``rows`` rows and ``cols`` columns.
``Matrix(rows, cols, a)``
    A matrix with ``rows`` rows and ``cols`` columns, each of whose elements is
    the float ``a``.
``Matrix(M)``
    A copy of the matrix ``M``.
``Matrix(DataIn)``
    A matrix with the elements of ``DataIn``, a sequence of rows, each of which
    is a sequence of floats.  The number of columns is the length of the first
    row, and the matrix is empty if ``DataIn`` or its first row is empty.)doc"},
    {"Matrix::set", R"doc(Set one element of this matrix.

Parameters
----------
row : int
    Index of the row.
col : int
    Index of the column.
value : float
    New value of the element.

Returns
-------
Matrix
    This object, modified in place.)doc"},
    {"Matrix::nrows", R"doc(Return the number of rows of this matrix.

Returns
-------
int
    The number of rows.)doc"},
    {"Matrix::ncols", R"doc(Return the number of columns of this matrix.

Returns
-------
int
    The number of columns.)doc"},
    {"Matrix::resize", R"doc(Change the shape of this matrix, discarding its contents.

The previous elements are not preserved; every element of the resized matrix
is set to ``value``.  If either new dimension is zero, the matrix becomes
empty.

Parameters
----------
newNRows : int
    New number of rows.
newNCols : int
    New number of columns.
value : float, default 0.0
    Value given to every element.)doc"},
    {"Matrix::clear", R"doc(Release the storage of this matrix, leaving it empty.

Afterward the matrix has zero rows and zero columns, and it may be given a new
size with `Matrix.resize`.  Clearing an empty matrix does nothing.)doc"},
    {"Unwrap", R"doc(Unwrap a phase so that it is (roughly) continuous.

Multiples of :math:`2\pi` are added to the elements of ``In`` so that the
difference between successive elements of the result lies between
:math:`-\pi` and :math:`\pi`.  This is the algorithm of MATLAB's ``unwrap``.

Parameters
----------
In : array_like of float
    The phase, in radians.

Returns
-------
numpy.ndarray of float
    The unwrapped phase, with the same length as ``In``.)doc"},
    {"ScalarIntegral", R"doc(Integrate a scalar function of time by the trapezoidal rule.

The integral is cumulative: element ``i`` of the result is the integral from
``t[0]`` to ``t[i]``, so that the first element is zero.

Parameters
----------
fdot : array_like of float
    Values of the integrand.
t : array_like of float
    Times corresponding to the elements of ``fdot``; the two must have the
    same length.

Returns
-------
numpy.ndarray of float
    The integral at each time.

See Also
--------
CumulativeScalarIntegral : The integral over the whole time span.)doc"},
    {"CumulativeScalarIntegral", R"doc(Integrate a scalar function over its whole time span by the trapezoidal rule.

Despite the name, only the total is returned: the result is the integral from
``t[0]`` to ``t[-1]``, which is the last element of the result of
`ScalarIntegral`.

Parameters
----------
fdot : array_like of float
    Values of the integrand.
t : array_like of float
    Times corresponding to the elements of ``fdot``; the two must have the
    same length.

Returns
-------
float
    The integral over the whole time span.)doc"},
    {"ScalarDerivative", R"doc(Differentiate a real function of time by finite differences.

The times need not be uniformly spaced.  When there are at least five points,
the fourth-order five-point formula of Eq. (A 5b) of Bowen and Smith,
"Derivative formulas and errors for non-uniformly spaced points", is used.
With three or four points a three-point formula is used, and with two points
the derivative is the slope of the line through them.  Fewer than two points
raise an error.

Parameters
----------
f : array_like of float
    Values of the function.
t : array_like of float
    Times corresponding to the elements of ``f``; the two must have the same
    length.

Returns
-------
numpy.ndarray of float
    The derivative at each time.)doc"},
    {"ComplexDerivative", R"doc(Differentiate a complex function of time by finite differences.

This is the complex counterpart of `ScalarDerivative`, and uses the same
formulas.

Parameters
----------
f : array_like of complex
    Values of the function.
t : array_like of float
    Times corresponding to the elements of ``f``; the two must have the same
    length.

Returns
-------
numpy.ndarray of complex
    The derivative at each time.)doc"},
    {"VectorIntegral", R"doc(Integrate a vector function of time by the trapezoidal rule.

Each component is integrated separately, as in `ScalarIntegral`.  The integral
is cumulative: row ``i`` of the result is the integral from ``t[0]`` to
``t[i]``, so that the first row is zero.

Parameters
----------
fdot : array_like of float, shape (N, M)
    Values of the integrand, with time along the first axis.
t : array_like of float, shape (N,)
    Times corresponding to the rows of ``fdot``.

Returns
-------
numpy.ndarray of float, shape (N, M)
    The integral at each time.)doc"},
    {"CumulativeVectorIntegral", R"doc(Integrate a vector function over its whole time span by the trapezoidal rule.

Despite the name, only the total is returned: the result is the integral from
``t[0]`` to ``t[-1]``, which is the last row of the result of
`VectorIntegral`.

Parameters
----------
fdot : array_like of float, shape (N, M)
    Values of the integrand, with time along the first axis.
t : array_like of float, shape (N,)
    Times corresponding to the rows of ``fdot``.

Returns
-------
numpy.ndarray of float, shape (M,)
    The integral over the whole time span.)doc"},
    {"Intersection", R"doc(Return a time series covering the intersection of two time series.

The result starts at the latest of ``t1[0]``, ``t2[0]``, and ``MinTime``, and
ends no later than the earliest of ``t1[-1]``, ``t2[-1]``, and ``MaxTime``.
The step from each time to the next is the smaller of the local steps of
``t1`` and ``t2`` at that time, or ``MinStep`` if that is larger.  Both inputs
are assumed to be strictly increasing.  An error is raised if either input is
empty or the intersection is empty.

Parameters
----------
t1 : array_like of float
    The first time series.
t2 : array_like of float
    The second time series.
MinStep : float, default 0.005
    Smallest time step in the result.
MinTime : float, default -1e300
    Earliest time allowed in the result.
MaxTime : float, default 1e300
    Latest time allowed in the result.

Returns
-------
numpy.ndarray of float
    The new time series.

See Also
--------
Union)doc"},
    {"Union", R"doc(Return a time series covering the union of two time series.

The result runs from the earlier of ``t1[0]`` and ``t2[0]`` to the later of
``t1[-1]`` and ``t2[-1]``.  Where only one input covers a time, the step from
that time to the next is the local step of that input; where both do, it is
the smaller of the two local steps.  In either case the step is never smaller
than ``MinStep``.  Both inputs are assumed to be strictly increasing.  If one
input is empty, the other is returned; if the two do not overlap, an error is
raised.

Parameters
----------
t1 : array_like of float
    The first time series.
t2 : array_like of float
    The second time series.
MinStep : float, default 0.005
    Smallest time step in the result.

Returns
-------
numpy.ndarray of float
    The new time series.

See Also
--------
Intersection)doc"},
    {"DominantPrincipalValue", R"doc(Return the largest eigenvalue of each of a sequence of symmetric matrices.

Each matrix must be a symmetric :math:`3\times3` `Matrix`; only the shape of
the first is checked.  The eigenvalues are computed with GSL from the
diagonal and lower triangle of each matrix.

Parameters
----------
M : list of Matrix
    The matrices.

Returns
-------
numpy.ndarray of float
    The largest eigenvalue of each matrix.)doc"},
    {"SubordinatePrincipalValue", R"doc(Return the smallest eigenvalue of each of a sequence of symmetric matrices.

Each matrix must be a symmetric :math:`3\times3` `Matrix`; only the shape of
the first is checked.  The eigenvalues are computed with GSL from the
diagonal and lower triangle of each matrix.

Parameters
----------
M : list of Matrix
    The matrices.

Returns
-------
numpy.ndarray of float
    The smallest eigenvalue of each matrix.)doc"},
    {"DominantPrincipalAxis", R"doc(Return the eigenvector of a symmetric matrix with the largest eigenvalue.

Parameters
----------
M : Matrix
    A symmetric :math:`3\times3` matrix.  It is modified; see the notes.

Returns
-------
numpy.ndarray of float, shape (3,)
    The normalized eigenvector belonging to the largest eigenvalue.

Notes
-----
The eigen-system is computed by GSL's ``gsl_eigen_symmv``, which uses only the
diagonal and lower triangle of ``M`` and overwrites them.  Pass a copy if the
matrix will be needed afterward.)doc"},
    {"SubordinatePrincipalAxis", R"doc(Return the eigenvector of a symmetric matrix with the smallest eigenvalue.

Parameters
----------
M : Matrix
    A symmetric :math:`3\times3` matrix.  It is modified; see the notes.

Returns
-------
numpy.ndarray of float, shape (3,)
    The normalized eigenvector belonging to the smallest eigenvalue.

Notes
-----
The eigen-system is computed by GSL's ``gsl_eigen_symmv``, which uses only the
diagonal and lower triangle of ``M`` and overwrites them.  Pass a copy if the
matrix will be needed afterward.)doc"},
    {"Eigenvalues", R"doc(Return the eigenvalues of a symmetric matrix in decreasing order.

Parameters
----------
M : Matrix
    A symmetric :math:`3\times3` matrix.  It is modified; see the notes.

Returns
-------
numpy.ndarray of float, shape (3,)
    The eigenvalues, from largest to smallest.

Notes
-----
The eigen-system is computed by GSL's ``gsl_eigen_symmv``, which uses only the
diagonal and lower triangle of ``M`` and overwrites them.  Pass a copy if the
matrix will be needed afterward.)doc"},
    {"Eigenvectors", R"doc(Return the eigenvectors of a symmetric matrix.

The eigenvectors are normalized and ordered by decreasing eigenvalue, and are
returned one after another in a flat array: elements ``3*j`` through
``3*j+2`` are the eigenvector with the ``j``-th largest eigenvalue.

Parameters
----------
M : Matrix
    A symmetric :math:`3\times3` matrix.  It is modified; see the notes.

Returns
-------
numpy.ndarray of float, shape (9,)
    The three eigenvectors, concatenated.

Notes
-----
The eigen-system is computed by GSL's ``gsl_eigen_symmv``, which uses only the
diagonal and lower triangle of ``M`` and overwrites them.  Pass a copy if the
matrix will be needed afterward.)doc"},
    {"Eigensystem", R"doc(Return the eigenvectors and eigenvalues of a symmetric matrix.

The eigenpairs are ordered by decreasing eigenvalue, and are returned one after
another in a flat array: elements ``4*j`` through ``4*j+2`` are the normalized
eigenvector with the ``j``-th largest eigenvalue, and element ``4*j+3`` is that
eigenvalue.

Parameters
----------
M : Matrix
    A symmetric :math:`3\times3` matrix.  It is modified; see the notes.

Returns
-------
numpy.ndarray of float, shape (12,)
    The three eigenpairs, concatenated.

Notes
-----
The eigen-system is computed by GSL's ``gsl_eigen_symmv``, which uses only the
diagonal and lower triangle of ``M`` and overwrites them.  Pass a copy if the
matrix will be needed afterward.)doc"},
    {"Determinant", R"doc(Return the determinant of a :math:`3\times3` matrix.

Parameters
----------
M : Matrix
    A :math:`3\times3` matrix; any other shape raises an error.

Returns
-------
float
    The determinant.)doc"},
    {"ScriFunctor", R"doc(Base class for real-valued functions of a rotor, used to fill a `DataGrid`.

To define such a function, subclass `ScriFunctor` in Python and give the
subclass a ``__call__(self, R)`` method that takes a `Quaternion` and returns a
float.  (A subclass that defines ``__init__`` must call
``super().__init__()``.)  An instance is then passed as ``f`` to
``DataGrid(Spin, N_theta, N_phi, v, f)``, which calls it once for each point of
a boosted grid.  The rotor ``R`` takes :math:`\hat{z}` into the point at which
the function is to be evaluated, and takes :math:`\hat{x} + i \hat{y}` into
the vector :math:`m` (up to normalization) needed for spin-weighted fields.
The ``__call__`` method of the base class returns zero.

A `ScriFunctor` is constructed as ``ScriFunctor()``, with no arguments.)doc"},
    {"Boost", R"doc(Return a rotor taking a direction into its boosted version.

The returned rotor :math:`R_b` takes the unit vector :math:`\hat{n}` (the
normalized ``n``), regarded as a point on the future null sphere, into its
image under a boost with velocity ``v``.  The rotor therefore depends on both
the direction and the boost.  It is a rotation about the axis
:math:`\hat{v} \times \hat{n}` that changes the angle :math:`\Theta` between
:math:`\hat{v}` and :math:`\hat{n}` into
:math:`\Theta' = 2 \arctan\left(e^{\alpha} \tan(\Theta/2)\right)`, where
:math:`\alpha` is the rapidity of the boost.  The identity is returned when
``v`` is negligibly small or ``n`` is parallel to ``v``.

Parameters
----------
v : array_like of float, shape (3,)
    Three-velocity of the new frame relative to the current frame, in units
    where :math:`c = 1`.
n : array_like of float, shape (3,)
    Direction to be boosted.  It is normalized, and must not be zero.

Returns
-------
Quaternion
    The rotor :math:`R_b`.)doc"},
    {"DataGrid", R"doc(Complex spin-weighted data on an equiangular grid on the sphere.

Given integers ``N_theta`` and ``N_phi``, the data are recorded on a
"rectangular" grid of ``N_theta * N_phi`` points.  The polar angle runs from 0
to :math:`\pi` inclusive in ``N_theta`` equal steps, so that there are
``N_phi`` points at each pole, and the azimuthal angle runs from 0 to
:math:`2\pi` exclusive in ``N_phi`` equal steps.  The point with polar index
``i_theta`` and azimuthal index ``i_phi`` is element
``i_theta * N_phi + i_phi``.  These objects are mainly a computational tool
for pointwise multiplication; otherwise `Modes` is the preferable
representation.

The arithmetic operators ``+``, ``-``, ``*``, and ``/`` act point by point on
two grids of the same dimensions; a product adds the spin weights, and a
quotient subtracts them.  A float ``a`` may also appear on the left, as in
``a * D``, ``a / D``, ``a + D``, and ``a - D``, each of which acts on every
point of ``D``.

A `DataGrid` may be constructed in any of these ways:

``DataGrid(size=0)``
    A grid of ``size`` zeros with spin weight 0, whose dimensions are both set
    to the integer part of the square root of ``size``.
``DataGrid(A)``
    A copy of the grid ``A``.
``DataGrid(Spin, N_theta, N_phi, D)``
    A grid with spin weight ``Spin``, dimensions ``N_theta`` and ``N_phi``,
    and the complex values ``D``, whose length must be ``N_theta * N_phi``.
``DataGrid(M, N_theta=0, N_phi=0)``
    The values on the grid of the function whose spin-weighted
    spherical-harmonic modes are the `Modes` object ``M``.  Each dimension is
    at least ``2 * M.EllMax() + 1``.
``DataGrid(M, v, N_theta=0, N_phi=0)``
    The values of ``M`` on the grid of a frame boosted by the three-velocity
    ``v``.  Each point of the grid is mapped by the rotor from `Boost` before
    ``M`` is evaluated there with `Modes.EvaluateAtPoint`.  Each dimension is
    at least ``2 * M.EllMax() + 1``.
``DataGrid(Spin, N_theta, N_phi, v, f)``
    The values of the `ScriFunctor` ``f`` on the grid of a frame boosted by the
    three-velocity ``v``, with spin weight ``Spin`` and dimensions ``N_theta``
    and ``N_phi``.)doc"},
    {"DataGrid::SetSpin", R"doc(Set the spin weight of these data, without changing the values.

Parameters
----------
ess : int
    The new spin weight.

Returns
-------
DataGrid
    This object, modified in place.)doc"},
    {"DataGrid::SetNTheta", R"doc(Set the number of polar points, without changing the values.

Parameters
----------
N_theta : int
    The new number of points.

Returns
-------
DataGrid
    This object, modified in place.)doc"},
    {"DataGrid::SetNPhi", R"doc(Set the number of azimuthal points, without changing the values.

Parameters
----------
N_phi : int
    The new number of points.

Returns
-------
DataGrid
    This object, modified in place.)doc"},
    {"DataGrid::size", R"doc(Return the number of values stored in this grid.

Returns
-------
int
    The number of values, which is also given by ``len``.)doc"},
    {"DataGrid::Spin", R"doc(Return the spin weight of these data.

Returns
-------
int
    The spin weight.)doc"},
    {"DataGrid::N_theta", R"doc(Return the number of points in the polar direction.

Returns
-------
int
    The number of points.)doc"},
    {"DataGrid::N_phi", R"doc(Return the number of points in the azimuthal direction.

Returns
-------
int
    The number of points.)doc"},
    {"DataGrid::Data", R"doc(Return a copy of the values on this grid.

Returns
-------
numpy.ndarray of complex
    The values, in the order described in `DataGrid`.)doc"},
    {"DataGrid::pow", R"doc(Raise the value at each point of this grid to an integer power.

The spin weight of the result is the same as that of this grid.

Parameters
----------
p : int
    The exponent.

Returns
-------
DataGrid
    The new grid.)doc"},
    {"ConformalFactorGrid", R"doc(Construct a grid with the conformal factor of a boost at each point.

The conformal factor is :math:`K = 1 / \left[\gamma (1 - \vec{v} \cdot
\hat{n})\right]`, where :math:`\hat{n}` is the direction of the grid point and
:math:`\gamma` is the Lorentz factor of ``v``.  The grid has spin weight 0.

Parameters
----------
v : array_like of float, shape (3,)
    Three-velocity of the boost.
n_theta : int
    Number of points in the polar direction.
n_phi : int
    Number of points in the azimuthal direction.

Returns
-------
DataGrid
    The conformal factor on the grid.)doc"},
    {"InverseConformalFactorGrid", R"doc(Construct a grid with the inverse conformal factor of a boost at each point.

The inverse conformal factor is :math:`1/K = \gamma (1 - \vec{v} \cdot
\hat{n})`, where :math:`\hat{n}` is the direction of the grid point and
:math:`\gamma` is the Lorentz factor of ``v``.  The grid has spin weight 0.

Parameters
----------
v : array_like of float, shape (3,)
    Three-velocity of the boost.
n_theta : int
    Number of points in the polar direction.
n_phi : int
    Number of points in the azimuthal direction.

Returns
-------
DataGrid
    The inverse conformal factor on the grid.

See Also
--------
InverseConformalFactorBoostedGrid, ConformalFactorGrid)doc"},
    {"InverseConformalFactorBoostedGrid", R"doc(Construct a boosted grid with the inverse conformal factor at each point.

This is the same function as in `InverseConformalFactorGrid`,
:math:`1/K = \gamma (1 - \vec{v} \cdot \hat{n})`, but evaluated on the grid
of the frame boosted by ``v``, as in ``DataGrid(Spin, N_theta, N_phi, v, f)``:
each grid point is first mapped by the rotor from `Boost`.  The grid has spin
weight 0.

Parameters
----------
v : array_like of float, shape (3,)
    Three-velocity of the boost.
n_theta : int
    Number of points in the polar direction.
n_phi : int
    Number of points in the azimuthal direction.

Returns
-------
DataGrid
    The inverse conformal factor on the boosted grid.)doc"},
    {"Modes", R"doc(Complex spin-weighted data on the sphere, as spherical-harmonic modes.

The data are the coefficients of the spin-weighted spherical harmonics
:math:`{}_{s}Y_{\ell m}` for every :math:`\ell` from 0 to `Modes.EllMax`,
even for :math:`\ell < |s|`, where they vanish.  They are stored in the order
:math:`(\ell, m) = (0, 0), (1, -1), (1, 0), (1, 1), (2, -2), \ldots`, so that
the mode :math:`(\ell, m)` is element :math:`\ell^2 + \ell + m`.

The arithmetic operators ``+``, ``-``, ``*``, and ``/`` act on two `Modes`
objects.  Sums and differences require equal spin weights, and treat the
missing higher modes of the shorter operand as zero.  Products and quotients
are computed point by point on a `DataGrid` large enough to resolve the
result, which is then truncated to the larger of the two values of
:math:`\ell_{\max}`; a product adds the spin weights, and a quotient subtracts
them.

A `Modes` object may be constructed in any of these ways:

``Modes(size=0)``
    ``size`` zero modes, with spin weight 0 and `Modes.EllMax` 0 (whatever
    the value of ``size``).
``Modes(A)``
    A copy of the `Modes` object ``A``.
``Modes(spin, Data)``
    The complex modes ``Data`` with spin weight ``spin``.  The length of
    ``Data`` must be :math:`(\ell_{\max}+1)^2` for some :math:`\ell_{\max}`.
``Modes(D, L=-1)``
    The modes of the `DataGrid` ``D``, computed with spinsfast.  The value of
    :math:`\ell_{\max}` is the larger of ``L`` and the largest value the grid
    resolves, which is the smaller of ``(D.N_theta() - 1) // 2`` and
    ``(D.N_phi() - 1) // 2``.)doc"},
    {"Modes::SetSpin", R"doc(Set the spin weight of these data, without changing the modes.

Parameters
----------
ess : int
    The new spin weight.

Returns
-------
Modes
    This object, modified in place.)doc"},
    {"Modes::SetEllMax", R"doc(Set the largest ell value of these data, without changing the modes.

Parameters
----------
ell : int
    The new value of :math:`\ell_{\max}`.

Returns
-------
Modes
    This object, modified in place.)doc"},
    {"Modes::size", R"doc(Return the number of modes stored.

Returns
-------
int
    The number of modes, which is also given by ``len``.)doc"},
    {"Modes::Spin", R"doc(Return the spin weight of these data.

Returns
-------
int
    The spin weight.)doc"},
    {"Modes::EllMax", R"doc(Return the largest ell value of these data.

Returns
-------
int
    The value of :math:`\ell_{\max}`.)doc"},
    {"Modes::Data", R"doc(Return a copy of the modes.

Returns
-------
numpy.ndarray of complex
    The modes, in the order described in `Modes`.)doc"},
    {"Modes::pow", R"doc(Raise the function represented by these modes to an integer power.

The power is taken point by point on a `DataGrid` with
``2 * EllMax() * p + 1`` points in each direction, and the result has
:math:`\ell_{\max}` equal to ``EllMax() * p``.  The spin weight of the result
is the same as that of these data.

Parameters
----------
p : int
    The exponent.

Returns
-------
Modes
    The new modes.)doc"},
    {"Modes::bar", R"doc(Return the complex conjugates of these modes, with the spin weight negated.

Each mode coefficient is replaced by its complex conjugate, without any
reordering of the modes.

Returns
-------
Modes
    The new modes.)doc"},
    {"Modes::edth", R"doc(Apply the spin-raising operator :math:`\eth` to these data.

This is the operator defined by Geroch, Held, and Penrose (1973).  It raises
the spin weight by 1, and leaves the boost weight unchanged.  Each mode with
:math:`\ell \geq |s+1|` is multiplied by
:math:`\sqrt{(\ell-s)(\ell+s+1)/2}`, and the others are set to zero.

Returns
-------
Modes
    The new modes, with spin weight :math:`s+1`.

See Also
--------
Modes.edthbar

Notes
-----
In this context the operator differs from the Newman--Penrose :math:`\eth`
only by normalization: the Newman--Penrose operator is :math:`\sqrt{2}` times
this one.

References
----------
R. Geroch, A. Held, and R. Penrose, "A space-time calculus based on pairs of
null directions", J. Math. Phys. 14, 874 (1973).)doc"},
    {"Modes::edthbar", R"doc(Apply the spin-lowering operator :math:`\bar{\eth}` to these data.

This is the operator defined by Geroch, Held, and Penrose (1973).  It lowers
the spin weight by 1, and leaves the boost weight unchanged.  Each mode with
:math:`\ell \geq |s-1|` is multiplied by
:math:`-\sqrt{(\ell+s)(\ell-s+1)/2}`, and the others are set to zero.

Returns
-------
Modes
    The new modes, with spin weight :math:`s-1`.

See Also
--------
Modes.edth

Notes
-----
In this context the operator differs from the Newman--Penrose
:math:`\bar{\eth}` only by normalization: the Newman--Penrose operator is
:math:`\sqrt{2}` times this one.

References
----------
R. Geroch, A. Held, and R. Penrose, "A space-time calculus based on pairs of
null directions", J. Math. Phys. 14, 874 (1973).)doc"},
    {"Modes::edth2edthbar2", R"doc(Apply the operator :math:`\eth^2 \bar{\eth}^2` to these data.

Each mode is multiplied by :math:`(\ell-1)\ell(\ell+1)(\ell+2)`, whatever the
spin weight, which is unchanged.  For data of spin weight 0, this is
:math:`\eth^2 \bar{\eth}^2` with the Newman--Penrose normalization, or four
times the result of applying `Modes.edthbar` twice and then `Modes.edth`
twice.

Returns
-------
Modes
    The new modes.)doc"},
    {"Modes::EvaluateAtPoint", R"doc(Evaluate the function represented by these modes at a point on the sphere.

The point may be given either by its polar and azimuthal angles or by a rotor
``R``.  The rotor might typically be thought of as the rotation taking
:math:`\hat{z}` into the point :math:`(\vartheta, \varphi)`.  More general
rotors are possible, however; they are used to construct a `DataGrid` from a
`Modes` object with a boost.

Parameters
----------
vartheta : float
    Polar angle of the point.
varphi : float
    Azimuthal angle of the point.
R : Quaternion
    Rotor taking :math:`\hat{z}` into the point.

Returns
-------
complex
    The value of the function at the point.)doc"},
    {"vFromOneOverK", R"doc(Derive the three-velocity of a boost from its inverse conformal factor.

Only the :math:`\ell = 0` and :math:`\ell = 1` modes of ``OneOverK`` are used.

Parameters
----------
OneOverK : Modes
    Modes of the inverse conformal factor :math:`1/K`.

Returns
-------
numpy.ndarray of float, shape (3,)
    The three-velocity.)doc"},
    {"SliceOfScriGrid", R"doc(The data on one slice of null infinity, each as a `DataGrid`.

A slice is the set of points of null infinity at one value of the retarded
time.  The data needed to describe its geometry are the Newman--Penrose
curvature scalars :math:`\psi_0, \ldots, \psi_4`, the complex shear
:math:`\sigma` of outgoing null rays, and its time derivative
:math:`\dot{\sigma}`.  Besides being available as attributes, these seven
grids may be accessed by index, in that order, so that ``S[0]`` is ``psi0``
and ``S[6]`` is ``sigmadot``; ``len(S)`` is always 7.  Indexing and attribute
access both return the stored grid, not a copy.  Objects of this class are
returned by `SliceModes.BMSTransformationOnSlice`.

A `SliceOfScriGrid` may be constructed in either of these ways:

``SliceOfScriGrid(size=0)``
    Seven grids constructed as ``DataGrid(size)``, with the spin weights 2,
    1, 0, -1, -2, 2, and 2 of the quantities above.
``SliceOfScriGrid(S)``
    A copy of the slice ``S``.

Attributes
----------
psi0, psi1, psi2, psi3, psi4 : DataGrid
    The Newman--Penrose curvature scalars.
sigma : DataGrid
    The complex shear.
sigmadot : DataGrid
    The time derivative of the complex shear.)doc"},
    {"SliceOfScriModes", R"doc(The data on one slice of null infinity, each as a `Modes` object.

A slice is the set of points of null infinity at one value of the retarded
time.  The data needed to describe its geometry are the Newman--Penrose
curvature scalars :math:`\psi_0, \ldots, \psi_4`, the complex shear
:math:`\sigma` of outgoing null rays, and its time derivative
:math:`\dot{\sigma}`.  Besides being available as attributes, these seven
`Modes` objects may be accessed by index, in that order, so that ``S[0]`` is
``psi0`` and ``S[6]`` is ``sigmadot``; ``len(S)`` is always 7.  Indexing and
attribute access both return the stored object, not a copy.  This is the base
class of `SliceModes`, which is usually the class to use.

A `SliceOfScriModes` may be constructed in either of these ways:

``SliceOfScriModes(size=0)``
    Seven objects constructed as ``Modes(size)``, with the spin weights 2, 1,
    0, -1, -2, 2, and 2 of the quantities above.
``SliceOfScriModes(S)``
    A copy of the slice ``S``.

Attributes
----------
psi0, psi1, psi2, psi3, psi4 : Modes
    The Newman--Penrose curvature scalars.
sigma : Modes
    The complex shear.
sigmadot : Modes
    The time derivative of the complex shear.)doc"},
    {"SliceModes", R"doc(The data on one slice of null infinity, as modes, with derived quantities.

This subclass of `SliceOfScriModes` adds methods for the supermomentum, the
Bondi four-momentum, and the mass on the slice, and for transforming the data
on the slice under the Bondi--Metzner--Sachs (BMS) group.

A `SliceModes` object may be constructed in either of these ways:

``SliceModes(ellMax=0)``
    Seven `Modes` objects of zeros, each with :math:`\ell_{\max}` equal to
    ``ellMax`` and with the spin weights 2, 1, 0, -1, -2, 2, and 2 of
    :math:`\psi_0, \ldots, \psi_4`, :math:`\sigma`, and :math:`\dot{\sigma}`.
``SliceModes(S)``
    A copy of the slice ``S``.)doc"},
    {"SliceModes::EllMax", R"doc(Return the largest ell value of the data on this slice.

Returns
-------
int
    The largest value of `Modes.EllMax` among the seven quantities.)doc"},
    {"SliceModes::Mass", R"doc(Calculate the mass of the system from the four-momentum.

Returns
-------
float
    The mass :math:`\sqrt{p_0^2 - p_1^2 - p_2^2 - p_3^2}`, where :math:`p` is
    the result of `SliceModes.FourMomentum`.)doc"},
    {"SliceModes::FourMomentum", R"doc(Calculate the Bondi four-momentum of the system from the supermomentum.

The four-momentum is given by the :math:`\ell = 0` and :math:`\ell = 1` modes
of the result of `SliceModes.SuperMomentum`.

Returns
-------
numpy.ndarray of float, shape (4,)
    The four-momentum, with the time component first.)doc"},
    {"SliceModes::SuperMomentum", R"doc(Calculate the Moreschi supermomentum on this slice.

The supermomentum is

.. math::

    \Psi = \psi_2 + \sigma \dot{\bar{\sigma}} + \eth^2 \bar{\sigma},

where :math:`\eth` is the operator of `Modes.edth`.

Returns
-------
Modes
    The supermomentum.)doc"},
    {"SliceModes::BMSTransformationOnSlice", R"doc(Apply a BMS transformation within this slice, without moving points in time.

A full BMS transformation needs data from several slices, because the
supertranslation moves points between slices.  This method only transforms
the data within the slice: it evaluates them on the grid of the frame boosted
by ``v``, with the factors of the conformal factor, the mixing of the
curvature scalars caused by the change of tetrad, and the shift of
:math:`\sigma` by :math:`\eth^2 \delta`.  Each point of the resulting grids
can then be interpolated in time to the supertranslated slice, as
`Scri.BMSTransformation` does.  The grids have ``2 * EllMax() + 1`` points in
each direction.

Parameters
----------
u : float
    Retarded time of this slice.
v : array_like of float, shape (3,)
    Three-velocity of the boost relative to the current frame.
delta : Modes
    Spherical-harmonic modes of the supertranslation.

Returns
-------
SliceOfScriGrid
    The transformed data on the boosted grid.)doc"},
    {"SliceModes::MoreschiIteration", R"doc(Find the next iteration of the BMS transformation by Moreschi's algorithm.

This method is meant for a slice that has already been transformed by the BMS
transformation given by :math:`K_i` and :math:`\delta_i`, and is meant to
compute :math:`K_{i+1}` and :math:`\delta_{i+1}` from the data on the slice,
storing them in the arguments.

Parameters
----------
OneOverK_ip1 : Modes
    Receives the inverse conformal factor for the next step.
delta_ip1 : Modes
    Receives the supertranslation for the next step.

Warnings
--------
This method is not finished.  It sets ``OneOverK_ip1`` to the :math:`\ell \leq
1` modes of the supermomentum divided by the mass, leaves ``delta_ip1``
unchanged, and then always raises an error.  `SuperMomenta.MoreschiIteration`
is a working implementation.)doc"},
    {"Scri", R"doc(The asymptotic geometry of null infinity, as a sequence of slices.

A `Scri` object contains the data needed to express the asymptotic geometry
of null infinity and its symmetry transformations.  This geometry is encoded
in the Newman--Penrose curvature scalars :math:`\psi_0, \ldots, \psi_4` and
the complex shear :math:`\sigma` of outgoing null rays (the strain in
gravitational-wave detectors).  The general symmetry transformation is an
element of the Bondi--Metzner--Sachs (BMS) group, which transforms these data
among themselves.

The data are stored as one `SliceModes` object at each time.  ``len(S)`` is
the number of times, and ``S[i]`` is a copy of the slice at time index ``i``,
which may be replaced by assigning to ``S[i]``.

A `Scri` object may be constructed in either of these ways:

``Scri(psi0, psi1, psi2, psi3, psi4, sigma)``
    The data in six `Waveform` objects, which must have the same number of
    times and the same value of `Waveform.EllMax`, and must contain every
    mode :math:`(\ell, m)` up to that value, including :math:`\ell = 0` and
    :math:`\ell = 1`.  The times are taken from ``psi0``, and
    :math:`\dot{\sigma}` is computed by differentiating ``sigma`` in time.
``Scri(S)``
    A copy of the `Scri` object ``S``.)doc"},
    {"Scri::BMSTransformation", R"doc(Apply a constant BMS transformation to the data on null infinity.

A general BMS transformation is a conformal transformation of the sphere
(which encompasses rotations and boosts) together with a supertranslation
(which affects only the time coordinate on null infinity).  Here the conformal
transformation is assumed to be a pure boost with velocity ``v``; the rotation
is taken to be zero.  The supertranslation is given by its (scalar)
spherical-harmonic modes, which are assumed to describe a real function on the
sphere.

The new and old time coordinates are related by :math:`u' = K(u - \delta)`,
where the conformal factor :math:`K` depends on the angle.  The slice
returned is the one on which :math:`u' = 0`, since any other choice can be
absorbed into a translation in time and space.  The data at each grid point
are therefore interpolated, with cubic splines in time, to the original time
:math:`u = u_0 + \delta`, which differs from point to point.  The data are
first evaluated on the grid of the final frame by
`SliceModes.BMSTransformationOnSlice`, and are finally transformed back into
modes.  An error is raised if these times are not all within the time span of
the data.

Parameters
----------
u0 : float
    Time of the slice to be transformed, to which the supertranslation is
    added.
v : array_like of float, shape (3,)
    Three-velocity of the boost relative to the current frame.
delta : Modes
    Spherical-harmonic modes of the supertranslation.

Returns
-------
SliceModes
    The transformed data on the new slice.)doc"},
    {"Scri::NTimes", R"doc(Return the number of time steps.

Returns
-------
int
    The number of slices, which is also given by ``len``.)doc"},
    {"Scri::T", R"doc(Return the times of the slices.

Returns
-------
numpy.ndarray of float
    The times.)doc"},
    {"SuperMomenta", R"doc(The supermomentum of a system as a function of time.

A `SuperMomenta` object stores a `Modes` object for the supermomentum at each
of a sequence of times, as computed by `SliceModes.SuperMomentum`, and is used
to find a BMS frame by Moreschi's algorithm.  ``len(S)`` is the number of
times, and ``S[i]`` is a copy of the supermomentum at time index ``i``, which
may be replaced by assigning to ``S[i]``.

A `SuperMomenta` object may be constructed in any of these ways:

``SuperMomenta(size)``
    ``size`` times, all zero, with empty `Modes` objects.
``SuperMomenta(S)``
    A copy of the `SuperMomenta` object ``S``.
``SuperMomenta(T, psi)``
    The times ``T``, a sequence of floats, and the supermomenta ``psi``, a
    sequence of `Modes` objects.
``SuperMomenta(scri)``
    The times of the `Scri` object ``scri``, and the supermomentum of each of
    its slices.)doc"},
    {"SuperMomenta::NTimes", R"doc(Return the number of time steps.

Returns
-------
int
    The number of times, which is also given by ``len``.)doc"},
    {"SuperMomenta::T", R"doc(Return the times at which the supermomentum is stored.

Returns
-------
numpy.ndarray of float
    The times.)doc"},
    {"SuperMomenta::BMSTransform", R"doc(Return the supermomentum on a constant-time slice of a BMS-transformed frame.

The boost of the transformation is found from ``OneOverK`` by
`vFromOneOverK`.  At each time, the quantity
:math:`(\Psi - \eth^2 \bar{\eth}^2 \delta)/K^3` is evaluated on the grid of
the boosted frame, and the values at each grid point are interpolated, with
cubic splines in time, to the time :math:`u = \delta` at that point.  This is
the slice :math:`u' = 0` of the new frame, which is centered on the
:math:`\ell = 0` mode of ``delta``.  An error is raised if these times are not
all within the time span of the data.

Parameters
----------
OneOverK : Modes
    Modes of the inverse conformal factor :math:`1/K` of the boost.
delta : Modes
    Spherical-harmonic modes of the supertranslation.

Returns
-------
Modes
    The transformed supermomentum.)doc"},
    {"SuperMomenta::MoreschiIteration", R"doc(Take one step of Moreschi's algorithm, updating a BMS transformation in place.

The supermomentum :math:`\Psi` is first transformed by
`SuperMomenta.BMSTransform` with the given ``OneOverK`` and ``delta``, and the
mass :math:`M` is computed from the four-momentum given by its
:math:`\ell \leq 1` modes.  The arguments are then replaced by the next step
of the algorithm.  Each mode of ``delta`` with :math:`\ell \geq 2` becomes
:math:`4/[(\ell-1)\ell(\ell+1)(\ell+2)]` times the corresponding mode of
:math:`\Psi + M K^3`; its :math:`\ell = 0` and :math:`\ell = 1` modes are
unchanged.  The :math:`\ell = 0` mode of ``OneOverK`` becomes
:math:`\Psi_{00} / M`, and its :math:`\ell = 1` modes become
:math:`-\Psi_{1m} / (3M)`, with the sign reversed so that the boost found is
removed by the next step rather than reported.

Parameters
----------
OneOverK : Modes
    Inverse conformal factor, with at least the :math:`\ell \leq 1` modes; it
    is replaced by the next step.
delta : Modes
    Supertranslation; it is replaced by the next step.)doc"},
    {"Waveform", R"doc(A gravitational waveform, decomposed into spin-weighted spherical harmonics.

A `Waveform` stores the complex mode data of a single waveform as functions of
time, together with the :math:`(\ell, m)` values of the modes and the frame
with respect to which the modes are decomposed.  It also records the spin
weight and boost weight of the data, the frame type and data type, whether the
radius and the mass have been scaled out of the data, and a history of the
operations applied to the object.  Its methods provide the manipulations needed
in waveform analysis: slicing and interpolation in time; rotations and
transformations to the coprecessing, corotating, and inertial frames;
alignment, comparison, and hybridization; and evaluation on the sphere.

The frame data are unit quaternions (rotors).  The rotor at each time takes the
standard basis :math:`(x, y, z)` into the basis :math:`(X, Y, Z)` with respect
to which the modes are decomposed.  If there are no frame data, the modes are
decomposed with respect to the standard basis; if there is a single rotor, the
frame is constant; otherwise there is one rotor for each time step.

Two waveforms with the same times and in the same frame may be added or
subtracted with ``+`` and ``-``, which act mode by mode, or multiplied or
divided with ``*`` and ``/``, which act pointwise on the sphere.  The spin
weight of a product is the sum of the spin weights of the factors, and that of
a quotient is their difference.  The result has every mode from
:math:`\ell = |s|`, where :math:`s` is its spin weight, up to the smaller of
the largest :math:`\ell` values of the two inputs; modes of the inputs above
that :math:`\ell` are ignored.  The pointwise values are computed on a grid
that resolves only those modes, so the result is exact only when the product
or quotient of the truncated inputs has no higher modes; otherwise the higher
modes are aliased into the result.  A waveform may also be multiplied or
divided by a float.

A `Waveform` may be constructed in any of these ways:

``Waveform()``
    An empty waveform.
``Waveform(W)``
    A copy of the waveform ``W``, including its history.
``Waveform(FileName, DataFormat)``
    The data in the text file ``FileName``.  Lines at the beginning of the
    file that start with ``#`` or ``%`` are stored in the history.  Each
    remaining line holds a time, followed by two columns for each mode: the
    real and imaginary parts if ``DataFormat`` contains ``"ReIm"``, or the
    magnitude and argument if it contains ``"MagArg"`` (ignoring case).  The
    modes are assumed to be ordered as :math:`(2, -2), (2, -1), \ldots,
    (2, 2), (3, -3), \ldots`; if they are not, the :math:`(\ell, m)` data
    must be corrected with `Waveform.SetLM`.
``Waveform(T, LM, Data)``
    A waveform with the times ``T``, an array of floats; the mode labels
    ``LM``, an array of integers with shape (NModes, 2), each row of which
    holds :math:`\ell` and :math:`m`; and the complex mode data ``Data``,
    with shape (NModes, NTimes).

Every constructor except the copy constructor sets the spin weight to -2 and
the boost weight to -1, and leaves the frame type and the data type unknown.)doc"},
    {"Waveform::CopyWithoutData", R"doc(Return a copy of this waveform without its times, frame, modes, or data.

The copy has the spin weight, boost weight, history, version history, frame
type, data type, and scaling flags of this waveform, but its times, frame data,
:math:`(\ell, m)` data, and mode data are empty.

Returns
-------
Waveform
    The copy.)doc"},
    {"Waveform::SliceOfTimeIndices", R"doc(Return a copy of the waveform restricted to a range of time indices.

The copy contains the time steps with indices from ``i_t_a`` up to, but not
including, ``i_t_b``.

Parameters
----------
i_t_a : int
    Index of the first time step to include.
i_t_b : int, default 0
    Index one past the last time step to include.  If it is 0, only the time
    step ``i_t_a`` is included.

Returns
-------
Waveform
    The slice of the waveform.)doc"},
    {"Waveform::SliceOfTimeIndicesWithEll2", R"doc(Return a copy of the waveform restricted to time indices and :math:`\ell = 2`.

The copy contains only the five modes with :math:`\ell = 2`, ordered by
:math:`m`, at the time steps with indices from ``i_t_a`` up to, but not
including, ``i_t_b``.

Parameters
----------
i_t_a : int
    Index of the first time step to include.
i_t_b : int, default 0
    Index one past the last time step to include.  If it is 0, only the time
    step ``i_t_a`` is included.

Returns
-------
Waveform
    The slice of the waveform.)doc"},
    {"Waveform::SliceOfTimeIndicesWithoutModes", R"doc(Return a copy of the waveform restricted to time indices, without mode data.

The copy contains the times and frame data at the time steps with indices from
``i_t_a`` up to, but not including, ``i_t_b``, but no modes.  Only frame data
with one rotor for each time step are copied; a constant frame is dropped.

Parameters
----------
i_t_a : int
    Index of the first time step to include.
i_t_b : int, default 0
    Index one past the last time step to include.  If it is 0, only the time
    step ``i_t_a`` is included.

Returns
-------
Waveform
    The slice of the waveform.)doc"},
    {"Waveform::SliceOfTimes", R"doc(Return a copy of the waveform restricted to an interval of time.

The copy begins at the last time step at or before ``t_a`` (or at the first
time step, if there is none) and ends at the last time step at or before
``t_b``.  With the default arguments, the whole waveform is copied.

Parameters
----------
t_a : float, default -1e300
    Beginning of the interval.
t_b : float, default 1e300
    End of the interval.

Returns
-------
Waveform
    The slice of the waveform.)doc"},
    {"Waveform::SliceOfTimesWithEll2", R"doc(Return a copy of the waveform restricted to a time span and :math:`\ell = 2`.

The copy contains only the five modes with :math:`\ell = 2`, ordered by
:math:`m`.  It begins at the last time step at or before ``t_a`` (or at the
first time step, if there is none) and ends at the last time step at or before
``t_b``.

Parameters
----------
t_a : float, default -1e300
    Beginning of the interval.
t_b : float, default 1e300
    End of the interval.

Returns
-------
Waveform
    The slice of the waveform.)doc"},
    {"Waveform::SliceOfTimesWithoutModes", R"doc(Return a copy of the waveform restricted to a time interval, without modes.

The copy contains the times and frame data, but no modes.  It begins at the
last time step at or before ``t_a`` (or at the first time step, if there is
none) and ends at the last time step at or before ``t_b``.  Only frame data
with one rotor for each time step are copied; a constant frame is dropped.

Parameters
----------
t_a : float, default -1e300
    Beginning of the interval.
t_b : float, default 1e300
    End of the interval.

Returns
-------
Waveform
    The slice of the waveform.)doc"},
    {"Waveform::Interpolate", R"doc(Interpolate the waveform to a new set of times.

The real and imaginary parts of each mode are interpolated with cubic splines,
and frame data with one rotor for each time step are interpolated with
`Quaternions.Squad`.  All other properties, including the spin weight, boost
weight, and version history, are copied.  This waveform is not changed.

Parameters
----------
NewTime : array_like of float
    New times to which the data are interpolated, in increasing order.
AllowTimesOutsideCurrentDomain : bool, default False
    If true, the modes are set to zero at times outside the current time
    domain, and frame data with one rotor for each time step are held at the
    rotor at the nearer end of that domain (the first rotor at earlier times,
    and the last at later times); if false, such times raise an error.

Returns
-------
Waveform
    The interpolated waveform.)doc"},
    {"Waveform::InterpolateInPlace", R"doc(Interpolate the waveform to a new set of times, in place.

The real and imaginary parts of each mode are interpolated with cubic splines,
and frame data with one rotor for each time step are interpolated with
`Quaternions.Squad`.  Times outside the current time domain raise an error.

Parameters
----------
NewTime : array_like of float
    New times to which the data are interpolated.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::DropTimesOutside", R"doc(Remove the data at times outside an interval, in place.

The time steps that remain are those with ``ta < t <= tb``, except that the
first time step is always removed, even when ``ta`` precedes it.

Parameters
----------
ta : float
    Beginning of the interval.
tb : float
    End of the interval.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::DropEllModes", R"doc(Remove the modes with the given values of :math:`\ell`, in place.

Parameters
----------
EllModesToDrop : array_like of int
    Values of :math:`\ell` whose modes are removed.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::KeepOnlyEllModes", R"doc(Remove all modes except those with the given values of :math:`\ell`, in place.

Parameters
----------
EllModesToKeep : array_like of int
    Values of :math:`\ell` whose modes are kept.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::KeepOnlyEll2", R"doc(Remove all modes except those with :math:`\ell = 2`, in place.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetSpinWeight", R"doc(Set the spin weight of the data.

Only the recorded spin weight changes; the data are not altered.

Parameters
----------
NewSpinWeight : int
    The new spin weight.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetBoostWeight", R"doc(Set the boost weight of the data.

Only the recorded boost weight changes; the data are not altered.

Parameters
----------
NewBoostWeight : int
    The new boost weight.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::AppendHistory", R"doc(Append text to the history of the waveform.

The text is appended exactly as given; no newline is added.

Parameters
----------
Hist : str
    Text to append.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetHistory", R"doc(Replace the history of the waveform.

Parameters
----------
Hist : str
    The new history.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetVersionHist", R"doc(Set the version history of the waveform.

Parameters
----------
VerHist : list of tuple of str
    The new version history, as a list of pairs of strings, in the form
    returned by `Waveform.VersionHist`.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetT", R"doc(Set the times of the waveform.

The mode data and frame data are not changed, so the number of times should
normally equal `Waveform.NTimes`.  This method is identical to
`Waveform.SetTime`.

Parameters
----------
a : array_like of float
    The new times.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetTime", R"doc(Set the times of the waveform.

The mode data and frame data are not changed, so the number of times should
normally equal `Waveform.NTimes`.  This method is identical to
`Waveform.SetT`.

Parameters
----------
a : array_like of float
    The new times.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetFrame", R"doc(Set the frame data of the waveform.

The mode data are not changed; to rotate them along with the frame, use
`Waveform.RotateDecompositionBasis` instead.

Parameters
----------
a : array_like of Quaternion
    Rotors taking the standard basis into the basis with respect to which the
    modes are decomposed: none, a single rotor for a constant frame, or one
    rotor for each time step.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetFrameType", R"doc(Set the frame type of the waveform.

Only the recorded frame type changes; the data are not altered.

Parameters
----------
Type : int
    One of the module constants ``GWFrames.UnknownFrameType``,
    ``GWFrames.Inertial``, ``GWFrames.Coprecessing``, ``GWFrames.Coorbital``,
    and ``GWFrames.Corotating``.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetDataType", R"doc(Set the data type of the waveform.

Only the recorded data type changes; the data are not altered.

Parameters
----------
Type : int
    One of the module constants ``GWFrames.UnknownDataType``, ``GWFrames.h``,
    ``GWFrames.hdot``, ``GWFrames.Psi4``, ``GWFrames.Psi3``,
    ``GWFrames.Psi2``, ``GWFrames.Psi1``, and ``GWFrames.Psi0``.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetRIsScaledOut", R"doc(Record whether the radial dependence has been scaled out of the data.

When it has, the data represent quantities such as :math:`r h` or
:math:`r \Psi_4`.  Only the recorded flag changes; the data are not altered.

Parameters
----------
Scaled : bool
    True if the radial dependence has been scaled out.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetMIsScaledOut", R"doc(Record whether the mass dependence has been scaled out of the data.

When it has, the data represent quantities such as :math:`h/M` or
:math:`M \Psi_4`, and the times represent :math:`t/M`.  Only the recorded flag
changes; the data are not altered.

Parameters
----------
Scaled : bool
    True if the mass dependence has been scaled out.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetLM", R"doc(Set the :math:`(\ell, m)` values of the modes.

The mode data are not changed, so the number of rows should normally equal
`Waveform.NModes`.

Parameters
----------
a : array_like of int, shape (N, 2)
    The new :math:`(\ell, m)` values, one row for each mode.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::SetData", R"doc(Set the mode data of the waveform, either entirely or at one point.

The times and the :math:`(\ell, m)` data are not changed.  When a single value
is set, the indices are not checked, so they must lie within the current shape
of the data.

Parameters
----------
i_Mode : int, optional
    Index of the mode to set.  It must be given together with ``i_Time``.
i_Time : int, optional
    Index of the time step to set.
a : array_like of complex or complex
    If ``i_Mode`` and ``i_Time`` are omitted, the complete mode data, with
    shape (NModes, NTimes), which replace the current data.  Otherwise, the
    single value to store for that mode and time step.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::ResizeData", R"doc(Change the shape of the array of mode data.

If the shape changes, the existing data are discarded and replaced by zeros.
The times and the :math:`(\ell, m)` data are not changed.

Parameters
----------
NModes : int
    The new number of modes.
NTimes : int
    The new number of time steps.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::swap", R"doc(Exchange the entire contents of this waveform with those of another.

All of the data, including the histories, are exchanged without being copied.

Parameters
----------
b : Waveform
    The other waveform.)doc"},
    {"Waveform::NTimes", R"doc(Return the number of time steps.

Returns
-------
int
    The number of time steps.)doc"},
    {"Waveform::NModes", R"doc(Return the number of modes.

Returns
-------
int
    The number of modes.)doc"},
    {"Waveform::SpinWeight", R"doc(Return the spin weight of the data.

Returns
-------
int
    The spin weight.)doc"},
    {"Waveform::BoostWeight", R"doc(Return the boost weight of the data.

Returns
-------
int
    The boost weight.)doc"},
    {"Waveform::HistoryStr", R"doc(Return the history of the waveform.

The history records the operations applied to the waveform since it was
created.

Returns
-------
str
    The history.)doc"},
    {"Waveform::VersionHist", R"doc(Return the version history of the waveform.

Returns
-------
list of tuple of str
    The version history, as a list of pairs of strings.)doc"},
    {"Waveform::FrameType", R"doc(Return the frame type of the waveform.

Returns
-------
int
    One of the module constants ``GWFrames.UnknownFrameType``,
    ``GWFrames.Inertial``, ``GWFrames.Coprecessing``, ``GWFrames.Coorbital``,
    and ``GWFrames.Corotating``.)doc"},
    {"Waveform::DataType", R"doc(Return the data type of the waveform.

Returns
-------
int
    One of the module constants ``GWFrames.UnknownDataType``, ``GWFrames.h``,
    ``GWFrames.hdot``, ``GWFrames.Psi4``, ``GWFrames.Psi3``,
    ``GWFrames.Psi2``, ``GWFrames.Psi1``, and ``GWFrames.Psi0``.)doc"},
    {"Waveform::FrameTypeString", R"doc(Return the name of the frame type of the waveform.

Returns
-------
str
    The name of the frame type, such as ``"Inertial"`` or ``"Corotating"``.)doc"},
    {"Waveform::DataTypeString", R"doc(Return the name of the data type of the waveform.

Returns
-------
str
    The name of the data type, such as ``"h"`` or ``"Psi4"``.)doc"},
    {"Waveform::DataTypeLaTeXString", R"doc(Return a LaTeX expression for the data type of the waveform.

Returns
-------
str
    The LaTeX expression, such as ``"\dot{h}"`` or ``"\Psi_4"``.)doc"},
    {"Waveform::DescriptorString", R"doc(Return a short description of the data, suitable for a file name.

The description combines the name of the data type with the factors of radius
and mass that have been scaled out of the data, as recorded by
`Waveform.RIsScaledOut` and `Waveform.MIsScaledOut`; examples are
``"rhOverM"``, ``"rMPsi4"``, and ``"rhdot"``.

Returns
-------
str
    The description.)doc"},
    {"Waveform::RIsScaledOut", R"doc(Return whether the radial dependence has been scaled out of the data.

Returns
-------
bool
    True if the data represent quantities such as :math:`r h` or
    :math:`r \Psi_4`.)doc"},
    {"Waveform::MIsScaledOut", R"doc(Return whether the mass dependence has been scaled out of the data.

Returns
-------
bool
    True if the data represent quantities such as :math:`h/M` or
    :math:`M \Psi_4`, and the times represent :math:`t/M`.)doc"},
    {"Waveform::T", R"doc(Return the times of the waveform, or the time at one time step.

Parameters
----------
TimeIndex : int, optional
    Index of the time step.  If omitted, all of the times are returned.

Returns
-------
float or numpy.ndarray of float
    The time at ``TimeIndex``, or the array of all times.)doc"},
    {"Waveform::Frame", R"doc(Return the frame data of the waveform, or the frame at one time step.

The frame data consist of no rotors, of a single rotor for a constant frame, or
of one rotor for each time step.  When ``TimeIndex`` is given for a constant
frame, that single rotor is returned; there must be at least one rotor.

Parameters
----------
TimeIndex : int, optional
    Index of the time step.  If omitted, all of the frame data are returned.

Returns
-------
Quaternion or numpy.ndarray of Quaternion
    The rotor at ``TimeIndex``, or the array of all frame data.)doc"},
    {"Waveform::LM", R"doc(Return the :math:`(\ell, m)` values of the modes, or of one mode.

Parameters
----------
Mode : int, optional
    Index of the mode.  If omitted, the values for all modes are returned.

Returns
-------
numpy.ndarray of int, shape (2,) or (N, 2)
    The values :math:`(\ell, m)` of the mode ``Mode``, or one row of such
    values for each mode.)doc"},
    {"Waveform::Re", R"doc(Return the real parts of the mode data.

Parameters
----------
Mode : int, optional
    Index of the mode.  If omitted, the data for all modes are returned.
TimeIndex : int, optional
    Index of the time step.  If omitted, the data at all times are returned.
    It may be given only together with ``Mode``.

Returns
-------
float or numpy.ndarray of float
    The real part of the data of one mode at one time step; or, if
    ``TimeIndex`` is omitted, the real parts for that mode at every time; or,
    if both are omitted, an array of shape (NModes, NTimes).)doc"},
    {"Waveform::Im", R"doc(Return the imaginary parts of the mode data.

Parameters
----------
Mode : int, optional
    Index of the mode.  If omitted, the data for all modes are returned.
TimeIndex : int, optional
    Index of the time step.  If omitted, the data at all times are returned.
    It may be given only together with ``Mode``.

Returns
-------
float or numpy.ndarray of float
    The imaginary part of the data of one mode at one time step; or, if
    ``TimeIndex`` is omitted, the imaginary parts for that mode at every time;
    or, if both are omitted, an array of shape (NModes, NTimes).)doc"},
    {"Waveform::Abs", R"doc(Return the absolute values of the mode data.

Parameters
----------
Mode : int, optional
    Index of the mode.  If omitted, the data for all modes are returned.
TimeIndex : int, optional
    Index of the time step.  If omitted, the data at all times are returned.
    It may be given only together with ``Mode``.

Returns
-------
float or numpy.ndarray of float
    The absolute value of the data of one mode at one time step; or, if
    ``TimeIndex`` is omitted, the absolute values for that mode at every
    time; or, if both are omitted, an array of shape (NModes, NTimes).)doc"},
    {"Waveform::Arg", R"doc(Return the arguments (phases) of the mode data.

The arguments are not unwrapped; they lie between :math:`-\pi` and
:math:`\pi`.  For a continuous phase, use `Waveform.ArgUnwrapped`.

Parameters
----------
Mode : int, optional
    Index of the mode.  If omitted, the data for all modes are returned.
TimeIndex : int, optional
    Index of the time step.  If omitted, the data at all times are returned.
    It may be given only together with ``Mode``.

Returns
-------
float or numpy.ndarray of float
    The argument of the data of one mode at one time step; or, if
    ``TimeIndex`` is omitted, the arguments for that mode at every time; or,
    if both are omitted, an array of shape (NModes, NTimes).)doc"},
    {"Waveform::ArgUnwrapped", R"doc(Return the unwrapped arguments (phases) of the mode data.

Multiples of :math:`2\pi` are added to the arguments of each mode to remove
the jumps between :math:`-\pi` and :math:`\pi`, so that the phase is
continuous in time.

Parameters
----------
Mode : int, optional
    Index of the mode.  If omitted, the data for all modes are returned.

Returns
-------
numpy.ndarray of float
    The unwrapped arguments for the mode ``Mode`` at every time, or, if
    ``Mode`` is omitted, an array of shape (NModes, NTimes).)doc"},
    {"Waveform::Data", R"doc(Return the complex mode data.

Parameters
----------
Mode : int, optional
    Index of the mode.  If omitted, the data for all modes are returned.
TimeIndex : int, optional
    Index of the time step.  If omitted, the data at all times are returned.
    It may be given only together with ``Mode``.

Returns
-------
complex or numpy.ndarray of complex
    The data of one mode at one time step; or, if ``TimeIndex`` is omitted,
    the data for that mode at every time; or, if both are omitted, an array of
    shape (NModes, NTimes).)doc"},
    {"Waveform::EllMax", R"doc(Return the largest value of :math:`\ell` present in the data.

Returns
-------
int
    The largest value of :math:`\ell`.)doc"},
    {"Waveform::FindModeIndex", R"doc(Return the index of the mode with the given :math:`(\ell, m)` values.

If no such mode is present, an `IndexError` is raised.

Parameters
----------
L : int
    The value of :math:`\ell`.
M : int
    The value of :math:`m`.

Returns
-------
int
    The index of the mode.)doc"},
    {"Waveform::FindModeIndexWithoutError", R"doc(Return the index of the mode with given :math:`(\ell, m)`, without raising.

Unlike `Waveform.FindModeIndex`, this method raises no error when the mode is
missing; instead, it returns ``NModes() + 1``, which is not a valid index.

Parameters
----------
L : int
    The value of :math:`\ell`.
M : int
    The value of :math:`m`.

Returns
-------
int
    The index of the mode, or ``NModes() + 1`` if it is not present.)doc"},
    {"Waveform::DataDot", R"doc(Return the time derivative of the data of one mode.

The derivative is computed with a fourth-order finite-difference formula for
nonuniformly spaced times.  Only the mode data are differentiated; the time
dependence of a rotating frame is not taken into account.

Parameters
----------
Mode : int
    Index of the mode.

Returns
-------
numpy.ndarray of complex
    The derivative at each time step.)doc"},
    {"Waveform::Differentiate", R"doc(Differentiate the waveform with respect to time, in place.

Each mode is differentiated as in `Waveform.DataDot`.  The waveform must be in
the inertial frame; for any other frame type, a `NotImplementedError` is
raised.  The boost weight is reduced by one, and the data type changes from
``h`` to ``hdot``, or from ``hdot`` to ``Psi4``; any other data type becomes
``UnknownDataType``.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::Norm", R"doc(Return the norm of the waveform at each time.

The norm is the sum over all modes of the squared magnitudes of the data,
:math:`\sum_{\ell,m} \lvert f^{\ell,m} \rvert^2`; it is called a norm in
analogy with ``std::norm`` for complex numbers in C++, which is the square of
the absolute value.  Optionally, the square root is taken, giving the usual
:math:`L^2` norm of the waveform.

Parameters
----------
TakeSquareRoot : bool, default False
    If true, the square root is taken at each time.

Returns
-------
numpy.ndarray of float
    The norm at each time step.)doc"},
    {"Waveform::MaxNormIndex", R"doc(Return the index of the time step at which the norm is largest.

The search for the maximum of `Waveform.Norm` begins at the time step with
index ``NTimes() // SkipFraction``, so that junk radiation near the beginning
of the data is skipped.  With ``SkipFraction = NTimes() + 1``, all of the data
are searched.

Parameters
----------
SkipFraction : int, default 4
    Inverse of the fraction of the data to skip before searching.

Returns
-------
int
    The index of the time step with the largest norm.)doc"},
    {"Waveform::MaxNormTime", R"doc(Return the time at which the norm is largest.

This is the time at the index returned by `Waveform.MaxNormIndex`, whose
documentation describes ``SkipFraction``.

Parameters
----------
SkipFraction : int, default 4
    Inverse of the fraction of the data to skip before searching.

Returns
-------
float
    The time with the largest norm.)doc"},
    {"Waveform::NormalizedAntisymmetry", R"doc(Return the normalized antisymmetry of the waveform at each time.

This is the quantity :math:`a` defined by Boyle et al. (2014): the difference
between the field at a point and the complex conjugate of the field at the
antipodal point, squared and integrated over the sphere, normalized by the
`Waveform.Norm` of the waveform (its total power at each instant), and with the
square root taken.  In terms of the modes,

.. math::

   a = \sqrt{\frac{\sum_{\ell,m} \lvert f^{\ell,m} - (-1)^{\ell+m}
   \bar{f}^{\ell,-m} \rvert^2}{4 \sum_{\ell,m} \lvert f^{\ell,m} \rvert^2}}.

Parameters
----------
LModesForAsymmetry : array_like of int, optional
    Values of :math:`\ell` to include in the numerator.  If empty or omitted,
    all :math:`\ell` values are used.  All modes are always used in the
    denominator.

Returns
-------
numpy.ndarray of float
    The normalized antisymmetry at each time step.)doc"},
    {"Waveform::DipoleMoment", R"doc(Return the dipole moment of the magnitude of the waveform at each time.

The dipole moment is
:math:`\vec{d} = \int \hat{n} \lvert f \rvert^2 \, d\Omega`.  Up to a
geometric factor, the dipole moment of :math:`\dot{h}` is the rate at which
momentum is emitted in gravitational waves.

Parameters
----------
ellMax : int, default 0
    Largest value of :math:`\ell` to include.  If it is 0, all modes up to
    `Waveform.EllMax` are used.

Returns
-------
numpy.ndarray of float, shape (N, 3)
    The dipole moment, one row for each time step.)doc"},
    {"Waveform::MinimalParityViolation", R"doc(Return the z-parity violation at each time, minimized over orientations.

At each time step, the decomposition basis is rotated so as to minimize the
violation of invariance under z parity (reflection across the x-y plane), as
measured by `Waveform.ZParityViolationNormalized`, and the minimum value is
returned.  The minimization is started from three different orientations, and
the smallest result is used.  Nonprecessing systems in a suitable frame should
have zero violation; precessing systems in any frame, and nonprecessing systems
in the wrong frame, show violations.

This method is slow, and it prints its progress to standard error.

Returns
-------
numpy.ndarray of float
    The minimal normalized violation at each time step.)doc"},
    {"Waveform::XParityInvolution", R"doc(Return the image of the waveform under x parity.

The x-parity operation reflects the x axis, :math:`x \to -x`, and takes the
complex conjugate.  In terms of the modes, it maps :math:`f^{\ell,m}` to
:math:`(-1)^m \bar{f}^{\ell,m}`.  The frame rotors are transformed with
`Quaternions.XParityConjugateSpinor`.  Every mode up to `Waveform.EllMax` must
be present.  This waveform is not changed.

Returns
-------
Waveform
    The transformed waveform.)doc"},
    {"Waveform::XParitySymmetricPart", R"doc(Return the part of the waveform that is symmetric under x parity.

This is :math:`(f + P_x f)/2`, where :math:`P_x f` is the image of the waveform
under x parity, as returned by `Waveform.XParityInvolution`.  The same
combination is applied to the frame rotors, which therefore are generally not
unit quaternions.  This waveform is not changed.

Returns
-------
Waveform
    The symmetric part.)doc"},
    {"Waveform::XParityAntisymmetricPart", R"doc(Return the part of the waveform that is antisymmetric under x parity.

This is :math:`(f - P_x f)/2`, where :math:`P_x f` is the image of the waveform
under x parity, as returned by `Waveform.XParityInvolution`.  The same
combination is applied to the frame rotors, which therefore are generally not
unit quaternions.  This waveform is not changed.

Returns
-------
Waveform
    The antisymmetric part.)doc"},
    {"Waveform::XParityViolationSquared", R"doc(Return the squared violation of x-parity symmetry at each time.

This is the sum over modes of the squared magnitude of the antisymmetric part
(see `Waveform.XParityAntisymmetricPart`),

.. math::

   \sum_{\ell,m} \frac{1}{4} \lvert f^{\ell,m} - (-1)^m \bar{f}^{\ell,m}
   \rvert^2.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sum.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float
    The squared violation at each time step.)doc"},
    {"Waveform::XParityViolationNormalized", R"doc(Return the normalized violation of x-parity symmetry at each time.

This is the square root of the ratio of `Waveform.XParityViolationSquared` to
`Waveform.Norm`.  The norm always includes all modes.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the violation.  If empty or omitted,
    all :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float
    The normalized violation at each time step.)doc"},
    {"Waveform::YParityInvolution", R"doc(Return the image of the waveform under y parity.

The y-parity operation reflects the y axis, :math:`y \to -y`, and takes the
complex conjugate.  In terms of the modes, it maps :math:`f^{\ell,m}` to
:math:`\bar{f}^{\ell,m}`.  The frame rotors are transformed with
`Quaternions.YParityConjugateSpinor`.  Every mode up to `Waveform.EllMax` must
be present.  This waveform is not changed.

Returns
-------
Waveform
    The transformed waveform.)doc"},
    {"Waveform::YParitySymmetricPart", R"doc(Return the part of the waveform that is symmetric under y parity.

This is :math:`(f + P_y f)/2`, where :math:`P_y f` is the image of the waveform
under y parity, as returned by `Waveform.YParityInvolution`.  The same
combination is applied to the frame rotors, which therefore are generally not
unit quaternions.  This waveform is not changed.

Returns
-------
Waveform
    The symmetric part.)doc"},
    {"Waveform::YParityAntisymmetricPart", R"doc(Return the part of the waveform that is antisymmetric under y parity.

This is :math:`(f - P_y f)/2`, where :math:`P_y f` is the image of the waveform
under y parity, as returned by `Waveform.YParityInvolution`.  The same
combination is applied to the frame rotors, which therefore are generally not
unit quaternions.  This waveform is not changed.

Returns
-------
Waveform
    The antisymmetric part.)doc"},
    {"Waveform::YParityViolationSquared", R"doc(Return the squared violation of y-parity symmetry at each time.

This is the sum over modes of the squared magnitude of the antisymmetric part
(see `Waveform.YParityAntisymmetricPart`),

.. math::

   \sum_{\ell,m} \frac{1}{4} \lvert f^{\ell,m} - \bar{f}^{\ell,m} \rvert^2.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sum.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float
    The squared violation at each time step.)doc"},
    {"Waveform::YParityViolationNormalized", R"doc(Return the normalized violation of y-parity symmetry at each time.

This is the square root of the ratio of `Waveform.YParityViolationSquared` to
`Waveform.Norm`.  The norm always includes all modes.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the violation.  If empty or omitted,
    all :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float
    The normalized violation at each time step.)doc"},
    {"Waveform::ZParityInvolution", R"doc(Return the image of the waveform under z parity.

The z-parity operation reflects the z axis, :math:`z \to -z`, and takes the
complex conjugate.  In terms of the modes, it maps :math:`f^{\ell,m}` to
:math:`(-1)^\ell \bar{f}^{\ell,-m}`.  The frame rotors are transformed with
`Quaternions.ZParityConjugateSpinor`.  Every mode up to `Waveform.EllMax` must
be present.  This waveform is not changed.

Returns
-------
Waveform
    The transformed waveform.)doc"},
    {"Waveform::ZParitySymmetricPart", R"doc(Return the part of the waveform that is symmetric under z parity.

This is :math:`(f + P_z f)/2`, where :math:`P_z f` is the image of the waveform
under z parity, as returned by `Waveform.ZParityInvolution`.  The same
combination is applied to the frame rotors, which therefore are generally not
unit quaternions.  This waveform is not changed.

Returns
-------
Waveform
    The symmetric part.)doc"},
    {"Waveform::ZParityAntisymmetricPart", R"doc(Return the part of the waveform that is antisymmetric under z parity.

This is :math:`(f - P_z f)/2`, where :math:`P_z f` is the image of the waveform
under z parity, as returned by `Waveform.ZParityInvolution`.  The same
combination is applied to the frame rotors, which therefore are generally not
unit quaternions.  This waveform is not changed.

Returns
-------
Waveform
    The antisymmetric part.)doc"},
    {"Waveform::ZParityViolationSquared", R"doc(Return the squared violation of z-parity symmetry at each time.

This is the sum over modes of the squared magnitude of the antisymmetric part
(see `Waveform.ZParityAntisymmetricPart`),

.. math::

   \sum_{\ell,m} \frac{1}{4} \lvert f^{\ell,m} - (-1)^\ell \bar{f}^{\ell,-m}
   \rvert^2.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sum.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float
    The squared violation at each time step.)doc"},
    {"Waveform::ZParityViolationNormalized", R"doc(Return the normalized violation of z-parity symmetry at each time.

This is the square root of the ratio of `Waveform.ZParityViolationSquared` to
`Waveform.Norm`.  The norm always includes all modes.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the violation.  If empty or omitted,
    all :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float
    The normalized violation at each time step.)doc"},
    {"Waveform::ParityInvolution", R"doc(Return the image of the waveform under parity.

The parity operation inverts all three axes, :math:`\vec{x} \to -\vec{x}`, and
takes the complex conjugate.  In terms of the modes, it maps
:math:`f^{\ell,m}` to :math:`(-1)^{\ell+m} \bar{f}^{\ell,-m}`.  The frame
rotors are transformed with `Quaternions.ParityConjugateSpinor`, which leaves
them unchanged.  Every mode up to `Waveform.EllMax` must be present.  This
waveform is not changed.  The result is the same as that of
`Waveform.ConjugateAntipodalEvaluation`.

Returns
-------
Waveform
    The transformed waveform.)doc"},
    {"Waveform::ParitySymmetricPart", R"doc(Return the part of the waveform that is symmetric under parity.

This is :math:`(f + P f)/2`, where :math:`P f` is the image of the waveform
under parity, as returned by `Waveform.ParityInvolution`.  The same combination
is applied to the frame rotors, which are therefore unchanged.  This waveform
is not changed.

Returns
-------
Waveform
    The symmetric part.)doc"},
    {"Waveform::ParityAntisymmetricPart", R"doc(Return the part of the waveform that is antisymmetric under parity.

This is :math:`(f - P f)/2`, where :math:`P f` is the image of the waveform
under parity, as returned by `Waveform.ParityInvolution`.  The same combination
is applied to the frame rotors, which are therefore all zero.  This waveform is
not changed.

Returns
-------
Waveform
    The antisymmetric part.)doc"},
    {"Waveform::ParityViolationSquared", R"doc(Return the squared violation of parity symmetry at each time.

This is the sum over modes of the squared magnitude of the antisymmetric part
(see `Waveform.ParityAntisymmetricPart`),

.. math::

   \sum_{\ell,m} \frac{1}{4} \lvert f^{\ell,m} - (-1)^{\ell+m}
   \bar{f}^{\ell,-m} \rvert^2.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sum.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float
    The squared violation at each time step.)doc"},
    {"Waveform::ParityViolationNormalized", R"doc(Return the normalized violation of parity symmetry at each time.

This is the square root of the ratio of `Waveform.ParityViolationSquared` to
`Waveform.Norm`.  The norm always includes all modes.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the violation.  If empty or omitted,
    all :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float
    The normalized violation at each time step.)doc"},
    {"Waveform::ConjugateAntipodalEvaluation", R"doc(Return the complex conjugate of the waveform at the antipodal point.

The result is the field :math:`\bar{f}(-\hat{n})`.  In terms of the modes, it
maps :math:`f^{\ell,m}` to :math:`(-1)^{\ell+m} \bar{f}^{\ell,-m}`.  This is
the same as `Waveform.ParityInvolution`.  This waveform is not changed.

Returns
-------
Waveform
    The transformed waveform.)doc"},
    {"Waveform::ConjugateAntipodalSymmetricPart", R"doc(Return the conjugate-antipodal symmetric part of the waveform.

This is :math:`[f(\hat{n}) + \bar{f}(-\hat{n})]/2`, where
:math:`\bar{f}(-\hat{n})` is returned by
`Waveform.ConjugateAntipodalEvaluation`.  This is the same as
`Waveform.ParitySymmetricPart`.  This waveform is not changed.

Returns
-------
Waveform
    The symmetric part.)doc"},
    {"Waveform::ConjugateAntipodalAntisymmetricPart", R"doc(Return the conjugate-antipodal antisymmetric part of the waveform.

This is :math:`[f(\hat{n}) - \bar{f}(-\hat{n})]/2`, where
:math:`\bar{f}(-\hat{n})` is returned by
`Waveform.ConjugateAntipodalEvaluation`.  This is the same as
`Waveform.ParityAntisymmetricPart`.  This waveform is not changed.

Returns
-------
Waveform
    The antisymmetric part.)doc"},
    {"Waveform::RotatePhysicalSystem", R"doc(Rotate the physical content of the waveform, in place.

The physical system is rotated by ``R_phys``, while the basis with respect to
which the modes are decomposed is left in place.  The modes are transformed
accordingly, and the frame data are multiplied on the right by
:math:`\bar{R}_{\mathrm{phys}}`, so that the frame data record the total
rotation.  The frame type is not changed.

Parameters
----------
R_phys : Quaternion or array_like of Quaternion
    Rotor by which to rotate the system, or one rotor for each time step.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::RotateDecompositionBasis", R"doc(Rotate the basis in which the modes are decomposed, in place.

The basis is rotated by ``R_frame``, while the physical system is left in
place.  The modes are transformed accordingly, and the frame data are
multiplied on the right by ``R_frame``, so that they record the total rotation
of the basis.  The frame type is not changed.

Parameters
----------
R_frame : Quaternion or array_like of Quaternion
    Rotor by which to rotate the basis, or one rotor for each time step.  If
    an empty array is given, nothing is done.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::LdtVector", R"doc(Return the vector :math:`\langle L \partial_t \rangle` at each time.

This quantity, defined in "Angular velocity of gravitational radiation from
precessing binaries and the corotating frame"
(https://arxiv.org/abs/1302.2919), is

.. math::

   \langle L \partial_t \rangle^a = \sum_{\ell,m,m'} \Im \left[
   \bar{f}^{\ell,m'} \langle \ell,m' | L_a | \ell,m \rangle \dot{f}^{\ell,m}
   \right].

The vector is given with respect to the (possibly rotating) mode frame
:math:`(X, Y, Z)`, rather than the inertial frame :math:`(x, y, z)`.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sum.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float, shape (N, 3)
    The vector, one row for each time step.)doc"},
    {"Waveform::LLMatrix", R"doc(Return the matrix :math:`\langle LL \rangle` at each time.

This quantity, defined in "Angular velocity of gravitational radiation from
precessing binaries and the corotating frame"
(https://arxiv.org/abs/1302.2919), is

.. math::

   \langle LL \rangle^{ab} = \sum_{\ell,m,m'} \bar{f}^{\ell,m'}
   \langle \ell,m' | L_a L_b | \ell,m \rangle f^{\ell,m}.

The matrix is given with respect to the (possibly rotating) mode frame
:math:`(X, Y, Z)`, rather than the inertial frame :math:`(x, y, z)`.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sum.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
list of Matrix
    The 3-by-3 matrix at each time step.)doc"},
    {"Waveform::LLDominantEigenvector", R"doc(Return the dominant eigenvector of :math:`\langle LL \rangle` at each time.

This is the principal axis of `Waveform.LLMatrix`, as prescribed by
O'Shaughnessy et al. (2011).  The eigenvectors are normalized, and their signs
are chosen so that they vary continuously in time, with the first one more
nearly parallel than antiparallel to ``RoughInitialEllDirection``.  For a
post-Newtonian system, the orbital angular-velocity direction is a good guess
for that direction; the default, the z axis, is usually correct.

The vectors are given with respect to the (possibly rotating) mode frame
:math:`(X, Y, Z)`, rather than the inertial frame :math:`(x, y, z)`.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sum.  If empty or omitted, all
    :math:`\ell` values are used.
RoughInitialEllDirection : Quaternion, default Quaternions.zHat
    Rough guess at the direction of the eigenvector at the first time step.

Returns
-------
numpy.ndarray of float, shape (N, 3)
    The unit eigenvector, one row for each time step.)doc"},
    {"Waveform::AngularVelocityVector", R"doc(Return the angular velocity of the waveform at each time.

This is the angular velocity defined in Sec. II of "Angular velocity of
gravitational radiation from precessing binaries and the corotating frame"
(https://arxiv.org/abs/1302.2919), found by solving
:math:`\langle LL \rangle \vec{\omega} = -\langle L \partial_t \rangle` at
each time step (see `Waveform.LLMatrix` and `Waveform.LdtVector`).  The vector
is given with respect to the (possibly rotating) mode frame :math:`(X, Y, Z)`,
rather than the inertial frame :math:`(x, y, z)`; for the latter, use
`Waveform.AngularVelocityVectorRelativeToInertial`.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sums.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float, shape (N, 3)
    The angular velocity, one row for each time step.)doc"},
    {"Waveform::AngularVelocityVectorRelativeToInertial", R"doc(Return the angular velocity of the waveform relative to the inertial frame.

This is the angular velocity defined in Sec. II of "Angular velocity of
gravitational radiation from precessing binaries and the corotating frame"
(https://arxiv.org/abs/1302.2919).  Unlike most methods of `Waveform`, this one
gives the vector with respect to the inertial frame :math:`(x, y, z)`, rather
than the (possibly rotating) mode frame :math:`(X, Y, Z)`: the result of
`Waveform.AngularVelocityVector` is rotated by the frame, and the angular
velocity of a time-dependent frame is added.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sums.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
numpy.ndarray of float, shape (N, 3)
    The angular velocity, one row for each time step.)doc"},
    {"Waveform::CorotatingFrame", R"doc(Return the corotating frame of the waveform.

The corotating frame is the frame whose angular velocity is the angular
velocity of the waveform, as returned by `Waveform.AngularVelocityVector`; it
is found with `Quaternions.FrameFromAngularVelocity`.  The rotors are given
relative to the (possibly rotating) mode frame :math:`(X, Y, Z)`, rather than
the inertial frame :math:`(x, y, z)`.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sums.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
numpy.ndarray of Quaternion
    The rotor of the corotating frame at each time step.)doc"},
    {"Waveform::TransformToCoprecessingFrame", R"doc(Transform the waveform to the coprecessing frame, in place.

In the coprecessing frame, the z axis is aligned with the dominant eigenvector
of :math:`\langle LL \rangle` (the vector :math:`V_h` of O'Shaughnessy et al.
2011; see `Waveform.LLDominantEigenvector`), and the rotation about the z axis
satisfies the minimal-rotation condition of Boyle, Owen, and Pfeiffer (2011).
The frame type is set to ``Coprecessing``.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sum.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::TransformToAngularVelocityFrame", R"doc(Transform the waveform to the frame aligned with its angular velocity.

In this frame, the z axis is aligned with the angular-velocity vector of
`Waveform.AngularVelocityVector`, and the rotation about the z axis satisfies
the minimal-rotation condition.  Note that this is not the corotating frame.
The frame type is set to ``Coprecessing``.  The transformation is done in
place.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sums.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::TransformToCorotatingFrame", R"doc(Transform the waveform to the corotating frame, in place.

The decomposition basis is rotated by the rotors of `Waveform.CorotatingFrame`,
and the frame type is set to ``Corotating``.  This leaves a constant of
integration unset; it may be fixed by aligning the frame to the modes, for
example with `Waveform.AlignDecompositionFrameToModes`.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to include in the sums.  If empty or omitted, all
    :math:`\ell` values are used.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::TransformToInertialFrame", R"doc(Transform the waveform to the inertial frame, in place.

The stored frame data are used to transform from whatever rotating frame the
waveform is in to a stationary, inertial frame.  This is the usual frame at
future null infinity, in which gravitational-wave observations are made.  After
the transformation, the frame data are empty and the frame type is
``Inertial``.

If the frame type is already ``Inertial`` or is ``UnknownFrameType``, or if
there are no frame data, a warning is printed and nothing is done.  Otherwise,
there must be one rotor for each time step.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::GetAlignmentsOfDecompositionFrameToModes", R"doc(Return the rotors aligning the decomposition frame to the modes at each time.

At each time step, this finds the rotor :math:`R_{\varepsilon}` by which the
decomposition basis would be rotated to align it with the waveform, as in
`Waveform.AlignDecompositionFrameToModes`.  Computing the rotors at every time
step at once is useful when many such alignments must be tried, because
interpolation is slow.  Unlike similar methods, this one does not use an
approximate direction of :math:`\hat{n}` to choose the direction of the X axis.

The waveform must be in the coprecessing, coorbital, or corotating frame, with
one frame rotor for each time step.

Parameters
----------
Lmodes : array_like of int, optional
    Values of :math:`\ell` to use in computing :math:`\langle LL \rangle`.  If
    empty or omitted, all :math:`\ell` values are used.

Returns
-------
numpy.ndarray of Quaternion
    The rotor :math:`R_{\varepsilon}` at each time step.)doc"},
    {"Waveform::GetAlignmentOfDecompositionFrameToModes", R"doc(Return the rotor aligning the decomposition frame to the modes at one time.

This finds the rotor :math:`R_{\varepsilon}` that
`Waveform.AlignDecompositionFrameToModes` applies, without applying it; see
that method for details.  The waveform must be in the coprecessing, coorbital,
or corotating frame, and ``t_fid`` must lie within its times.

Parameters
----------
t_fid : float
    Fiducial time at which the alignment is done.
nHat_t_fid : Quaternion, default Quaternions.xHat
    Approximate direction of :math:`\hat{n}` (the vector from one black hole
    to the other) at ``t_fid``, which sets the direction of the X axis.
Lmodes : array_like of int, optional
    Values of :math:`\ell` to use in computing :math:`\langle LL \rangle`.  If
    empty or omitted, all :math:`\ell` values are used.

Returns
-------
Quaternion
    The rotor :math:`R_{\varepsilon}`.)doc"},
    {"Waveform::AlignDecompositionFrameToModes", R"doc(Fix the orientation of the corotating frame by aligning it to the modes.

The corotating frame is defined only up to a constant rotor
:math:`R_{\varepsilon}`: if :math:`R_{\mathrm{corot}}` is corotating, then so
is :math:`R_{\mathrm{corot}} R_{\varepsilon}`.  This method uses that freedom
to align the frame with the modes at the fiducial time ``t_fid``.  The Z axis
of the decomposition basis is placed along the dominant eigenvector of
:math:`\langle LL \rangle` (as suggested by O'Shaughnessy et al.), with its
sign chosen to be closer to the angular velocity than not; the rotation about
that axis is fixed by the phases of the :math:`(2, \pm 2)` modes; and the
remaining choice of sign of the X axis is made so that X lies closer to
``nHat_t_fid`` than not.  The decomposition basis is then rotated by
:math:`R_{\varepsilon}`.

The waveform must be in the coprecessing, coorbital, or corotating frame, and
``t_fid`` must lie within its times.

Parameters
----------
t_fid : float
    Fiducial time at which the alignment is done.
nHat_t_fid : Quaternion, default Quaternions.xHat
    Approximate direction of :math:`\hat{n}` (the vector from one black hole
    to the other) at ``t_fid``.
Lmodes : array_like of int, optional
    Values of :math:`\ell` to use in computing :math:`\langle LL \rangle`.  If
    empty or omitted, all :math:`\ell` values are used.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::Compare", R"doc(Return the difference between another waveform and this one.

The mode data of the result are the data of ``B`` minus the data of this
waveform, both interpolated with cubic splines to the times common to the two
(as found by `Intersection`).  The frame of the result is the rotation taking
the frame of this waveform into the frame of ``B``.  The data are compared as
complex numbers, rather than as moduli and phases.  A warning is printed if the
two waveforms have different frame types.

Parameters
----------
B : Waveform
    The waveform to compare with this one.
MinTimeStep : float, default 0.005
    Intended lower limit on the time step of the result.  It is currently
    ignored; the default of `Intersection`, 0.005, is always used.
MinTime : float, default -3.0e300
    Intended earliest time of the result.  It is currently ignored.

Returns
-------
Waveform
    The difference.)doc"},
    {"Waveform::Hybridize", R"doc(Return a hybrid of this waveform and another.

The result has exactly the data of this waveform before ``t1`` and exactly the
data of ``B`` after ``t2``, with a smooth (:math:`C^\infty`) blend between the
two in the interval from ``t1`` to ``t2``; the frame data are blended in the
same way.  The waveforms are not aligned first; that must already have been
done (for example, with `AlignWaveforms`).  Both waveforms must have data
throughout the transition interval.  This waveform is not changed.

Parameters
----------
B : Waveform
    The waveform with which to hybridize this one.
t1 : float
    Beginning of the transition interval.
t2 : float
    End of the transition interval.
tMinStep : float, default 0.005
    Lower limit on the time step of the result.

Returns
-------
Waveform
    The hybrid waveform.)doc"},
    {"Waveform::EvaluateAtPoint", R"doc(Evaluate the waveform at a point on the sphere.

The mode data are summed with the spin-weighted spherical harmonics at the
given point, taking the frame data into account, so that the angles are
measured with respect to the inertial frame.  A warning is printed if the
frame type is ``UnknownFrameType``.

Parameters
----------
vartheta : float
    Polar angle of the point.
varphi : float
    Azimuthal angle of the point.
i_0 : int, default 0
    Index of the first time step at which to evaluate.
i_1 : int, default -1
    Index one past the last time step at which to evaluate.  If it is -1,
    the waveform is evaluated through the last time step.

Returns
-------
numpy.ndarray of complex
    The value of the waveform at each time step from ``i_0`` up to, but not
    including, ``i_1``.)doc"},
    {"Waveform::InterpolateToPoint", R"doc(Evaluate the waveform at a point on the sphere and an instant of time.

The waveform is evaluated at the point, as in `Waveform.EvaluateAtPoint`, at
the four time steps nearest ``t_i``, and the result is interpolated in time
with a cubic spline.  The waveform must have at least four time steps, and
``t_i`` must lie within its times.

Parameters
----------
vartheta : float
    Polar angle of the point.
varphi : float
    Azimuthal angle of the point.
t_i : float
    Time at which to evaluate.

Returns
-------
complex
    The value of the waveform.)doc"},
    {"Waveform::Translate", R"doc(Return the waveform translated in space by a time-dependent displacement.

At each time :math:`t` and direction :math:`\hat{n}`, the field of the result
is the field of this waveform at the time :math:`t - \hat{n} \cdot \delta
\vec{x}(t)`.  The result is decomposed into modes again, and it contains only
the times at which data are available in every direction.  The displacement is
given relative to the inertial frame.  This waveform may be in any frame, and
need not be transformed first; the result is always in the inertial frame,
without frame data.  This waveform is not changed.

This method is slow, because it interpolates the waveform to many points at
each time step.

Parameters
----------
deltax : array_like of float, shape (N, 3)
    The displacement at each time step.

Returns
-------
Waveform
    The translated waveform.)doc"},
    {"Waveform::BoostPsi4", R"doc(Apply a boost to :math:`\Psi_4` data, in place.

This method does three things at each time step.  First, it evaluates the
waveform on what will become an equiangular grid after transformation by the
boost.  Second, at each point of that grid, it combines the present value of
:math:`\Psi_4` and its conjugate to give the value of :math:`\Psi_4` observed
in the boosted frame.  Finally, it decomposes the result into modes using the
new equiangular grid.  Time steps with negligible velocity are left unchanged.

Parameters
----------
v : array_like of float, shape (N, 3)
    The three-velocity of the boosted frame relative to the present frame, at
    each time step.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::BoostHFaked", R"doc(Apply an approximate boost to :math:`h` data, in place.

This method does three things at each time step.  First, it evaluates the
waveform on what will become an equiangular grid after transformation by the
boost.  Second, at each point of that grid, it combines the present value of
:math:`h` and its conjugate to give the value of :math:`h` observed in the
boosted frame.  Finally, it decomposes the result into modes using the new
equiangular grid.  Time steps with negligible velocity are left unchanged.

The transformation is that of :math:`\Psi_4`, multiplied at each point by the
appropriate power (-2) of the Lorentz factor.  This assumes that
:math:`\ddot{h} = \Psi_4` in both frames, which need not be true; hence
"Faked" in the name of this method, which prints a warning to that effect.

Parameters
----------
v : array_like of float, shape (N, 3)
    The three-velocity of the boosted frame relative to the present frame, at
    each time step.

Returns
-------
Waveform
    This object, modified in place.)doc"},
    {"Waveform::Output", R"doc(Write the waveform to a text file.

The file begins with the history of the waveform and a header naming the
columns.  Each following line holds a time, followed by the real and imaginary
parts of each mode.  Such a file can be read by the constructor
``Waveform(FileName, "ReIm")`` if its modes are in the order that the
constructor assumes.

Parameters
----------
FileName : str
    Path of the file to write.
precision : int, default 14
    Number of significant digits written.

Returns
-------
Waveform
    This object, unchanged.)doc"},
    {"AlignWaveforms", R"doc(Align two waveforms in time and orientation, in place.

If either waveform is in the inertial frame, it is first transformed to its
corotating frame; otherwise, it must already be in the corotating frame (the
coorbital frame, for example, raises an error).  The frame of each waveform is
then aligned to its modes at the time ``t_mid = (t_1 + t_2) / 2``, as in
`Waveform.AlignDecompositionFrameToModes`.  Finally, ``B`` is shifted in time
and rotated by a constant rotor so that it matches ``A`` as closely as possible
over the interval from ``t_1`` to ``t_2``; during this optimization, the
alignment of the frame of ``B`` to its modes is adjusted to account for the
changing meaning of ``t_mid``.

The times ``t_1`` and ``t_2`` refer to the time axis of ``A``.  The algorithm
assumes that the waveforms are already reasonably well aligned in time.  In
particular, the final value of ``t_mid + deltat`` for ``B`` must lie somewhere
in the interval from ``t_1`` to ``t_2``, and after the time shift, ``B`` must
have data over all of that interval.  The initial, coarse search considers
time shifts of at most ``(t_2 - t_1) / 2`` in either direction, further
limited so that the shifted interval stays at least one time step inside the
data of ``B``.

As long as this last condition is satisfied, and the waveforms are even
remotely well sampled, and assuming the code's logic does not fail due to
otherwise weird or unexpected inputs, I (Mike Boyle) do hereby guarantee that
this algorithm will find the optimal alignment in both time and attitude.  Or
your money back.

Parameters
----------
A : Waveform
    The fixed waveform.  Its frame is aligned to its modes, but it is not
    shifted in time.
B : Waveform
    The adjusted waveform.  Its frame is aligned to its modes, and it is
    shifted in time and rotated.
t_1 : float
    Beginning of the alignment interval.
t_2 : float
    End of the alignment interval.
InitialEvaluations : int, default 0
    Number of time offsets at which the objective function is evaluated in
    the initial, coarse optimization.  If it is 0, or is not smaller than the
    number of available time steps of ``B``, every one of those time steps is
    used.
nHat_A : array_like of float, shape (3,), optional
    Approximate direction of :math:`\hat{n}` (the vector from black hole A to
    black hole B) in ``A`` at ``t_mid``, which sets the direction of the x axis
    of the rotating frame.  If empty or omitted, the x axis is used.
Debug : bool, default False
    If true, intermediate values of the objective function are written to
    the files ``XiIntegral.dat`` and ``UpsilonIntegral.dat`` in the current
    directory.)doc"},
    {"PNWaveform", R"doc(Post-Newtonian waveform of a precessing binary black hole.

A `PNWaveform` is a `Waveform` computed by evolving the orbital dynamics and
spins of a quasicircular binary in post-Newtonian (PN) theory, and then
evaluating the PN expressions for the modes of :math:`h`.  In addition to the
data of a `Waveform`, it stores at each time step the two spin vectors
``chi1`` and ``chi2``, the orbital angular velocity ``Omega_orb``, the
precessional angular velocity ``Omega_prec``, the orbital angular momentum
``L``, and the orbital phase ``Phi_orb``, which are returned by the methods of
the same names.  All of these vectors are expressed in the inertial basis.
Methods also return their magnitudes (`PNWaveform.Omega_orbMag`, and so on)
and their directions (`PNWaveform.OmegaHat_orb`, and so on).

A `PNWaveform` may be constructed in any of these ways:

``PNWaveform()``
    An empty waveform.
``PNWaveform(W)``
    A copy of the `PNWaveform` ``W``, including its history.
``PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i, ...)``
    The PN waveform of a binary with the given parameters, which are described
    below.  The optional parameters ``Omega_orb_0``, ``R_frame_i``,
    ``MinStepsPerOrbit``, ``PNWaveformModeOrder``, and
    ``PNOrbitalEvolutionOrder`` follow, in that order.

The PN system is defined with respect to an inertial basis
:math:`(\hat{x}, \hat{y}, \hat{z})`, in which the spin vectors are given.  A
second basis :math:`(\hat{X}, \hat{Y}, \hat{Z})` is obtained by rotating the
inertial basis by ``R_frame_i``.  At the initial time, black hole 1 lies on
the positive :math:`X` axis, black hole 2 lies on the negative :math:`X` axis,
and the orbital angular velocity points along the positive :math:`Z` axis
with magnitude ``Omega_orb_i``.  The orbital evolution is integrated forward
in time from there until the PN equations stop making progress near merger
(when the PN parameter :math:`v` reaches 1, stops increasing, or requires
very small time steps), and, if ``Omega_orb_0`` is smaller than
``Omega_orb_i``, also backward in time until the orbital frequency falls to
``Omega_orb_0``.  The time :math:`t = 0` corresponds to the moment at which
the orbital frequency is ``Omega_orb_i``.  The time steps are chosen by the
adaptive integrator, so they are not uniform.

The parameters of the third constructor are these:

``Approximant`` : str
    The PN approximant used to evolve the orbital frequency: ``"TaylorT1"``,
    ``"TaylorT4"``, or ``"TaylorT5"``.  Any other value raises an error.
``delta`` : float
    The normalized mass difference :math:`(M_1 - M_2) / (M_1 + M_2)`, which
    lies between -1 and 1 and is 0 for equal masses.  The total mass
    :math:`M_1 + M_2` is 1, so that all quantities are in units of the total
    mass.
``chi1_i``, ``chi2_i`` : array_like of float, shape (3,)
    The initial dimensionless spin vectors :math:`\vec{S}_1 / M_1^2` and
    :math:`\vec{S}_2 / M_2^2` of the two black holes, with components in the
    inertial basis.
``Omega_orb_i`` : float
    The orbital angular frequency at :math:`t = 0`, in units of :math:`1/M`.
``Omega_orb_0`` : float, default -1.0
    The orbital angular frequency at which the waveform begins.  If it is
    smaller than ``Omega_orb_i``, the system is also evolved backward in time
    to that frequency.  A value that is not positive (the default), or that
    is not smaller than ``Omega_orb_i``, makes the waveform begin at
    ``Omega_orb_i``.
``R_frame_i`` : Quaternion, default Quaternion(1, 0, 0, 0)
    The rotor taking the inertial basis to the initial basis of the binary,
    as described above.  The default places the binary in the standard
    position, with its orbital angular velocity along :math:`\hat{z}`.
``MinStepsPerOrbit`` : int, default 32
    The minimum number of time steps per orbit.  The step size of the
    integrator is limited to a fraction of the instantaneous orbital period,
    so that at least this many steps are taken in each orbit.
``PNWaveformModeOrder`` : float, default 3.5
    The PN order at which the waveform modes are computed, a multiple of 0.5
    from 0 to 3.5.
``PNOrbitalEvolutionOrder`` : float, default 4.0
    The PN order at which the orbital evolution is computed, a multiple of 0.5
    from 0 to 6.0.

Both PN orders are rounded down to a multiple of 0.5, and an order outside
the stated range raises an error.

The waveform has modes with :math:`2 \leq \ell \leq 8` and all values of
:math:`m`.  Its data type is :math:`h`, with the radius and the total mass
scaled out.  The modes are computed in the co-orbital frame, whose rotor at
each time is stored as the waveform's `Waveform.Frame`, so the waveform
initially has frame type ``Coorbital``; this frame is nearly, but not quite,
corotating.  To obtain the waveform in an inertial frame, apply
`Waveform.TransformToInertialFrame`; to obtain it in the corotating frame,
apply `Waveform.TransformToInertialFrame` and then
`Waveform.TransformToCorotatingFrame`.

Examples
--------
>>> import GWFrames
>>> W = GWFrames.PNWaveform("TaylorT4", 0.0, [0.0, 0.0, 0.9],
...                          [0.0, 0.0, 0.9], 0.01))doc"},
    {"PNWaveform::chi1", R"doc(Return the spin vector of black hole 1.

The dimensionless spin vector :math:`\vec{S}_1 / M_1^2` is given in the
inertial basis.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the vectors at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The spin vector at the given time step, or at each of the ``N`` time
    steps.)doc"},
    {"PNWaveform::chi2", R"doc(Return the spin vector of black hole 2.

The dimensionless spin vector :math:`\vec{S}_2 / M_2^2` is given in the
inertial basis.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the vectors at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The spin vector at the given time step, or at each of the ``N`` time
    steps.)doc"},
    {"PNWaveform::Omega_orb", R"doc(Return the orbital angular-velocity vector.

The orbital angular velocity is :math:`v^3 \hat{\ell}`, where :math:`v` is the
PN parameter and :math:`\hat{\ell}` is the direction normal to the orbital
plane.  It is given in the inertial basis, in units of :math:`1/M`.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the vectors at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The orbital angular velocity at the given time step, or at each of the
    ``N`` time steps.)doc"},
    {"PNWaveform::Omega_prec", R"doc(Return the precessional angular-velocity vector.

The precessional angular velocity is the angular velocity of the co-orbital
frame minus the orbital angular velocity, so that it describes the motion of
the orbital plane.  It is given in the inertial basis, in units of
:math:`1/M`.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the vectors at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The precessional angular velocity at the given time step, or at each of
    the ``N`` time steps.)doc"},
    {"PNWaveform::Omega_tot", R"doc(Return the total angular-velocity vector of the binary.

The total angular velocity is the sum of `PNWaveform.Omega_orb` and
`PNWaveform.Omega_prec`, which is the angular velocity of the co-orbital
frame.  It is given in the inertial basis, in units of :math:`1/M`.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the vectors at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The total angular velocity at the given time step, or at each of the
    ``N`` time steps.)doc"},
    {"PNWaveform::L", R"doc(Return the orbital angular-momentum vector.

The orbital angular momentum is computed from its PN expression, at the PN
order of the orbital evolution, and is given in the inertial basis in units
of :math:`M^2`.  Its leading term is :math:`M^2 \nu \hat{\ell} / v`, where
:math:`\nu` is the symmetric mass ratio; the spin-orbit terms make it
deviate slightly from the direction :math:`\hat{\ell}` of the orbital
angular velocity.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the vectors at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The orbital angular momentum at the given time step, or at each of the
    ``N`` time steps.)doc"},
    {"PNWaveform::Omega_orbMag", R"doc(Return the magnitude of the orbital angular velocity.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the magnitudes at all times are
    returned.

Returns
-------
float or numpy.ndarray of float, shape (N,)
    The magnitude of `PNWaveform.Omega_orb` at the given time step, or at
    each of the ``N`` time steps.)doc"},
    {"PNWaveform::Omega_precMag", R"doc(Return the magnitude of the precessional angular velocity.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the magnitudes at all times are
    returned.

Returns
-------
float or numpy.ndarray of float, shape (N,)
    The magnitude of `PNWaveform.Omega_prec` at the given time step, or at
    each of the ``N`` time steps.)doc"},
    {"PNWaveform::Omega_totMag", R"doc(Return the magnitude of the total angular velocity.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the magnitudes at all times are
    returned.

Returns
-------
float or numpy.ndarray of float, shape (N,)
    The magnitude of `PNWaveform.Omega_tot` at the given time step, or at
    each of the ``N`` time steps.)doc"},
    {"PNWaveform::LMag", R"doc(Return the magnitude of the orbital angular momentum.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the magnitudes at all times are
    returned.

Returns
-------
float or numpy.ndarray of float, shape (N,)
    The magnitude of `PNWaveform.L` at the given time step, or at each of the
    ``N`` time steps.)doc"},
    {"PNWaveform::chiHat1", R"doc(Return the direction of the spin of black hole 1.

The result is `PNWaveform.chi1` divided by its magnitude, so it is undefined
(NaN) when the spin is zero.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the directions at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The unit vector at the given time step, or at each of the ``N`` time
    steps.)doc"},
    {"PNWaveform::chiHat2", R"doc(Return the direction of the spin of black hole 2.

The result is `PNWaveform.chi2` divided by its magnitude, so it is undefined
(NaN) when the spin is zero.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the directions at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The unit vector at the given time step, or at each of the ``N`` time
    steps.)doc"},
    {"PNWaveform::OmegaHat_orb", R"doc(Return the direction of the orbital angular velocity.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the directions at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The unit vector along `PNWaveform.Omega_orb` at the given time step, or
    at each of the ``N`` time steps.)doc"},
    {"PNWaveform::OmegaHat_prec", R"doc(Return the direction of the precessional angular velocity.

The result is undefined (NaN) wherever the precessional angular velocity is
zero, as it is for a nonprecessing system.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the directions at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The unit vector along `PNWaveform.Omega_prec` at the given time step, or
    at each of the ``N`` time steps.)doc"},
    {"PNWaveform::OmegaHat_tot", R"doc(Return the direction of the total angular velocity.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the directions at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The unit vector along `PNWaveform.Omega_tot` at the given time step, or
    at each of the ``N`` time steps.)doc"},
    {"PNWaveform::LHat", R"doc(Return the direction of the orbital angular momentum.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the directions at all times are
    returned.

Returns
-------
numpy.ndarray of float, shape (3,) or (N, 3)
    The unit vector along `PNWaveform.L` at the given time step, or at each
    of the ``N`` time steps.)doc"},
    {"PNWaveform::Phi_orb", R"doc(Return the orbital phase.

The orbital phase is the integral of the orbital angular frequency, and is
zero at :math:`t = 0`, where the orbital frequency is ``Omega_orb_i``.

Parameters
----------
iTime : int, optional
    Index of the time step.  If omitted, the phases at all times are
    returned.

Returns
-------
float or numpy.ndarray of float, shape (N,)
    The orbital phase at the given time step, or at each of the ``N`` time
    steps.)doc"},
    {"PNWaveform::chi1Mag", R"doc(Return the magnitude of the spin of black hole 1.

Parameters
----------
iTime : int
    Index of the time step.

Returns
-------
float
    The magnitude of `PNWaveform.chi1` at the given time step.)doc"},
    {"PNWaveform::chi2Mag", R"doc(Return the magnitude of the spin of black hole 2.

Parameters
----------
iTime : int
    Index of the time step.

Returns
-------
float
    The magnitude of `PNWaveform.chi2` at the given time step.)doc"},
    {"WaveformAtAPointFT", R"doc(Fourier transform of a waveform measured by a detector at one point.

A `WaveformAtAPointFT` stores the Fourier transform of the real signal that a
detector would measure from a `Waveform` evaluated in one direction on the
sphere, in physical units.  It is used to compute signal-to-noise ratios,
inner products, and matches weighted by a detector's noise curve.

A `WaveformAtAPointFT` is constructed as

``WaveformAtAPointFT(W, Dt, Vartheta, Varphi, TotalMass, ...)``
    The Fourier transform of the waveform ``W`` at the point
    (``Vartheta``, ``Varphi``), computed as described below.  The optional
    parameters ``WindowNCycles``, ``DetectorResponseAmp``,
    ``DetectorResponsePhase``, and ``ExtraZeroPadPowers`` follow, in that
    order.

The construction takes these steps:

1. ``W`` is interpolated to times beginning with its first time and spaced
   uniformly by ``Dt``.  The number of samples is a power of 2, large enough
   to cover the duration of ``W``, and is multiplied further by
   :math:`2^{\mathtt{ExtraZeroPadPowers}}`; samples after the end of ``W``
   are zero.
2. The waveform is evaluated at the point (``Vartheta``, ``Varphi``) with
   `Waveform.EvaluateAtPoint`, giving a complex time series :math:`h(t)`.
3. The detector response :math:`F = F_+ + i F_\times`, given by its amplitude
   ``DetectorResponseAmp`` and phase ``DetectorResponsePhase``, is applied,
   and the real part :math:`\mathrm{Re}(F h)` is kept as the signal.
4. The signal is set to zero before its first zero crossing and is tapered
   smoothly to its full amplitude over the following ``WindowNCycles``
   cycles.
5. The times are converted to seconds using ``TotalMass``, and the signal is
   Fourier transformed.

The transform is normalized as an approximation to the continuous Fourier
transform (the discrete sum is multiplied by the time step in seconds), and
only the nonnegative frequencies are stored, from zero up to the Nyquist
frequency.  The values at the Nyquist frequency and the imaginary part at
zero frequency are set to zero.  No rescaling of the amplitude, by distance
or otherwise, is applied.

Parameters
----------
W : Waveform
    The waveform, whose times are in units of the total mass.
Dt : float
    The uniform time step, in units of the total mass, to which ``W`` is
    interpolated.
Vartheta : float
    The polar angle of the point at which the waveform is evaluated.
Varphi : float
    The azimuthal angle of the point at which the waveform is evaluated.
TotalMass : float
    The total mass of the system in solar masses, which sets the conversion
    of time to seconds.
WindowNCycles : int, default 1
    The number of cycles of the signal over which it is tapered at the start.
DetectorResponseAmp : float, default 1.0
    The amplitude of the complex detector response :math:`F_+ + i F_\times`.
DetectorResponsePhase : float, default 0.0
    The phase of the complex detector response :math:`F_+ + i F_\times`.
ExtraZeroPadPowers : int, default 0
    The number of additional factors of 2 by which the data are padded with
    zeros before the transform.)doc"},
    {"WaveformAtAPointFT::F", R"doc(Return the frequencies of the Fourier transform, in hertz.

The frequencies are nonnegative and uniformly spaced, starting at zero.

Parameters
----------
f : int, optional
    Index of the frequency.  If omitted, all frequencies are returned.

Returns
-------
float or numpy.ndarray of float
    The frequency with the given index, or all frequencies.)doc"},
    {"WaveformAtAPointFT::NFreq", R"doc(Return the number of frequencies stored.

Returns
-------
int
    The number of frequencies.)doc"},
    {"WaveformAtAPointFT::IsNormalized", R"doc(Return whether the transform has been normalized.

Returns
-------
bool
    True if `WaveformAtAPointFT.Normalize` has been applied to this
    object.)doc"},
    {"WaveformAtAPointFT::Dt", R"doc(Return the time step used before the Fourier transform.

Returns
-------
float
    The uniform time step, in units of the total mass, to which the waveform
    was interpolated.)doc"},
    {"WaveformAtAPointFT::Vartheta", R"doc(Return the polar angle of the point at which the waveform was evaluated.

Returns
-------
float
    The polar angle.)doc"},
    {"WaveformAtAPointFT::Varphi", R"doc(Return the azimuthal angle of the point at which the waveform was evaluated.

Returns
-------
float
    The azimuthal angle.)doc"},
    {"WaveformAtAPointFT::Re", R"doc(Return the real part of the Fourier transform.

Parameters
----------
f : int, optional
    Index of the frequency.  If omitted, the values at all frequencies are
    returned.

Returns
-------
float or numpy.ndarray of float
    The real part at the frequency with the given index, or at each of the
    frequencies `WaveformAtAPointFT.F`.)doc"},
    {"WaveformAtAPointFT::Im", R"doc(Return the imaginary part of the Fourier transform.

Parameters
----------
f : int, optional
    Index of the frequency.  If omitted, the values at all frequencies are
    returned.

Returns
-------
float or numpy.ndarray of float
    The imaginary part at the frequency with the given index, or at each of
    the frequencies `WaveformAtAPointFT.F`.)doc"},
    {"WaveformAtAPointFT::InversePSD", R"doc(Return the inverse noise power spectral density of a detector.

The inverse PSD :math:`1/S_n(f)` is evaluated at each of the frequencies
`WaveformAtAPointFT.F`, and may be passed as the ``InversePSD`` argument of
`WaveformAtAPointFT.SNR`, `WaveformAtAPointFT.InnerProduct`,
`WaveformAtAPointFT.Match`, and `WaveformAtAPointFT.Normalize`.  It is zero
(or negligibly small) outside the sensitive band of the detector.

Parameters
----------
Detector : str, default "AdvLIGO_ZeroDet_HighP"
    The name of the noise curve, which must be one of these:

    * ``"AdvLIGO_ZeroDet_HighP"``: Advanced LIGO, zero detuning and high
      power, from LIGO-T0900288-v3, between 10 Hz and 16384 Hz.
    * ``"AdvLIGO_ZeroDet_LowP"``: Advanced LIGO, zero detuning and low power,
      from LIGO-T0900288-v3, between 10 Hz and 16384 Hz.
    * ``"AdvLIGO_NSNSOptimal"``: a fit to the Advanced LIGO noise curve
      optimized for neutron-star binaries, between 10 Hz and 8192 Hz.
    * ``"IniLIGO_Approx"``: an approximation to the design goal of initial
      LIGO, between 40 Hz and 8192 Hz.
    * ``"Flat"``: a PSD equal to 1 at every frequency.

    Any other name raises an error.

Returns
-------
numpy.ndarray of float
    The inverse PSD at each frequency.)doc"},
    {"WaveformAtAPointFT::SNR", R"doc(Compute the signal-to-noise ratio of the waveform.

The signal-to-noise ratio is

.. math::

   \rho = \sqrt{4\, \Delta f \sum_k \frac{|\tilde{h}(f_k)|^2}{S_n(f_k)}},

where the sum runs over the stored frequencies :math:`f_k` and
:math:`\Delta f` is their spacing.

Parameters
----------
InversePSD : array_like of float, optional
    The inverse noise PSD at each of the frequencies `WaveformAtAPointFT.F`,
    as returned by `WaveformAtAPointFT.InversePSD`.  If it is given,
    ``Detector`` may not be.
Detector : str, default "AdvLIGO_ZeroDet_HighP"
    The name of the noise curve to use when ``InversePSD`` is not given; see
    `WaveformAtAPointFT.InversePSD` for the accepted names.

Returns
-------
float
    The signal-to-noise ratio.)doc"},
    {"WaveformAtAPointFT::InnerProduct", R"doc(Compute the noise-weighted inner product with another waveform.

The inner product is

.. math::

   (A|B) = 4\, \Delta f \sum_k
   \frac{\mathrm{Re}[\tilde{A}(f_k) \tilde{B}^*(f_k)]}{S_n(f_k)},

where :math:`A` is this waveform, the sum runs over the stored frequencies
:math:`f_k`, and :math:`\Delta f` is their spacing.  No maximization over
time or phase offsets is done; for that, use `WaveformAtAPointFT.Match`.

Parameters
----------
B : WaveformAtAPointFT
    The other waveform, which must have the same number of frequencies.
InversePSD : array_like of float
    The inverse noise PSD at each of the frequencies `WaveformAtAPointFT.F`,
    as returned by `WaveformAtAPointFT.InversePSD`.

Returns
-------
float
    The inner product.)doc"},
    {"WaveformAtAPointFT::Match", R"doc(Compute the match between this waveform and another.

The complex noise-weighted overlap of the two waveforms is computed as a
function of the time offset between them by an inverse Fourier transform,
and its magnitude is maximized over that offset; maximizing the magnitude
also maximizes over the phase offset.  The match is

.. math::

   4\, \Delta f \max_t \left| \sum_k
   \frac{\tilde{A}(f_k) \tilde{B}^*(f_k)}{S_n(f_k)}
   e^{2\pi i f_k t} \right|.

This is the conventional match, between 0 and 1, only if both waveforms have
been normalized with `WaveformAtAPointFT.Normalize`; otherwise, a warning is
printed and the unnormalized overlap is returned.  The two waveforms must have
the same number of frequencies and the same frequency spacing, or an error is
raised.

The method may be called as ``Match(B, InversePSD)`` or ``Match(B,
Detector)``, which return only the match, or as ``Match(B, InversePSD,
timeOffset, phaseOffset, match)`` or ``Match(B, timeOffset, phaseOffset,
match, Detector)``, which return the tuple ``(timeOffset, phaseOffset,
match)``.  In the latter forms, the input values of ``timeOffset``,
``phaseOffset``, and ``match`` are ignored; they serve only to select the
form that returns the offsets.

Parameters
----------
B : WaveformAtAPointFT
    The waveform with which to compute the match.
InversePSD : array_like of float, optional
    The inverse noise PSD at each of the frequencies `WaveformAtAPointFT.F`,
    as returned by `WaveformAtAPointFT.InversePSD`.  If it is given,
    ``Detector`` may not be.
timeOffset : float, optional
    Ignored; see above.
phaseOffset : float, optional
    Ignored; see above.
match : float, optional
    Ignored; see above.
Detector : str, default "AdvLIGO_ZeroDet_HighP"
    The name of the noise curve to use when ``InversePSD`` is not given; see
    `WaveformAtAPointFT.InversePSD` for the accepted names.

Returns
-------
timeOffset : float
    The time offset between the waveforms, in seconds, at which the overlap
    is maximized.  It is returned only in the tuple forms.
phaseOffset : float
    Half the phase of the complex overlap at that time offset.  It is
    returned only in the tuple forms.
match : float
    The match between the two waveforms.)doc"},
    {"WaveformAtAPointFT::Normalize", R"doc(Normalize the waveform to unit signal-to-noise ratio.

The Fourier transform is divided by the signal-to-noise ratio computed by
`WaveformAtAPointFT.SNR` with the given noise curve.  If the waveform has
already been normalized, it is left unchanged, even if a different noise
curve is given.

Parameters
----------
InversePSD : array_like of float, optional
    The inverse noise PSD at each of the frequencies `WaveformAtAPointFT.F`,
    as returned by `WaveformAtAPointFT.InversePSD`.  If it is given,
    ``Detector`` may not be.
Detector : str, default "AdvLIGO_ZeroDet_HighP"
    The name of the noise curve to use when ``InversePSD`` is not given; see
    `WaveformAtAPointFT.InversePSD` for the accepted names.

Returns
-------
WaveformAtAPointFT
    This object, modified in place.)doc"},
    {"WaveformAtAPointFT::ZeroAbove", R"doc(Set the Fourier transform to zero above a frequency.

Parameters
----------
Frequency : float
    The frequency, in hertz, above which the data are set to zero.

Returns
-------
WaveformAtAPointFT
    This object, modified in place.)doc"},
  };
  return docstrings;
}

}  // namespace GWFramesBindings

#endif  // GWFRAMES_BINDINGS_DOCSTRINGS_HPP
