#!/usr/bin/env bash
# Build static, position-independent copies of GSL and FFTW, so that the
# wheels built by cibuildwheel do not depend on shared libraries.  The
# installation prefix is taken from DEPS_PREFIX; point CMAKE_PREFIX_PATH at
# the same directory when building GWFrames.
#
# GSL is compiled with -fexceptions because the GSL error handlers installed
# by GWFrames throw C++ exceptions, which must unwind through GSL's C frames.

set -euo pipefail

PREFIX="${DEPS_PREFIX:-/tmp/gwframes-deps}"

GSL_VERSION=2.8
GSL_SHA256=6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190
FFTW_VERSION=3.3.11
FFTW_SHA256=5630c24cdeb33b131612f7eb4b1a9934234754f9f388ff8617458d0be6f239a1

if [[ -f "${PREFIX}/lib/libgsl.a" && -f "${PREFIX}/lib/libfftw3.a" ]]; then
  echo "GSL and FFTW are already installed in ${PREFIX}"
  exit 0
fi

JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

sha256() {
  if command -v sha256sum >/dev/null; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

fetch() {
  local url="$1" expected="$2" file="${WORK}/$(basename "$1")"
  curl -fsSL --retry 5 -o "${file}" "${url}"
  local actual
  actual="$(sha256 "${file}")"
  if [[ "${actual}" != "${expected}" ]]; then
    echo "Checksum mismatch for ${url}: expected ${expected}, got ${actual}" >&2
    exit 1
  fi
  tar -xzf "${file}" -C "${WORK}"
}

fetch "https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz" "${GSL_SHA256}"
(
  cd "${WORK}/gsl-${GSL_VERSION}"
  CFLAGS="-O2 -fPIC -fexceptions" ./configure --prefix="${PREFIX}" --disable-shared --enable-static --with-pic
  make -j"${JOBS}"
  make install
)

fetch "https://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz" "${FFTW_SHA256}"
(
  cd "${WORK}/fftw-${FFTW_VERSION}"
  CFLAGS="-O2 -fPIC" ./configure --prefix="${PREFIX}" --disable-shared --enable-static --with-pic --disable-fortran
  make -j"${JOBS}"
  make install
)
