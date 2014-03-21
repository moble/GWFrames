#include "fft.hpp"

#include "NumericalRecipes.hpp"
#include "VectorFunctions.hpp"
#include "Utilities.hpp"

using namespace std;
namespace WU = WaveformUtilities;

vector<double> WU::TimeToFrequency(const vector<double>& Time) {
  /// This returns the double-sided frequency-space equivalent of a time vector
  const unsigned int N = Time.size();
  if(N&(N-1)) {
    cerr << "\nN=" << N << " is not a power of 2." << endl;
    Throw1WithMessage("Bad size");
  }
  const double df = 1.0 / (N*(Time[1]-Time[0]));
  vector<double> Freq(N, 0.0);
  for(unsigned int i=0; i<N/2;  ++i) {
    Freq[i] = i*df;
  }
  for(unsigned int i=N/2; i<N; ++i) {
    Freq[i] = i*df - N*df;
  }
  return Freq;
}

vector<double> WU::TimeToPositiveFrequencies(const vector<double>& Time) {
  /// This returns the single-sided frequency-space equivalent of a time vector
  const unsigned int N = Time.size();
  if(N&(N-1)) {
    cerr << "\nN=" << N << " is not a power of 2." << endl;
    Throw1WithMessage("Bad size");
  }
  const unsigned int n = 1 + (N/2);
  const double df = 1.0 / (N*(Time[1]-Time[0]));
  vector<double> Freq(n, 0.0);
  for(unsigned int i=0; i<=N/2;  ++i) {
    Freq[i] = i*df;
  }
  return Freq;
}

void WU::WrapVecDoub::validate() {
  if (n&(n-1)) {
    cerr << "\nn=" << n << endl;
    Throw1WithMessage("vector size n must be power of 2");
  }
}

void four1(vector<double>& data, const int isign);

void WU::dft(vector<double>& data) {
  four1(data, -1);
  return;
}

void WU::idft(vector<double>& data) {
  four1(data, 1);
  return;
}

void realft(VecDoub_IO &Data, const Int isign);

void  WU::realdft(std::vector<double>& data) {
  realft(data, 1);
  for(unsigned int i=3; i<data.size(); ++i) {
    data[i++] *= -1.0;
  }
  return;
}


//// Numerical Recipes routines
void four1(Doub *data, const Int n, const Int isign) {
  Int nn,mmax,m,j,istep,i;
  Doub wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
  if (n<2 || n&(n-1)) Throw1WithMessage("n must be power of 2 in four1");
  nn = n << 1;
  j = 1;
  for (i=1;i<nn;i+=2) {
    if (j > i) {
      SWAP(data[j-1],data[i-1]);
      SWAP(data[j],data[i]);
    }
    m=n;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (nn > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=nn;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j-1]-wi*data[j];
	tempi=wr*data[j]+wi*data[j-1];
	data[j-1]=data[i-1]-tempr;
	data[j]=data[i]-tempi;
	data[i-1] += tempr;
	data[i] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
void four1(VecDoub_IO &data, const Int isign) {
  four1(&data[0],data.size()/2,isign);
}

void realft(VecDoub_IO &data, const Int isign) {
  Int i,i1,i2,i3,i4,n=data.size();
  Doub c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp;
  Doub theta=3.141592653589793238/Doub(n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  for (i=1;i<(n>>2);i++) {
    i2=1+(i1=i+i);
    i4=1+(i3=n-i1);
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r= -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4]= -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[0] = (h1r=data[0])+data[1];
    data[1] = h1r-data[1];
  } else {
    data[0]=c1*((h1r=data[0])+data[1]);
    data[1]=c1*(h1r-data[1]);
    four1(data,-1);
  }
}

void convlv(VecDoub_I &data, VecDoub_I &respns, const Int isign, VecDoub_O &ans) {
  Int i,no2,n=data.size(),m=respns.size();
  Doub mag2,tmp;
  VecDoub temp(n);
  temp[0]=respns[0];
  for (i=1;i<(m+1)/2;i++) {
    temp[i]=respns[i];
    temp[n-i]=respns[m-i];
  }
  for (i=(m+1)/2;i<n-(m-1)/2;i++)
    temp[i]=0.0;
  for (i=0;i<n;i++)
    ans[i]=data[i];
  realft(ans,1);
  realft(temp,1);
  no2=n>>1;
  if (isign == 1) {
    for (i=2;i<n;i+=2) {
      tmp=ans[i];
      ans[i]=(ans[i]*temp[i]-ans[i+1]*temp[i+1])/no2;
      ans[i+1]=(ans[i+1]*temp[i]+tmp*temp[i+1])/no2;
    }
    ans[0]=ans[0]*temp[0]/no2;
    ans[1]=ans[1]*temp[1]/no2;
  } else if (isign == -1) {
    for (i=2;i<n;i+=2) {
      if ((mag2=SQR(temp[i])+SQR(temp[i+1])) == 0.0)
	Throw1WithMessage("Deconvolving at response zero in convlv");
      tmp=ans[i];
      ans[i]=(ans[i]*temp[i]+ans[i+1]*temp[i+1])/mag2/no2;
      ans[i+1]=(ans[i+1]*temp[i]-tmp*temp[i+1])/mag2/no2;
    }
    if (temp[0] == 0.0 || temp[1] == 0.0)
      Throw1WithMessage("Deconvolving at response zero in convlv");
    ans[0]=ans[0]/temp[0]/no2;
    ans[1]=ans[1]/temp[1]/no2;
  } else Throw1WithMessage("No meaning for isign in convlv");
  realft(ans,-1);
}

vector<double> convlv(const vector<double>& data, const vector<double>& respns, const int isign) {
  vector<double> ans(data.size());
  convlv(data, respns, isign, ans);
  return ans;
}
