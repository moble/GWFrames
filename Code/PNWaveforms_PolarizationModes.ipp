// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

// Note: This header should only be included from PNWaveforms.cpp, and
// should therefore be invisible everywhere else.  You should not need
// to include this header directly.

#ifndef DOXYGEN

const int ellMax_PNWaveforms = 8;

// Construct PN waveform given dynamics, assuming standard orientation
class PNWaveformFromDynamics {
private:
  // These constants will be fixed once the mass-difference is known
  const double delta, nu, pownu2, pownu3;
  const double HhatL2M0Rev0, HhatL2M1Imv1, HhatL2M1Imv3, HhatL2M1Imv4, HhatL2M1Imv5, HhatL2M1Rev6, HhatL2M1Imv6, HhatL2M2Rev0, HhatL2M2Rev2, HhatL2M2Rev5, HhatL2M2Imv5, HhatL2M2Rev6, HhatL2M2Imv6, HhatL2M2Rev6lnv, HhatL2M2Rev7, HhatL2M2Imv7, HhatL3M0Imv5, HhatL3M1Imv1, HhatL3M1Imv3, HhatL3M1Imv4, HhatL3M1Imv5, HhatL3M1Rev6, HhatL3M1Imv6, HhatL3M2Rev2, HhatL3M2Rev4, HhatL3M2Rev5, HhatL3M2Imv5, HhatL3M2Rev6, HhatL3M3Imv1, HhatL3M3Imv3, HhatL3M3Imv4, HhatL3M3Imv5, HhatL3M3Rev6, HhatL3M3Imv6, HhatL4M0Rev0, HhatL4M1Imv3, HhatL4M1Imv5, HhatL4M1Rev6, HhatL4M1Imv6, HhatL4M2Rev2, HhatL4M2Rev4, HhatL4M2Rev5, HhatL4M2Imv5, HhatL4M2Rev6, HhatL4M3Imv3, HhatL4M3Imv5, HhatL4M3Rev6, HhatL4M3Imv6, HhatL4M4Rev2, HhatL4M4Rev4, HhatL4M4Rev5, HhatL4M4Imv5, HhatL4M4Rev6, HhatL5M1Imv3, HhatL5M1Imv5, HhatL5M1Rev6, HhatL5M1Imv6, HhatL5M2Rev4, HhatL5M2Rev6, HhatL5M3Imv3, HhatL5M3Imv5, HhatL5M3Rev6, HhatL5M3Imv6, HhatL5M4Rev4, HhatL5M4Rev6, HhatL5M5Imv3, HhatL5M5Imv5, HhatL5M5Rev6, HhatL5M5Imv6, HhatL6M1Imv5, HhatL6M2Rev4, HhatL6M2Rev6, HhatL6M3Imv5, HhatL6M4Rev4, HhatL6M4Rev6, HhatL6M5Imv5, HhatL6M6Rev4, HhatL6M6Rev6, HhatL7M1Imv5, HhatL7M2Rev6, HhatL7M3Imv5, HhatL7M4Rev6, HhatL7M5Imv5, HhatL7M6Rev6, HhatL7M7Imv5, HhatL8M2Rev6, HhatL8M4Rev6, HhatL8M6Rev6, HhatL8M8Rev6;
public:
  PNWaveformFromDynamics(const double idelta) :
    delta(idelta), nu((1.0-delta*delta)/4.0), pownu2(nu*nu), pownu3(nu*pownu2),
    HhatL2M0Rev0(-0.1458029608799511),
    HhatL2M1Imv1(0.3333333333333333*delta),
    HhatL2M1Imv3(0.011904761904761904*delta*(-17. + 20.*nu)),
    HhatL2M1Imv4(1.0471975511965976*delta),
    HhatL2M1Imv5(-0.0006613756613756613*delta*(172. + 2036.*nu - 237.*pownu2)),
    HhatL2M1Rev6(0.005952380952380952*delta*(-64.13400827807628 + 722.6355323334387*nu)),
    HhatL2M1Imv6(0.03739991254273563*delta*(-17. + 6.*nu)),
    HhatL2M2Rev0(1.),
    HhatL2M2Rev2(0.023809523809523808*(-107. + 55.*nu)),
    HhatL2M2Rev5(0.1495996501709425*(-107. + 34.*nu)),
    HhatL2M2Imv5(-24.*nu),
    HhatL2M2Rev6(32.35880116100775 - 4.147801379146346*nu - 7.309163059163059*pownu2 + 1.1487393779060446*pownu3),
    HhatL2M2Imv6(12.80573005463268),
    HhatL2M2Rev6lnv(-8.152380952380952),
    HhatL2M2Rev7(-0.004155545838081737*(2173. + 4918.*nu - 1120.*pownu2)),
    HhatL2M2Imv7(88.4753086419753*nu - 4.302645502645503*pownu2),
    HhatL3M0Imv5(-0.3703280399090206*nu),
    HhatL3M1Imv1(0.0222717701593687*delta),
    HhatL3M1Imv3(-0.014847846772912466*delta*(4. + nu)),
    HhatL3M1Imv4(0.06996882951511307*delta),
    HhatL3M1Imv5(-0.00011248368767357928*delta*(-607. + 272.*nu + 247.*pownu2)),
    HhatL3M1Rev6(0.0014847846772912468*delta*(-111.45177444479563 - 23.260151319598084*nu)),
    HhatL3M1Imv6(-0.01166147158585218*delta*(16. + 7.*nu)),
    HhatL3M2Rev2(0.2817180849095055*(1. - 3.*nu)),
    HhatL3M2Rev4(0.00313020094343895*(-193. + 725.*nu - 365.*pownu2)),
    HhatL3M2Rev5(1.7700869318701762*(1. - 3.*nu)),
    HhatL3M2Imv5(0.1690308509457033*(-5. + 22.*nu)),
    HhatL3M2Rev6(-0.00007114093053270341*(1451. + 17387.*nu - 100026.*pownu2 + 16023.*pownu3)),
    HhatL3M3Imv1(-0.7763237542601484*delta),
    HhatL3M3Imv3(-1.5526475085202969*delta*(-2. + nu)),
    HhatL3M3Imv4(-7.316679009572791*delta),
    HhatL3M3Imv5(0.002352496225030753*delta*(-369. + 3676.*nu - 887.*pownu2)),
    HhatL3M3Rev6(0.0006389495919836613*delta*(8588.637447565925 - 43669.239042837224*nu)),
    HhatL3M3Imv6(3.6583395047863956*delta*(8. - 3.*nu)),
    HhatL4M0Rev0(-0.0014029896452114037),
    HhatL4M1Imv3(0.0037646162621052135*delta*(1. - 2.*nu)),
    HhatL4M1Imv5(0.00002851982016746374*delta*(-404. + 1011.*nu - 332.*pownu2)),
    HhatL4M1Rev6(0.00012548720873684045*delta*(105.58883083359672 - 1744.1776616671934*nu)),
    HhatL4M1Imv6(0.011826890792614406*delta*(1. - 2.*nu)),
    HhatL4M2Rev2(0.03549314249999666*(1. - 3.*nu)),
    HhatL4M2Rev4(0.00010755497727271716*(4025.*nu - 57.*(23. + 5.*pownu2))),
    HhatL4M2Rev5(0.22300999146161035*(1. - 3.*nu)),
    HhatL4M2Imv5(0.14907119849998596*(-1. + 4.*nu)),
    HhatL4M2Rev6(1.9698713786211933e-8*(9.342351e6 - 3.8225313e7*nu + 2.803171e7*pownu2 + 2.707215e6*pownu3)),
    HhatL4M3Imv3(0.26892643710023856*delta*(-1. + 2.*nu)),
    HhatL4M3Imv5(-0.0020373214931836254*delta*(-468. + 1267.*nu - 524.*pownu2)),
    HhatL4M3Rev6(0.0003320079470373316*delta*(-3213.439574594321 + 12359.879149188642*nu)),
    HhatL4M3Imv6(2.534571957450561*delta*(-1. + 2.*nu)),
    HhatL4M4Rev2(0.751248226425348*(-1. + 3.*nu)),
    HhatL4M4Rev4(-0.002276509777046509*(6365.*nu - 3.*(593. + 875.*pownu2))),
    HhatL4M4Rev5(9.44046363664094*(-1. + 3.*nu)),
    HhatL4M4Imv5(0.0187812056606337*(114.1929022208175 - 527.5787066624525*nu)),
    HhatL4M4Rev6(4.1694318260925076e-7*(-9.618039e6 + 6.8551497e7*nu - 1.1309683e8*pownu2 + 2.3740185e7*pownu3)),
    HhatL5M1Imv3(0.00017696083036028665*delta*(1. - 2.*nu)),
    HhatL5M1Imv5(4.537457188725299e-6*delta*(-179. + 352.*nu - 4.*pownu2)),
    HhatL5M1Rev6(2.528011862289809e-6*delta*(278.0406052783923 - 8958.081210556786*nu)),
    HhatL5M1Imv6(0.0005559388446330261*delta*(1. - 2.*nu)),
    HhatL5M2Rev4(-0.009988146110566549*(-1. + 5.*nu - 5.*pownu2)),
    HhatL5M2Rev6(0.000010975984736886318*(-3911. + 21553.*nu - 28910.*pownu2 + 8085.*pownu3)),
    HhatL5M3Imv3(0.046446908841268335*delta*(-1. + 2.*nu)),
    HhatL5M3Imv5(-0.001190946380545342*delta*(-207. + 464.*nu - 88.*pownu2)),
    HhatL5M3Rev6(9.101882978888565e-7*delta*(-271701.69319944223 + 923537.3863988845*nu)),
    HhatL5M3Imv6(0.4377518027930503*delta*(-1. + 2.*nu)),
    HhatL5M4Rev4(0.2767996245907637*(-1. + 5.*nu - 5.*pownu2)),
    HhatL5M4Rev6(-0.0003041754116382019*(-4451. + 25333.*nu - 36470.*pownu2 + 11865.*pownu3)),
    HhatL5M5Imv3(-0.8013768943966975*delta*(-1. + 2.*nu)),
    HhatL5M5Imv5(0.020548125497351216*delta*(-263. + 688.*nu - 256.*pownu2)),
    HhatL5M5Rev6(0.0000183171861576388*delta*(164747.8048050571 - 679921.6096101142*nu)),
    HhatL5M5Imv6(-12.587998820966341*delta*(-1. + 2.*nu)),
    HhatL6M1Imv5(0.000023582988833355464*delta*(1. - 4.*nu + 3.*pownu2)),
    HhatL6M2Rev4(-0.0008352507379744677*(-1. + 5.*nu - 5.*pownu2)),
    HhatL6M2Rev6(0.000059660766998176277*(-81. + 413.*nu - 448.*pownu2 + 49.*pownu3)),
    HhatL6M3Imv5(0.01630976217812644*delta*(-1. + 4.*nu - 3.*pownu2)),
    HhatL6M4Rev4(0.058558165806265966*(-1. + 5.*nu - 5.*pownu2)),
    HhatL6M4Rev6(-0.004182726129018997*(-93. + 497.*nu - 616.*pownu2 + 133.*pownu3)),
    HhatL6M5Imv5(-0.2993579796535637*delta*(-1. + 4.*nu - 3.*pownu2)),
    HhatL6M6Rev4(0.9031413708076581*(1. - 5.*nu + 5.*pownu2)),
    HhatL6M6Rev6(0.06451009791483273*(-113. + 637.*nu - 896.*pownu2 + 273.*pownu3)),
    HhatL7M1Imv5(8.175930333399787e-7*delta*(1. - 4.*nu + 3.*pownu2)),
    HhatL7M2Rev6(0.0001922578318979773*(1. - 7.*nu + 14.*pownu2 - 7.*pownu3)),
    HhatL7M3Imv5(0.0018582230503756002*delta*(-1. + 4.*nu - 3.*pownu2)),
    HhatL7M4Rev6(0.023085290615904162*(-1. + 7.*nu - 14.*pownu2 + 7.*pownu3)),
    HhatL7M5Imv5(-0.07338616249054006*delta*(-1. + 4.*nu - 3.*pownu2)),
    HhatL7M6Rev6(-0.3352043015692001*(-1. + 7.*nu - 14.*pownu2 + 7.*pownu3)),
    HhatL7M7Imv5(1.0542244493439168*delta*(-1. + 4.*nu - 3.*pownu2)),
    HhatL8M2Rev6(0.000012039652448587866*(1. - 7.*nu + 14.*pownu2 - 7.*pownu3)),
    HhatL8M4Rev6(0.0032325872683449304*(-1. + 7.*nu - 14.*pownu2 + 7.*pownu3)),
    HhatL8M6Rev6(-0.0921843838337631*(-1. + 7.*nu - 14.*pownu2 + 7.*pownu3)),
    HhatL8M8Rev6(1.2608629579851318*(-1. + 7.*nu - 14.*pownu2 + 7.*pownu3))
  { }
  // This is the method to call to replace 'modes' with the waveform polarization data
  void operator()(const double v, const double chis, const double chia, std::vector<std::complex<double> >& modes) {
    double Re,Im;
    
    // Compute various v- or spin-dependent terms
    const double lnv = std::log(v);
    const double powv2 = v*v;
    const double powv3 = v*powv2;
    const double powv4 = powv2*powv2;
    const double powv5 = powv2*powv3;
    const double powv6 = powv3*powv3;
    const double NormalizationFactor = 2*nu*sqrt(16*M_PI/5.0)*powv2;
    const double HhatL2M1Imv2 = 0.5*(chia + chis*delta);
    const double HhatL2M1Rev4 = 0.023809523809523808*(chia*(7. - 205.*nu) + delta*(26.408121055678468 + 7.*chis - 33.*chis*nu));
    const double HhatL2M2Rev3 = 0.3333333333333333*(18.84955592153876 - 4.*chia*delta + 4.*chis*(-1. + nu));
    const double HhatL2M2Rev4 = 0.0006613756613756613*(-2173. - 7.*(1069. + 432.*chia*chia - 432.*chis*chis)*nu + 2047.*pownu2);
    const double HhatL3M1Rev4 = 0.0022271770159368698*(delta*(27.862943611198908 + chis*(20. - 65.*nu)) + chia*(20. - 55.*nu));
    const double HhatL3M2Rev3 = 1.126872339638022*chis*nu;
    const double HhatL3M3Rev4 = -0.07763237542601484*(5.*chia*(-4. + 19.*nu) + delta*(17.672093513510138 + 5.*chis*(-4. + 5.*nu)));
    const double HhatL4M1Rev4 = 0.009411540655263034*(chia - 1.*chis*delta)*nu;
    const double HhatL4M3Rev4 = 0.6723160927505963*(chia - 1.*chis*delta)*nu;
    
    // Make sure 'modes' is the right size
    modes.resize(ellMax_PNWaveforms*(2+ellMax_PNWaveforms) - 3);
    
    ///////////////////////////////////////////////
    // Step through 'modes' calculating the data //
    ///////////////////////////////////////////////
    
    // (2,+-0)
    Re = HhatL2M0Rev0;
    Im = 0;
    modes[2] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (2,+-1)
    Re = (HhatL2M1Rev4 + HhatL2M1Rev6*powv2)*powv4;
    Im = v*(HhatL2M1Imv1 + v*(HhatL2M1Imv2 + v*(HhatL2M1Imv3 + v*(HhatL2M1Imv4 + v*(HhatL2M1Imv5 + HhatL2M1Imv6*v)))));
    modes[3] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[1] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (2,+-2)
    Re = HhatL2M2Rev0 + powv2*(HhatL2M2Rev2 + v*(HhatL2M2Rev3 + v*(HhatL2M2Rev4 + v*(HhatL2M2Rev5 + v*(HhatL2M2Rev6 + HhatL2M2Rev6lnv*lnv + HhatL2M2Rev7*v)))));
    Im = powv5*(HhatL2M2Imv5 + v*(HhatL2M2Imv6 + HhatL2M2Imv7*v));
    modes[4] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[0] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (3,+-0)
    Re = 0;
    Im = HhatL3M0Imv5*powv5;
    modes[8] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (3,+-1)
    Re = (HhatL3M1Rev4 + HhatL3M1Rev6*powv2)*powv4;
    Im = v*(HhatL3M1Imv1 + powv2*(HhatL3M1Imv3 + v*(HhatL3M1Imv4 + v*(HhatL3M1Imv5 + HhatL3M1Imv6*v))));
    modes[9] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[7] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (3,+-2)
    Re = powv2*(HhatL3M2Rev2 + v*(HhatL3M2Rev3 + v*(HhatL3M2Rev4 + v*(HhatL3M2Rev5 + HhatL3M2Rev6*v))));
    Im = HhatL3M2Imv5*powv5;
    modes[10] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[6] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (3,+-3)
    Re = (HhatL3M3Rev4 + HhatL3M3Rev6*powv2)*powv4;
    Im = v*(HhatL3M3Imv1 + powv2*(HhatL3M3Imv3 + v*(HhatL3M3Imv4 + v*(HhatL3M3Imv5 + HhatL3M3Imv6*v))));
    modes[11] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[5] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (4,+-0)
    Re = HhatL4M0Rev0;
    Im = 0;
    modes[16] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (4,+-1)
    Re = (HhatL4M1Rev4 + HhatL4M1Rev6*powv2)*powv4;
    Im = powv3*(HhatL4M1Imv3 + powv2*(HhatL4M1Imv5 + HhatL4M1Imv6*v));
    modes[17] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[15] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (4,+-2)
    Re = powv2*(HhatL4M2Rev2 + powv2*(HhatL4M2Rev4 + v*(HhatL4M2Rev5 + HhatL4M2Rev6*v)));
    Im = HhatL4M2Imv5*powv5;
    modes[18] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[14] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (4,+-3)
    Re = (HhatL4M3Rev4 + HhatL4M3Rev6*powv2)*powv4;
    Im = powv3*(HhatL4M3Imv3 + powv2*(HhatL4M3Imv5 + HhatL4M3Imv6*v));
    modes[19] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[13] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (4,+-4)
    Re = powv2*(HhatL4M4Rev2 + powv2*(HhatL4M4Rev4 + v*(HhatL4M4Rev5 + HhatL4M4Rev6*v)));
    Im = HhatL4M4Imv5*powv5;
    modes[20] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[12] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (5,+-0)
    Re = 0;
    Im = 0;
    modes[26] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (5,+-1)
    Re = HhatL5M1Rev6*powv6;
    Im = powv3*(HhatL5M1Imv3 + powv2*(HhatL5M1Imv5 + HhatL5M1Imv6*v));
    modes[27] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[25] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (5,+-2)
    Re = (HhatL5M2Rev4 + HhatL5M2Rev6*powv2)*powv4;
    Im = 0;
    modes[28] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[24] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (5,+-3)
    Re = HhatL5M3Rev6*powv6;
    Im = powv3*(HhatL5M3Imv3 + powv2*(HhatL5M3Imv5 + HhatL5M3Imv6*v));
    modes[29] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[23] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (5,+-4)
    Re = (HhatL5M4Rev4 + HhatL5M4Rev6*powv2)*powv4;
    Im = 0;
    modes[30] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[22] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (5,+-5)
    Re = HhatL5M5Rev6*powv6;
    Im = powv3*(HhatL5M5Imv3 + powv2*(HhatL5M5Imv5 + HhatL5M5Imv6*v));
    modes[31] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[21] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (6,+-0)
    Re = 0;
    Im = 0;
    modes[38] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (6,+-1)
    Re = 0;
    Im = HhatL6M1Imv5*powv5;
    modes[39] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[37] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (6,+-2)
    Re = (HhatL6M2Rev4 + HhatL6M2Rev6*powv2)*powv4;
    Im = 0;
    modes[40] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[36] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (6,+-3)
    Re = 0;
    Im = HhatL6M3Imv5*powv5;
    modes[41] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[35] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (6,+-4)
    Re = (HhatL6M4Rev4 + HhatL6M4Rev6*powv2)*powv4;
    Im = 0;
    modes[42] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[34] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (6,+-5)
    Re = 0;
    Im = HhatL6M5Imv5*powv5;
    modes[43] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[33] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (6,+-6)
    Re = (HhatL6M6Rev4 + HhatL6M6Rev6*powv2)*powv4;
    Im = 0;
    modes[44] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[32] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (7,+-0)
    Re = 0;
    Im = 0;
    modes[52] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (7,+-1)
    Re = 0;
    Im = HhatL7M1Imv5*powv5;
    modes[53] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[51] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (7,+-2)
    Re = HhatL7M2Rev6*powv6;
    Im = 0;
    modes[54] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[50] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (7,+-3)
    Re = 0;
    Im = HhatL7M3Imv5*powv5;
    modes[55] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[49] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (7,+-4)
    Re = HhatL7M4Rev6*powv6;
    Im = 0;
    modes[56] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[48] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (7,+-5)
    Re = 0;
    Im = HhatL7M5Imv5*powv5;
    modes[57] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[47] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (7,+-6)
    Re = HhatL7M6Rev6*powv6;
    Im = 0;
    modes[58] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[46] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (7,+-7)
    Re = 0;
    Im = HhatL7M7Imv5*powv5;
    modes[59] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[45] = std::complex<double>(-NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (8,+-0)
    Re = 0;
    Im = 0;
    modes[68] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    
    // (8,+-1)
    Re = 0;
    Im = 0;
    modes[69] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[67] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (8,+-2)
    Re = HhatL8M2Rev6*powv6;
    Im = 0;
    modes[70] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[66] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (8,+-3)
    Re = 0;
    Im = 0;
    modes[71] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[65] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (8,+-4)
    Re = HhatL8M4Rev6*powv6;
    Im = 0;
    modes[72] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[64] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (8,+-5)
    Re = 0;
    Im = 0;
    modes[73] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[63] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (8,+-6)
    Re = HhatL8M6Rev6*powv6;
    Im = 0;
    modes[74] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[62] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (8,+-7)
    Re = 0;
    Im = 0;
    modes[75] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[61] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
    // (8,+-8)
    Re = HhatL8M8Rev6*powv6;
    Im = 0;
    modes[76] = std::complex<double>(NormalizationFactor*Re, NormalizationFactor*Im);
    modes[60] = std::complex<double>(NormalizationFactor*Re, -NormalizationFactor*Im);
    
  }
};


#endif // DOXYGEN
