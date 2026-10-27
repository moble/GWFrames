#! /usr/bin/perl

$NoiseCurve_txt = $ARGV[0];
($NoiseCurve) = ($NoiseCurve_txt =~ /(.*)\.txt/);
${NOISECURVE_HPP} = "${NoiseCurve}_HPP";
print STDERR "Trying to open |${NoiseCurve_txt}|\n";
print STDERR "NoiseCurve = |${NoiseCurve}|\n";
print STDERR "NOISECURVE_HPP = |${NOISECURVE_HPP}|\n";
open FILE, "${NoiseCurve_txt}" or die $!;

print "#ifndef ${NOISECURVE_HPP}\n";
print "#define ${NOISECURVE_HPP}\n";
print "\n";
print "#include <vector>\n";
print "\n";
print "std::vector<double> ${NoiseCurve}LogF(3000);\n";
print "std::vector<double> ${NoiseCurve}LogPSD(3000);\n";
print "\n";

$i = 0;
while ($Line = <FILE>) {
    chomp ($Line);
    ($F, $PSD) = ($Line =~ / *([\d\.e+-]*) *([\d\.e+-]*)/);
    $LogF = log($F);
    $LogPSD = 2.0*log($PSD);
    print "${NoiseCurve}LogF[${i}]=${LogF}; ${NoiseCurve}LogPSD[${i}]=${LogPSD};\n";
    $i = $i+1;
}

print "\n";
print "#endif // ${NOISECURVE_HPP}\n";

close FILE;
