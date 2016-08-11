#!/usr/bin/env perl 

use strict;
use warnings;
use utf8;
use Math::Complex;

my ($ns,$tau,$file1,$wmax,$dw,$w0,$ST)=@ARGV;

my $Tfin;
my $pi = acos(-1);

my $USAGE = "Calculates the Fourier transform of a time signal provided in file.";

defined($ns) or die    "$USAGE USAGE: $0 L tau file wmax dw w0 [window]\n";
defined($tau) or die   "$USAGE USAGE: $0 L tau file wmax dw w0 [window]\n";
defined($file1) or die "$USAGE USAGE: $0 L tau file wmax dw w0 [window]\n";
defined($wmax) or die  "$USAGE USAGE: $0 L tau file wmax dw w0 [window]\n";
defined($dw) or die    "$USAGE USAGE: $0 L tau file wmax dw w0 [window]\n";
defined($w0) or die    "$USAGE USAGE: $0 L tau file wmax dw w0 [window]\n";

open(my $fh, '<:encoding(UTF-8)', $file1) or die "Could not open file '$file1' $!";

my $i = 0;
my $time = 0.0;

my @valr;
my @vali;

while(<$fh>){
        chomp;
        my @temp = split;
        next if (/^#/);
	
        ($valr[$i],$vali[$i]) = ($temp[1],$temp[2]);

        $time = $temp[0];
	++$i;
}

my $tot = $i;
$Tfin = $time;

my $totw = int($wmax/$dw);

#print "tau=$tau, Tfin=$Tfin, tot= $tot, totw = $totw\n";

for (my $jp = 0; $jp < $totw; ++$jp){

        my $omega = $w0 + $jp * $dw;

	my $fre = 0.0;
        my $fim = 0.0;

	for (my $ip = 0; $ip < $tot; ++$ip){

		
		my $xt = $ip * $tau;
		my $filter = 0.5*(1.0+cos($xt*$pi/$Tfin));
		#$filter = exp(-$xt*$xt/(2.0*$ST*$ST));
		#$filter = 1.0;
		$fre += $tau*$filter*$valr[$ip]*cos($omega*$xt); 
                $fim += $tau*$filter*$vali[$ip]*cos($omega*$xt);

	}

	print " $omega $fre $fim\n";

}
