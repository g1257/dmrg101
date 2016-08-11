use strict;
use warnings;
use utf8;
use Math::Complex;

my ($i,$j,$dt,$tot)=@ARGV;

defined($tot) or die "USAGE: $0 site(i)[apply] site(j)[measure] dt tot\n";

if($i<0 || $i>5)  {die "USAGE: $0 site(i) should be 0 <= i <= 5\n";}
if($j<0 || $j>5)  {die "USAGE: $0 site(j) should be 0 <= j <= 5\n";}

if($dt<0 || $tot<0)  {die "USAGE: $0 dt is the time step length and tot is the total number of steps.\n Use for instance dt = 0.1 and tot = 1000\n";}

my $L = 6;
my $N = 4;
my @k;
my @en;
my @dens;
my $pi = acos(-1);

for (my $ip = 0; $ip<$L+1; $ip++){
                $k[$ip] = $ip*$pi/($L+1);
		$en[$ip] = -2*cos($k[$ip]);
}

my $norm = 2/($L+1);

for (my $jp = 1; $jp<$L+1; $jp++){
	my $sum = 0.0;
	for (my $ip = 1; $ip<$N+1; $ip++){
    		$sum += $norm*sin($k[$ip]*$jp)*sin($k[$ip]*$jp);
	}
        $dens[$jp] = $sum;
}

$i++;
$j++;

my $denstzero = $norm*(sin($k[$N+2]*$i)*sin($k[$N+2]*$i)+sin($k[$N+1]*$i)*sin($k[$N+1]*$i));

for (my $time = 0; $time<$tot; ++$time){

        my $t = $time*$dt;
	my $part1 = $norm*sin($k[$N+1]*$i)*sin($k[$N+1]*$i)*(1.0-$norm*sin($k[$N+2]*$j)*sin($k[$N+2]*$j));
        my $part2 = $norm*sin($k[$N+2]*$i)*sin($k[$N+2]*$i)*(1.0-$norm*sin($k[$N+1]*$j)*sin($k[$N+1]*$j));
	my $part3 = 2.0*$norm*$norm*sin($k[$N+2]*$i)*sin($k[$N+2]*$j)*sin($k[$N+1]*$i)*sin($k[$N+1]*$j)*cos(($en[$N+1]-$en[$N+2])*$t);
	my $de = ($part1 + $part2 + $part3)/$denstzero;
        print "$t $de\n";
}

