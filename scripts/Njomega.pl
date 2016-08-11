use strict;
use warnings;
use utf8;
use Math::Complex;

my ($site,$eta,$dw,$tot)=@ARGV;

defined($tot) or die "USAGE: $0 site eta domega tot\n";

if($site<0 || $site>5)  {die "USAGE: $0 site should be 0 <= i <= 5\n";}

$site++;

my $L = 6;
my $N = 4;
my @k;
my @en;
my @dens;
my $pi = acos(-1);
#print STDERR "PI=$pi\n";

for (my $ip = 0; $ip<$L+1; $ip++){
                $k[$ip] = $ip*$pi/($L+1);
		$en[$ip] = -2*cos($k[$ip]);
#		print STDERR "$k[$ip]\n";
}

my $norm = 2/($L+1);

for (my $jp = 1; $jp<$L+1; $jp++){
	my $sum = 0.0;
	for (my $ip = 1; $ip<$N+1; $ip++){
    		$sum += $norm*sin($k[$ip]*$jp)*sin($k[$ip]*$jp);
	}
        $dens[$jp] = $sum;
        #print STDERR "$jp $dens[$jp]\n";
}


my $center = ($L/2.0);

for (my $om = 0; $om<$tot; ++$om){

        my $omega = $om*$dw;

	for (my $ip = 1; $ip<$L+1; $ip++){
		
	        my $sumre = 0.0;
        	my $sumim = 0.0;

		for (my $jp = 1; $jp<$L+1; $jp++){

                        my $lorentz5re = ($omega-$en[$N+1])/(($omega-$en[$N+1])*($omega-$en[$N+1])+$eta*$eta);
                        my $lorentz6re = ($omega-$en[$N+2])/(($omega-$en[$N+2])*($omega-$en[$N+2])+$eta*$eta);

			my $lorentz5im = $eta/(($omega-$en[$N+1])*($omega-$en[$N+1])+$eta*$eta);
                        my $lorentz6im = $eta/(($omega-$en[$N+2])*($omega-$en[$N+2])+$eta*$eta);

                        $sumre += $norm*cos($k[$ip]*($jp-$center))*(sin($k[$N+1]*$center)*sin($k[$N+1]*$jp)*$lorentz5re+sin($k[$N+2]*$center)*sin($k[$N+2]*$jp)*$lorentz6re);
			$sumim += $norm*cos($k[$ip]*($jp-$center))*(sin($k[$N+1]*$center)*sin($k[$N+1]*$jp)*$lorentz5im+sin($k[$N+2]*$center)*sin($k[$N+2]*$jp)*$lorentz6im);
		}
                #print STDERR "$k[$ip] $omega $sumre $sumim\n";
	}
}

for (my $om = 0; $om<$tot; ++$om){

        my $omega = $om*$dw;

        #for (my $ip = 1; $ip<$L+1; $ip++){

                my $sumre = 0.0;
                my $sumim = 0.0;

                #for (my $jp = 1; $jp<$L+1; $jp++){

		my $jp =$L/2;
                        my $lorentz5re = ($omega-$en[$N+1])/(($omega-$en[$N+1])*($omega-$en[$N+1])+$eta*$eta);
                        my $lorentz6re = ($omega-$en[$N+2])/(($omega-$en[$N+2])*($omega-$en[$N+2])+$eta*$eta);

                        my $lorentz5im = $eta/(($omega-$en[$N+1])*($omega-$en[$N+1])+$eta*$eta);
                        my $lorentz6im = $eta/(($omega-$en[$N+2])*($omega-$en[$N+2])+$eta*$eta);

                        $sumre += $norm*(sin($k[$N+1]*$site)*sin($k[$N+1]*$site)*$lorentz5re+sin($k[$N+2]*$site)*sin($k[$N+2]*$site)*$lorentz6re);
                        $sumim += $norm*(sin($k[$N+1]*$site)*sin($k[$N+1]*$site)*$lorentz5im+sin($k[$N+2]*$site)*sin($k[$N+2]*$site)*$lorentz6im);

                #}
                print "$omega $sumim\n";
        #}
}

