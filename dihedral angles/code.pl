use Math::Trig;
use Data::Dumper;
#rad2deg
open $file,("diala.sd") or die $!;

#data extraction
while(<$file>){
chomp;
 #1 N          38.4030   60.6080   35.4500 N.4       1 ALA    0.0000
if($_ =~/(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+/){

	push @x_cordinate, $1;
	push @y_cordinate, $2;
	push @z_cordinate, $3;
	}

if($_ =~/(\d+)\s+(\d+)\s+(\d+)\s+/){
	#print $1,"\t","\t",$2,"\t",$3,"\n";
	push @connections,"$1,$2";
	;
	}


}
#print @connections;
#data procssing 

$len = @connections;
#print"$len\n";
foreach $i(0..$len-1){

	@con=split/,/,$connections[$i];
	
	print ($con[0],"--",$con[1],"\n");
	
	$euclidean_dis = (($x_cordinate[$con[0]-1]-$x_cordinate[$con[1]-1])**2+($y_cordinate[$con[0]-1]-$y_cordinate[$con[1]-1])**2+($z_cordinate[$con[0]-1]-$z_cordinate[$con[1]-1])**2)**0.5;
	
	#print $euclidean_dis,"\n";

}

foreach $i (0..$len-2){

	my @coni=split/,/,@connections[$i];
	
	foreach $j ($i+1..$len){
	
		my @conj=split/,/,@connections[$j];
		
		if(($coni[0] == $conj[1]) or ($coni[1] == $conj[0]) or ($coni[0] == $conj[0]) or ($coni[1] == $conj[1]) ){
		
			my @con = (@coni,@conj);
			
			my @value = (unique(@con));
			
			my $v = join//,@con;
			
			if (!$angles{$v}){
			
				@{$angles{$v}}=@value;
				
				}
			}
		}
		
		
	}
#calculating bond angles
for $i(values %angles){
	
	$a = @{$i}[0]-1;
	$b = @{$i}[1]-1;
	$c = @{$i}[2]-1;
	#print $a+1,'--',$c+1,'--',$b+1,"\n";
	$ac = euc_distance($x_cordinate[$a],$x_cordinate[$c],$y_cordinate[$a],$y_cordinate[$c],$z_cordinate[$a],$z_cordinate[$c]);
	$bc = euc_distance($x_cordinate[$b],$x_cordinate[$c],$y_cordinate[$b],$y_cordinate[$c],$z_cordinate[$b],$z_cordinate[$c]);
	$ab = euc_distance($x_cordinate[$a],$x_cordinate[$b],$y_cordinate[$a],$y_cordinate[$b],$z_cordinate[$a],$z_cordinate[$b]);
	$angle = acos((($ac)**2+($bc)**2-($ab)**2)/(2*$ac*$bc));
	#$angle = acos(((@{$i}[0])**2+(@{$i}[2])**2-(@{$i}[1])**2)/2*(@{$i}[0])*(@{$i}[2]));

	#print rad2deg($angle),"\n";


}
#calculating dihedral angles
$num = 0;
for $i(values %angles){
	
	$a = @{$i}[0];
	$c = @{$i}[1];
	$b = @{$i}[2];
	
	for $j (values %angles){
	
		$aa = @{$j}[0];
		$cc = @{$j}[1];
		$bb = @{$j}[2];
		
		if (($a==$bb and $b==$cc)){
			$num++;
			push @{$dihedral[$num]},($aa,$a,$b,$c);
		}
		if(($aa==$b and $bb == $a)){
			$num++;
			push @{$dihedral[$num]},($cc,$a,$b,$c);
		}
		if (($c==$bb and $b==$aa)){
			$num++;
			push @{$dihedral[$num]},($a,$aa,$bb,$cc);
		}
		if(($cc==$b and $bb == $c)){
			$num++;
			push @{$dihedral[$num]},($aa,$c,$b,$a);
		}
	}
}


print Dumper(@dihedral);



sub unique {

	@words = @_;
	my @unique;
	my %seen;

	foreach my $value (@words) {
	
	  if ($seen{$value}!=1) {
	  
	    push @unique, $value;
	    
	    $seen{$value} = 1;
	  }
	  elsif($seen{$value} = 1){
	     $repeat = $value;
	  }
	}
	
	my $index = 0;
	$index++ until $unique[$index]== $repeat;
	splice(@unique, $index, 1);
	push @unique,$repeat;
	return @unique;
	
}

sub euc_distance{

	$x1= $_[0];
	$x2= $_[1];
	$y1= $_[2];
	$y2= $_[3];
	$z1= $_[4];
	$z2= $_[5];

	$dis = (($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2)**0.5;
	return $dis;	
}


