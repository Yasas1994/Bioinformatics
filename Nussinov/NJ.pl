#NJ algorithm
$sequence="CGACUUAGUCG";
@seq=split//,uc($sequence);

require List::Util;

#initiazes the similarity matrix 
foreach $i(1..length $sequence){
	foreach $j(1..$i){
		#print "\t",$seq[$i-1],"\n\n";
		#print $seq[$j-1],"\n";
		$smatrix[$i][$j]=0;
		if ($i==$j){
			$matrix[$i][$j]=0;
			}
			
		elsif (($seq[$i-1] eq "G") and ($seq[$j-1] eq "C") or ($seq[$i-1] eq "C") and ($seq[$j-1] eq "G")){
			$matrix[$i][$j]=3;}
		
		elsif (($seq[$i-1] eq "A") and ($seq[$j-1] eq "U") or ($seq[$i-1] eq "U") and ($seq[$j-1] eq "A")){
			$matrix[$i][$j]=2;}
			
		elsif (($seq[$i-1] eq "G") and ($seq[$j-1] eq "U")or ($seq[$i-1] eq "U") and ($seq[$j-1] eq "G")){
			$matrix[$i][$j]=1;}
		}}

#filling S matrix through diagnals
foreach $ii(reverse 1..((length $sequence) -1)){
	$x++;
	foreach $jj(1..$ii){
		
			#print $jj+$x,",", $jj,",",$x,"\n";
			# D(i-1,j+1)+w(i,j)
			$diagonal = $smatrix[($jj+$x)-1][($jj)+1]+ $matrix[$jj+$x][$jj];
			#D(i-1,j)
			$row = $smatrix[($jj+$x)-1][$jj];
			#D(i,j+1)
			$col=$smatrix[$jj+$x][($jj)+1];

	

		        my @max=();
			foreach my $k(reverse (($jj)+1)..(($jj+$x)-1)){
				#print $k,"\n";
				$v=$smatrix[$jj+$x][$k]+$smatrix[($k-1)][$jj],"\n";
				push @max, $v;
				@max=sort{$b<=>$a}@max;
				$out=$max[0];
				}
				

		$smatrix[$jj+$x][$jj]=List::Util::max($diagonal,$row,$col,$out);
              

		}}

#printing s-matrix
print "S-Matrix for Nussinov algorithm\n\n";
foreach $r(1..length $sequence){
		foreach $s(1..$r){
			
			printf ("%3d",$smatrix[$r][$s]);}
		print"\n\n";}
		