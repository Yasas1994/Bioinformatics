#Needleman and Wunsch algorithm
#input sequence 1&2,missmatch and match score 
require List::Util;
open F,"BLOSUM62.txt" or die "File not found";
print "Needleman and Wunsch\n\n";
print"input sequence 1:"; chomp($seq1=<>);#row
printf"\ninput sequence 2:"; chomp($seq2=<>);#col
print"\nmissmatch score:";chomp($miss=<>);
print"(M) match score:";chomp($mat=<>);
print"(G) gap penalty:";chomp($gap=<>);print"\n";
$seq1=uc $seq1;
$seq2=uc $seq2;
@seq1=split//,$seq1;
@seq2=split//,$seq2;
 unshift @seq1,0;
 unshift @seq2,0;
#print @seq1;
#similarity matrix


foreach $i(1..$#seq1){
	foreach $j(1..$#seq2){
		if ($seq1[$i] eq $seq2[$j])
			{$matsim[$i][$j]=$mat;}
				else{$matsim[$i][$j]=$miss;}
		
		}}
=head		
#print similarity matrix	
foreach $r(1..$#seq1){
	foreach $s(1..$#seq2){
		
		printf ("%3d",$matsim[$r][$s]);}
		print"\n";}
=cut

#fmatrix; fills up first row and first column 
foreach $xi(0..$#seq1){
	$fmat[$xi][0]=$gap*$xi;
	$xmat[$xi][0]=$gap*$xi;}#rows
foreach $yj(0..$#seq2){
	$fmat[0][$yj]=$gap*$yj;
	$xmat[0][$yj]=$gap*$yj;}#col
	print"\n\n";

#fills up the matrix based on matrix formular 

foreach $I(1..$#seq1){
	foreach $J(1..$#seq2){
	   if($seq1[$I] eq $seq2[$J]){
	   	$score=$mat;
	        $match=$fmat[($I)-1][($J)-1]+$score;#matrix
	        }
		   else{
			$score=$miss;
			$match=$fmat[($I)-1][($J)-1]+$score;#matrix
			}
		
		$gap1=$fmat[$I-1][$J]+$gap;#adds gaps in seq2
		$gap2=$fmat[$I][$J-1]+$gap;#adds gaps in seq1
		$fmat[$I][$J]=List::Util::max($match,$gap1,$gap2);
			
			if($fmat[$I][$J]==$match)
			    {$xmat[$I][$J]=0;
			    #print "$xmat[$I][$J]\n";
			    }
				if($fmat[$I][$J]==$gap1)
				  {$xmat[$I][$J]=1;
				  #print "$xmat[$I][$J]\n";
				  }
					if($fmat[$I][$J]==$gap2)
					   {$xmat[$I][$J]=-1;
					   #print "$xmat[$I][$J]\n";
					   }
						}}
=head
#print fmatrix	
foreach $r(0..$#seq1){
	foreach $s(0..$#seq2){
		
		printf ("%3d",$fmat[$r][$s]);}
		print"\n\n";}
		print"\n\n";

#print xmatrix(traceback)
foreach $r(0..$#seq1){
	foreach $s(0..$#seq2){
		
		printf ("%3d",$xmat[$r][$s]);}
		print"\n\n";}
=cut		

#back-trasing to find the alignment 
#for ($IB=scalar @seq1;$IB>=0;$IB--){
	#for ($IB=$#seq1;$IB>=0;$IB--){

$m=$#seq2;
$n=$#seq1;
#extends the aligned sequences usinf the traceback matrix	
		while($m>0 or $n>0){
			#print"$xmat[$n][$m]:$n N $m\n ";
			$sc=$sc+$fmat[$n+1][$m+1];
			if($xmat[$n][$m]==0){
				$ali1=$ali1.$seq2[$m];
				$ali2=$ali2.$seq1[$n];
				$n=$n-1;
				$m=$m-1;
					if($m<-1 or $n<-1){
					  print "Sequence extension terminated.Please tweak the parameters!\n";last;}}
				
				elsif($xmat[$n][$m]==1){
					$ali2=$ali2.$seq1[$n];
					$ali1=$ali1."-";
					$n=$n-1;
						if($m<-1 or $n<-1)
						{print "Sequence extension terminated.Please tweak the parameters!\n";last;}}
					
					elsif($xmat[$n][$m]==-1) {
						$ali2=$ali2."-";
						$ali1=$ali1.$seq2[$m];
						$m=$m-1;
							if($m<-1 or $n<-1)
							{print "Sequence extension terminated.Please tweak the parameters!\n";last;}}
				}
			
			
					
#calculates the Alignment score using BLOSUM62 matrix
				
#outputs aligned sequences and score
$d= $#seq1;
$c=$#seq2;
			$a1=reverse $ali1; $a2= reverse $ali2;
			@a11=split//,$a1;
			@a22=split//,$a2;
			foreach $abc(0..$#a11){
						if ($a11[$abc] eq $a22[$abc]){
						$ali3=$ali3.'*';}
							else{$ali3=$ali3." ";}
							}
					
		        print"Sequence 1:$a1\nSequence 2:$a2\n          :$ali3\n\n\nAlignment score:$fmat[$d][$c]\n";
			
	
#this block calculates alignment score using BLOSUM62 matrix

#takes input for processing
#print"input 1:"; chomp($in1=<>);
#print"input 2:"; chomp($in2=<>);
#print"\n";
#print"Gap penalty:";chomp($gap=<>);
#@in1=split//,$in1;
#@in2=split//,$in2;

$iii=0;
$jjj=0;
while(<F>){

	if($_=~/^[A-Z\*]\s{1,2}(.+)$/){
		
		@line=split/\s{1,2}/,$1;
		foreach $jjj(0..23){
			$blo[$iii][$jjj]=@line[$jjj];}
			$iii++;}
			         }
=head
			foreach $ii(0..23){
				foreach $jj(0..23){
					printf("%3s",$blo[$ii][$jj]);}
					print "\n"; }
	
=cut
$n=0;		
%ami=("A",0,"R",1,"N",2,"D",3,"C",4,"Q",5,"E",6,"G",7,"H",8,"I",9,"L",10,"K",11,"M",12,"F",13,"P",14,"S",15,"T",16,"W",17,"Y",18,"V",19,"B",20,"Z",21,"X",22,"*",23);			
	
	
 #calculates alignment score

       if ($#a11==$#a22){
       	foreach $p(0..$#a22){
       ;
       			if(exists $ami{$a11[$p]} and exists $ami{$a22[$p]}){
       				$sco=$sco + ($blo[$ami{$a11[$p]}][$ami{$a22[$p]}]);
       				}
			else{$sco+=$gap;}
			
			}
		
       	}
       
       		 print"score(BLOSUM 62):$sco\n\n";
			 
