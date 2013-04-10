#!/usr/bin/perl

use Perl::utils;

$n=10;				# number of species
$m=100;				# number of segments
$M=1000;			# segment length
$k=8;				# seed length
$INF=65535*65535;		# infinity
$th=0.99;			# threshold t_2
$dir  = "benchmark/files/";
$file = "test1";
$N_iter = 10;

open FILE, ">$dir$file.cps";
for($j=1;$j<=$m;$j++) {
    print FILE "chr\t",$j*$M,"\t1\t1\t",$j,"\n";
}   
for($j=1;$j<=$m;$j++) {
    print FILE "chr\t",$j*$M,"\t1\t1\t",$j+$m,"\n";
}   
close FILE;

open FILE, ">$dir$file.cfg";
print FILE "halfsize\t",$k/2,"\ngap\t2\nthreshold\t$th\ngccontent\t0\n\npath $dir\n";
for($i=0;$i<$n;$i++) {
    print  FILE "species $file\_$i 1\n";
}
print FILE "extention\t.sus\nannotation\t$file.cps\n";
close FILE;


for($I=0;$I<$N_iter;$I++) {
    print STDERR "Iteration $I: ";
    for($j=1;$j<=$m;$j++) {
    	$mem[$j] = randseq($k);
    	$mem[$j] =~ tr/[a-z]/[A-Z]/;
    }
    for($i=0;$i<$n;$i++) {
    	open FILE, ">$dir$file\_$i.sus";
    	for($j=1;$j<=$m;$j++) {
	    $seq  = randseq($M);
	    $pos = int(rand()*($M-$k-1));
	    substr($seq, $pos, $k) = $mem[$j];
	    print FILE "$j\tchr\t",$j*$M,"\t$M\t1\t$INF\t$seq\n";
    	}
    	for($j=1;$j<=$m;$j++) {
            $seq  = randseq($M);
	    $pos = int(rand()*($M-$k-1));
            substr($seq, $pos, $k) = revcomp($mem[$j]) if($i<=$th*$n);
            print FILE $m+$j,"\tchr\t",($j+$m)*$M,"\t$M\t1\t$INF\t$seq\n";
    	}
        close FILE;
    }
    `Progs/trim -i $dir$file.cfg -o $dir$file.met -v`;
    ($num) = split / /, `Progs/irbis -l $dir$file.met -r $dir$file.met -u 1 -o $dir$file.out -v`;
    print STDERR "$num\n";
    print "$num\n";
}

sub spoil {
    my $x=substr(@_[0],0,$k/2);
    my $y=substr(@_[0],$k/2);
    return($x.randseq(int(rand()*2)).$y);
}


