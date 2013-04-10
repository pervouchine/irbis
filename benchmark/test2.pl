#!/usr/bin/perl

use Perl::utils;
use Perl::align;

$n=10;                          # number of species
$m=100;                         # number of segments
$M=500;                         # segment length
$k=8;                           # seed length
$INF=65535*65535;               # infinity
$th=0.7;                        # threshold t_2
$dir  = "benchmark/files/";
$file = "test2";
$N_iter = 10;

my $percent = 0.02;

my @species=();
my @seq =();
my %h = ();

$wc{A}{T}=$wc{T}{A}=$wc{C}{G}=$wc{G}{C}=1;

$N=0;

for(my $j=1;$j<=$m;$j++) {
    print STDERR "[$j]";

    @species=();
    open FILE, ">$dir$file$$.fa";
    $seq = randseq($M);
    evolve($seq);
    for(my $i=0;$i<@species;$i++) {
	$seq[$i][$j]=$species[$i];
	print FILE ">$i\n$species[$i]\n";
    }
    close FILE;
    system("muscle -in $dir$file$$.fa -out $dir$file$$.aln -quiet");
    my %f = fasta2hash("$dir$file$$.aln");
    system("rm -f $dir$file$$.fa $dir$file$$.aln");

    open FILE, ">$dir$file.clw";
    print FILE "CLUSTAL W (1.81) multiple sequence alignment\n\n";
    foreach my $z(sort{$a<=>$b} keys(%f)) {
    	print FILE "$z\t$f{$z}\n";
	push @{$output{$z}}, "$j\tchr\t".($j*$M)."\t$M\t1\t$INF\t$f{$z}";
    }
    close FILE;

    ($trash, $out) = split /\s/, `benchmark/ViennaRNA-1.8.5/Progs/RNAalifold $dir$file.clw  2>>/dev/null`;
    @stack = ();
    @res = ();

    for($i=0;$i<length($out);$i++) {
	push @stack, $i if(substr($out,$i,1) eq "(");
	if(substr($out,$i,1) eq ")") {
	    $q = pop(@stack);
	    push @res, [$q,$i];
	}	
    }

    @res = sort{$a->[0] <=> $b->[0]} @res;
    $QQ = 0;
    for($i=0;$i<@res-$k;$i++) {
	$flag = 1;
	for($q=0;$q<$k-1;$q++) {
	    $flag = undef unless($res[$i+$q+1]->[0] - $res[$i+$q]->[0] == 1 && $res[$i+$q+1]->[1] - $res[$i+$q]->[1] == -1); 
	}
	if($flag) {
	    $RR = 0;
	    foreach my $z(sort{$a<=>$b} keys(%f)) {
		$s = $f{$z};
		$c = 0;
		for($q=0;$q<$k;$q++) {
		    $c+=$wc{substr($s,$res[$i+$q]->[0],1)}{substr($s,$res[$i+$q]->[1],1)};
		}
		$RR++ if($c>=$k-1);
	    }
	    $QQ=1 if($RR>=keys(%f)*$th);
	}
    }
    $N++ if($QQ>0);
}

foreach $z(keys(%output)) {
    open FILE, ">$dir$file\_$z.sus";
    print FILE join("\n", @{$output{$z}});
    close FILE;
}

open FILE, ">$dir$file.cps";
for(my $j=1;$j<=$m;$j++) {
    print FILE "chr\t",$j*$M,"\t1\t1\t",$j,"\n";
}
close FILE;

open FILE, ">$dir$file.cfg";
print FILE "halfsize\t",$k/2,"\ngap\t0\nthreshold\t$th\ngccontent\t0\n\npath $dir\n";
for(my $i=0;$i<$n;$i++) {
    print  FILE "species $file\_$i 1\n";
}
print FILE "extention\t.sus\nannotation\t$file.cps\n";
close FILE;


open FILE, ">$dir$file.a2a";
for(my $j=1;$j<=$m;$j++) {
    print FILE "$j\t$j\n";
}
close FILE;

system("Progs/trim  -i $dir$file.cfg -o $dir$file.met -v");
($Z) = split /\s/, `Progs/irbis -l $dir$file.met -r $dir$file.met -u 1 -g 1 -o $dir$file.out -b $dir$file.a2a -v`;
$W = `cut -f2,3 $dir$file.out | uniq | wc -l`;
print "$N\t$Z\t$W\n";
exit;


sub evolve {
    my $level = @_[1];
    if($level==4) {
	push @species, @_[0];
    }
    else {
    	evolve(mutate(@_[0]), $level+1);
    	evolve(mutate(@_[0]), $level+1);
    }
}

sub mutate {
    my $x = @_[0];
    my @arr = split //, "acgt";
    for(my $i=0;$i<$percent*length($x);$i++) {
	my $s = $arr[int(rand()*@arr)];
	if($s eq "-") {
	    my $len = int(rand()*5);
	    my $pos = int(rand()*(length($x)-$len));
	    substr($x,$pos,$len) = "";
	
	}
	else {
	    my $pos = int(rand()*length($x));
	    substr($x,$pos,1) = $s;
	}
    }
    return($x);
}
