#!/usr/bin/perl

use Perl::utils;

@levels = (0.7, 0.74, 0.78, 0.82, 0.86, 0.9);
@length = (8, 9, 10, 11, 12);
$params = '-t 0.7 -g 1 -u 0';

($metadata, $output, $db, $config) = @ARGV;

read_configuration($config);

$dir  = "benchmark/files/";
$file = "test3_$db"; 

`perl Perl/all2all.pl -in $metadata/$db.cps -out $dir$db.a2a`;
`Progs/irbis -l $output/ncpcg.met -r $output/ncpcg.met -o $dir$file.tab $params -b $dir$db.a2a` unless(-e "$dir$file.tab");

test();

for($n=0; $n<=9; $n++) {
    unless(-e "$dir$file.tab$n") {
    	$line = `perl benchmark/shuffle.pl $dir$db.a2a $output/ncpcg.met.cns $metadata/ncpcg.cps 50 | sort -k1,1n -k2,2n > $dir$db.shuffled.srt`;
	$line = `Progs/irbis -l $output/ncpcg.met -r $output/ncpcg.met -o $dir$file.tab$n $params -b $dir$db.shuffled.srt`;
    }
    print STDERR "[$n]";
    test($n);
}


sub test {
    my $n = @_[0];
    my $m = $n + (($n eq undef) ? 0 : 1);
    foreach $T(@levels) {
	foreach $L(@length) {
	    open FILE, "$dir$file.tab$n";
	    %tot_weight=();
	    while($line=<FILE>) {
    		chomp($line);
    		($specie, $left_site, $right_site, $len, $beg1, $end1, $beg2, $end2, $s1, $s2) = split /\t/, $line;
    		next unless($len>=$L);
    		$tot_weight{$left_site}{$right_site} += $weigh{$species[$specie]};
	    }
	    $number = 0;
	    foreach $left_site(keys(%tot_weight)) {
		foreach $right_site(keys(%{$tot_weight{$left_site}})) {
	    	    $number++ if($tot_weight{$left_site}{$right_site}>=$T*$sumweigh);
		}
	    }
	    close FILE;
	    print "$m\t$T\t$L\t$number\n";
	}
    }
}
