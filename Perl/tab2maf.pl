#!/usr/bin/perl
#	Copyright 2011,2012 Dmitri Pervouchine (dp@crg.eu)
#	This file is a part of the IRBIS package.
#	IRBIS package is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#	
#	IRBIS package is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	
#	You should have received a copy of the GNU General Public License
#	along with IRBIS package.  If not, see <http://www.gnu.org/licenses/>.

use Perl::utils;
use Perl::align;

$extention_program = "Progs/extend ";
$MAXALIGNMENTS = 1800;

if(@ARGV==0) {
    print STDERR "Creates tabular summary\n";
    print STDERR "Usage: $0 -l [left config] -r [right config] -i [tab output] -o [output file]\n";
    exit(1);
}

for(my $i=0;$i<@ARGV;$i++) {
    $rc = 1 if($ARGV[$i] eq "-rc");
    $left_file  = $ARGV[++$i] if($ARGV[$i] eq "-l");
    $right_file = $ARGV[++$i] if($ARGV[$i] eq "-r");
    $tab_file   = $ARGV[++$i] if($ARGV[$i] eq "-i");
    $out_file   = $ARGV[++$i] if($ARGV[$i] eq "-o");
    $cons_thr   = $ARGV[++$i] if($ARGV[$i] eq "-C");
    $maf_file_name = $ARGV[++$i] if($ARGV[$i] eq "-maf");
    $MAXALIGNMENTS = $ARGV[++$i] if($ARGV[$i] eq "-max");
}

die("Out file not specified, exiting\n") unless($out_file);

#########################################################################################################

print STDERR "[Tab-input:\t$tab_file";
open FILE, $tab_file || die(" can't open 'tab_file'\n");
%left=%right=%site_pairs=();
if($cons_thr) {
    while($line=<FILE>) {
    	chomp($line);
    	($specie, $left_site, $right_site, $len, $beg1, $end1, $beg2, $end2, $s1, $s2) = split /\t/, $line;
    	$cons_score{$left_site}{$right_site}+=$s1+$s2;
    }
    seek FILE,0,0;
    print STDERR ", pass2";
}
while($line=<FILE>) {
    chomp($line);
    ($specie, $left_site, $right_site, $len, $beg1, $end1, $beg2, $end2, $s1, $s2) = split /\t/, $line;
    next unless($cons_score{$left_site}{$right_site}>=$cons_thr);
    push @{$site_pairs{$left_site}{$right_site}}, $line;
    $left{$left_site}=$right{$right_site}=1;
}
close FILE;
print STDERR "]\n";

$number = 0;
foreach $left_site(keys(%site_pairs)) {
    $number+=keys(%{$site_pairs{$left_site}});
}
print STDERR "Conservation threshold = $cons_thr\n" if($cons_thr);
print STDERR "Found $number site pairs\n";

unless($number>0) {
    print STDERR "Can't proceed, exiting without $out_file\n";
    exit(0);
}

if($number>$MAXALIGNMENTS) {
    print STDERR "[WARNING: will only do first $MAXALIGNMENTS alignments]\n";
    $number=$MAXALIGNMENTS;
}

read_configuration($left_file);
read_annotation();
read_introns();
read_sequences(%left);
@LEFT_DATA = @data;
@LEFT_ATTR = @attr;

read_configuration($right_file);
read_annotation();
read_introns();
read_sequences(%right);
@RIGHT_DATA = @data;
@RIGHT_ATTR = @attr;

#########################################################################################################
open OUTPUT, ">$out_file" || die("Output file can't be opened, exiting\n");
print STDERR "[Save output to $out_file]\n";

foreach $left_site(keys(%site_pairs)) {
    foreach $right_site(keys(%{$site_pairs{$left_site}})) {
	last if($counter>$number);
	progressbar(++$counter,$number,"Alignment ");
	%lseq = %rseq = ();
	$ld = $rd = undef;
	foreach $line(@{$site_pairs{$left_site}{$right_site}}) {
	    ($org, $x, $y, $l, $p1, $p2, $q1, $q2, $lc, $rc, $ls, $rs) = split /\t/, $line;
	    $loopsize = 2*$gapsize;
	    ($lseq{$org}, $rseq{$org}) = split /\n/, `$extention_program -i $LEFT_DATA[$org][$x] $p1 $p2 $RIGHT_DATA[$org][$y] $q1 $q2 -helix $halfsize -loop $loopsize`;
	    #write2log("$left_site:$right_site\t$extention_program -i $LEFT_DATA[$org][$x] $p1 $p2 $RIGHT_DATA[$org][$y] $q1 $q2 -helix $halfsize -loop $loopsize");
	    $ld = $p1 if($org == 0);
            $rd = $q2 if($org == 0);
	}
    	%lseq = qalign(%lseq);
    	%rseq = qalign(%rseq);

        print OUTPUT "a id=$x\n";
        for($i=0;$i<$nsp;$i++) {
            print OUTPUT $lseq{$i} ? "s $LEFT_ATTR[$i][$x] $lseq{$i}\n" : "";
    	}
    	print OUTPUT "\n";
    	print OUTPUT "a id=$y\n";
    	for($i=0;$i<$nsp;$i++) {
            print OUTPUT $rseq{$i} ? "s $RIGHT_ATTR[$i][$y] $rseq{$i}\n" : "";
    	}
    	print OUTPUT "\n\n";
    }
}
close OUTPUT;
print STDERR "[WARNING: $nwarnings sequence trimmed to $MAXLENMUSCLE]\n" if($nwarnings);
print STDERR "Done\n";
exit;
