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

# This script filters out the lncRNA from ARGV[0] (Rory's file) which overlap with protein coding genes (either strand)

if(@ARGV==0) {
    print STDERR "This script filters out the lncRNA which overlap with protein coding genes (either strand)\n";
    print STDERR "Usage: $0 -s cps_file -i lncRNA_file (Rory's table) -o output_cps [-z overlapping_file]\n";
    exit(1);
}

for(my $i=0;$i<@ARGV;$i++) {
    $cps_file   = $ARGV[++$i] if($ARGV[$i] eq "-s");
    $inp_file   = $ARGV[++$i] if($ARGV[$i] eq "-i");
    $lnc_file   = $ARGV[++$i] if($ARGV[$i] eq "-z");
}
die("No database, exiting\n") unless($cps_file);
die("No input, exiting\n")    unless($inp_file);

open FILE, $cps_file || die("Can't open $cps_file, exiting\n");
print STDERR "[Reading $cps_file";
while($line=<FILE>) {
    chomp($line);
    ($chr, $pos, $str, $gene, $site, $type, $u, $d, $name, $bt) = split /\t/, $line;
    if($bt eq "protein_coding") {
	push @{$position{$name}}, $pos;
	$location{$chr}{$name}=1;
    }
    push @{$data{$name}}, $line;
}
print STDERR "]\n";
close FILE;

print STDERR "[Processing ";
foreach $name(keys(%position)) {
    @arr = sort {$a<=>$b} @{$position{$name}};
    $min{$name} = shift(@arr);
    $max{$name} = pop(@arr);
}
print STDERR 0+keys(%max),"]\n";

open FILE, $inp_file || die("Can't open $inp_file, exiting\n");
print STDERR "[$inp_file -> $out_file";
while($line=<FILE>) {
    chomp($line);
    ($ensg, $blah, $loc) = split /\s+/, $line;
    ($chr,  $pos1, $pos2) = split /[\:\-]/, $loc;
    ($pos1, $pos2) = sort {$a<=>$b} ($pos1, $pos2);
    $flag = undef;
    foreach $name(keys(%{$location{$chr}})) {
	$flag=$name unless($max{$name} < $pos1 || $pos2 < $min{$name}); 
	last if($flag);
    }
    print STDERR $flag ? "!" : "";
    ($name, $blah) = split /\./, $ensg;
    unless($flag) {
	foreach $line(@{$data{$name}}) {
	    ($chr, $pos, $str, $gene, $site, $type, $u, $d, $name, $bt) = split /\t/, $line;
	    print "$line\n" if($pos1<=$pos+1 && $pos-1<=$pos2);
	}
    }
    else {
	$rest.="$line\t::$flag\n";
    }
}
print STDERR "]\n";
close FILE;

if($lnc_file) {
    open OUT, ">$lnc_file" || die("Can't open $lnc_file, exiting\n");
    print STDERR "[Saving rests to $lnc_file";
    print OUT $rest;
    print STDERR "]";
    close OUT;
}
