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

if(@ARGV==0) {
    print STDERR "Takes CPS input and creates binary file with relation all-2-all sites within gene\n";
    print STDERR " -in  input cps\n";
    print STDERR " -out output binary file\n";
    print STDERR " -bin binary mode\n";
    exit(0);
}

for(my $i=0;$i<@ARGV;$i++) {
    $cps     = $ARGV[++$i] if($ARGV[$i] eq "-in");
    $outfile = $ARGV[++$i] if($ARGV[$i] eq "-out");
    $binmode = 1           if($ARGV[$i] eq "-bin");
}

die("No input provided, exiting\n") unless($cps);
die("No output, exiting\n") unless($outfile);

print STDERR "[Reading $cps";
open FILE, $cps;
while($line=<FILE>) {
    ($chr, $pos, $str, $gene, $site) = split /\t/, $line;
    $list{$gene}.="$site\t";
}
print STDERR "]\n";


open OUT, ">$outfile";
binmode(OUT) if($binmode);
@g = sort {$a<=>$b} keys(%list);
foreach $gene(@g) {
    progressbar(++$n,0+@g,"All-2-all ");
    @arr = sort {$a<=>$b} split /\t/, $list{$gene};
    for($i=0;$i<@arr;$i++) {
	for($j=$i;$j<@arr;$j++) {
	    if($binmode) {
		print OUT pack("L", $arr[$i]);
		print OUT pack("L", $arr[$j]);
	    }
	    else {
	    	print OUT "$arr[$i]\t$arr[$j]\n";
	    }
	}
    }
}
close OUT;
