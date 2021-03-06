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



if(@ARGV==0) {
    print STDERR "This script takes .out output of the aligner and creates segments in bed format\n";
    print STDERR "Usage: $0 -i [.out file] -o [bed file] -cis [report in cis]\n";
    exit(1);
}

for(my $i=0;$i<@ARGV;$i++) {
    $tab_file   = $ARGV[++$i] if($ARGV[$i] eq "-in");
    $bed_file   = $ARGV[++$i] if($ARGV[$i] eq "-out");
    $cis = 1 		      if($ARGV[$i] eq "-cis");
}
die("No input, exiting\n") unless($tab_file);
die("No output, exiting\n") unless($bed_file);

print STDERR "[<$tab_file";
@list = split /\n/, `cat $tab_file`;
foreach $line(@list) {
    if($cis) {
	($chr2, $pos2, $str2, $gene, $site, $type) = split /\t/,$line;
    }
    else {
    	($chr1, $pos1, $str1, $chr2, $pos2, $str2, $gene, $site, $type) = split /\t/,$line;
    }
    push @{$data{$gene}}, "$chr2\t$pos2\t$str2";
    push @{$site{$gene}}, $site;
}
print STDERR "]\n";

print STDERR "[>$bed_file";
open FILE, ">$bed_file";
@keys = sort {$a<=>$b} (keys(%data));
foreach $gene(@keys) {
    for($i=0;$i<@{$data{$gene}}-1;$i++) {
	($chr1, $pos1, $str1) = split /\t/, ${$data{$gene}}[$i];
        ($chr2, $pos2, $str2) = split /\t/, ${$data{$gene}}[$i+1];
	unless($chr1 eq $chr2 && $str1 eq $str2) {
	    $problem++;
	}
	$site1 = ${$site{$gene}}[$i];
	$site2 = ${$site{$gene}}[$i+1];
	($site1, $site2) = sort {$a<=>$b} ($site1, $site2);
	($pos1, $pos2) = sort {$a<=>$b} ($pos1, $pos2);
	next unless($site2-$site1==1);
	print FILE "$chr1\t$pos1\t$pos2\t.\t$site1\t$str1\n" unless($been{$site1});
	print STDERR "!" if $been{$site1};
	$been{$site1}=1;
    }
}
close FILE;
print STDERR "]\n";
print STDERR "[WARNING: there were $problem segments with incomartible chromosome or strand]\n" if($problem);

