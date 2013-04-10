#      Copyright 2011,2012 Dmitri Pervouchine (dp@crg.eu)
#      This file is a part of the IRBIS package.
#      IRBIS package is free software: you can redistribute it and/or modify
#      it under the terms of the GNU General Public License as published by
#      the Free Software Foundation, either version 3 of the License, or
#      (at your option) any later version.
#      
#      IRBIS package is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.
#      
#      You should have received a copy of the GNU General Public License
#      along with IRBIS package.  If not, see <http:#www.gnu.org/licenses/>.


for(my $i=0;$i<@ARGV;$i++) {
    $aln = $ARGV[++$i] if($ARGV[$i] eq "-in");
    $net = $ARGV[++$i] if($ARGV[$i] eq "-net");
}

die("[ERROR: cannot access input file]") unless($aln);
die("[ERROR: cannot access net file]") unless($net);

print STDERR "[<$aln";
open FILE, $aln || die("[ERROR: cannot access input file]");
while($line=<FILE>) {
    chomp $line;
    ($chr1, $pos1, $str1, $chr2, $pos2, $str2, $chain) = split /\t/, $line;
    $data{$chr1}{$pos1}{$chain} = $line;
}
print STDERR "]\n";

print STDERR "[<$net";
open FILE, $net || die("[ERROR: cannot access input file]");
while($line=<FILE>) {
    $chr = $1 if($line=~/^net (\w+)/);
    if($line=~/fill (\d+) (\d+)/) {
	$pos = $1;
	$len = $2;	
	next unless($chr && $pos && $len);
	if($line=~/id (\d+)/) {
	    $chain = $1;
	    foreach $x(keys(%{$data{$chr}})) {
		print "$data{$chr}{$x}{$chain}\n" if($pos<=$x && $x< $pos+$len && $data{$chr}{$x}{$chain});
	    }
	}
    }
}
