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

for(my $i=0;$i<@ARGV;$i++) {
    $cfg  = $ARGV[++$i] 	if($ARGV[$i] eq "-i");
    push @names, $ARGV[++$i] 	if($ARGV[$i] eq "-n");
}

die unless($cfg || @names==0);
read_configuration($cfg);


foreach $name(@names) {
    @arr = split /\n/, `awk '\$11=="$name"' $annotation_file_name | cut -f5`;
    foreach $site(@arr) {
	$id{$site}=$name;
    }
    push @sites, @arr;
}


$n=0;
foreach $z(keys(%index)) {
    print STDERR "$n $path$z$extention.ind\n";
    %a = split /[\t\n]/, `cat $path$z$extention.ind`;
    foreach $z(@sites) {
	$f{$z}{$n}=($a{$z}? "*" : " ");
    }
    $n++;
}

foreach $z(sort{$a<=>$b}keys(%f)) {
    print "$z\t";
    for($i=0;$i<$n;$i++) {
	print $f{$z}{$i};
    }
    print "\t$id{$z}\n";
}
