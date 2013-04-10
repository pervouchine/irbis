#!/usr/bin/perl
# Binary relations library
# Binary relation = associative array with two keys; values don't matter

#       Copyright 2011,2012 Dmitri Pervouchine (dp@crg.eu)
#       This file is a part of the IRBIS package.
#       IRBIS package is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       IRBIS package is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with IRBIS package.  If not, see <http://www.gnu.org/licenses/>.
#

return(1);

##############################################################################################################
# I/O functions
##############################################################################################################

sub print_relation {
    my $x;
    my $y;
    foreach $x (keys(%{@_[0]})) {
	foreach $y (keys(%{${@_[0]}{$x}})) {
	    print "$x\t$y\n";
	}
    }
    print "\n";
}

##############################################################################################################
# Algebraic functions
##############################################################################################################

sub symmetrify {
    my $x;
    my $y;
    my %aux = %{@_[0]};
    foreach $x (keys(%aux)) {
        foreach $y (keys(%{$aux{$x}})) {
	    ${@_[0]}{$y}{$x}=${@_[0]}{$x}{$y};
	}
	${@_[0]}{$x}{$x}=1
    }
}

sub product {
    my %res=();
    foreach my $x (keys(%{@_[0]})) {
	foreach my $y (keys(%{${@_[0]}{$x}})) {
	    foreach my $z (keys(%{${@_[1]}{$y}})) {
		next unless($z);
		$res{$x}{$z}=1;
	    }
	}
    }
    return(%res);
}

sub adjoint {
    foreach my $x (keys(%{@_[0]})) {
	foreach my $y (keys(%{${@_[0]}{$x}})) {
	    ${@_[1]}{$x}{$y}=1;
	}
    }
}

##############################################################################################################
# Transitive closure 
##############################################################################################################

sub getall {
    my $x = @_[1];
    my $y;
    my @res = ($x);
    $flag{$x}=1;
    foreach $y (keys(%{${@_[0]}{$x}})) {
	next unless($y);
	push @res, getall(@_[0],$y) unless($flag{$y});
    }
    return(@res);
}

sub transitive_closure {
    %flag=();
    my $x;
    my %res=();
    foreach $x (keys(%{@_[0]})) {
	next unless($x);
	unless($flag{$x}) {
	    foreach $y(getall(@_[0],$x)) {
		next unless($y);
		$res{$x}{$y}=1;
	    }
	}
    }
    return(%res);
}

sub cliques {
    %flag=();
    my $x;
    my %res=();
    my $n=1;
    foreach $x (keys(%{@_[0]})) {
        next unless($x);
	unless($flag{$x}) {
	    my @q = getall(@_[0],$x);
	    my @s = ();
	    foreach $y(@q) {
                next unless($y);
		push @s, keys(%{$transcript_gene{$y}});
     	    }
	    @s = sort(@s);
	    $id = shift(@s);
	    foreach $y(@q) {
		next unless($y);
		$res{$y} = $id;
	    }
	    $n++;
	}
    }
    return(%res);
}
