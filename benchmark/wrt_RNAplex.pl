@alignments = split /\n\n\n/, `cat $ARGV[0]`;

foreach $aln(@alignments) {
    ($left, $right) = split /\n\n/, $aln;
    $x = out($left,  "query.aln");
    $y = out($right, "target.aln");
    @array = split /\s+/, `benchmark/RNAplex-0.2/Progs/RNAplex -A query.aln target.aln 2>>/dev/null`;
    my %a = plex_parse(@array);
    my %b = irbis_parse($x, $y);
    print join("\t", compare(\%a,\%b)), "\n";

}

sub out {
    my @array = split /\n/, @_[0];
    shift(@array);
    open FILE, ">@_[1]";
    print FILE "CLUSTAL\n\n";
    my $res= undef;
    foreach $line(@array) {
	($s, $chr, $pos, $beg, $str, $len, $seq) = split /\s+/, $line;
	$res = $seq unless($res);
	$z = (' ' x 40).$seq;
	substr($z, 0, length($chr)) = $chr;
	print FILE "$z\n";
    }
    close FILE;
    return($res);
}


sub plex_parse {
    my @array = @_;
    my %res = ();
    ($b1,$e1) = split /\,/, $array[1];
    ($b2,$e2) = split /\,/, $array[3];
    ($l, $r)  = split /\&/, $array[0];
    while($l && $r) {
	while($l=~s/^\.//) {$b1++;}
	while($r=~s/\.$//) {$e2--;}
	if($l=~s/^\(// && $r=~s/\)$//) {
	    #print "$b1 $e2\t",substr($x,$b1-1,1), " ", substr($y,$e2-1,1),"\n";
	    $res{"$b1;$e2"}++;
	    $b1++;$e2--;
	}
    }
    return(%res);
}

sub irbis_parse {
    my $x = @_[0];
    my $y = @_[1];
    my $i = 1;
    my $j = length($y);
    my %res = ();
    while($x && $y) {
	while($x=~s/^[acgt\-\d\.\:]//) {$i++;}
	while($y=~s/[acgt\-\d\.\:]$//) {$j--;}
	if($x=~s/^[A-Z]// && $y=~s/[A-Z]$//) {
	    #print "$i $j\t",substr(@_[0],$i-1,1), " ", substr(@_[1],$j-1,1),"\n";
	    $res{"$i;$j"}++;
	    $i++;
	    $j--;
	}
    }
    return(%res);
}



sub compare {
    my %a = %{@_[0]};
    my %b = %{@_[1]};
    my $n=0;
    foreach $z(keys(%a)) {
	$n++ if($b{$z});
    }
    return(0+keys(%a), 0+keys(%b), $n);
}


