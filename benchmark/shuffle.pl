@ids = split /\n/, `cut -f5 $ARGV[2]`;
@shuffle{@ids} = @ids;
@gene{@ids} = split /\n/, `cut -f4 $ARGV[2]`;
%conservation = split /[\t\n]/, `cat $ARGV[1]`;

foreach $id(@ids) {
    $bin = int($conservation{$id}*9999/10000/1000);
    push @{$arr{$bin}}, $id;
}

foreach $bin(sort {$a<=>$b} keys(%arr)) {
    for($i=0;$i<@{$arr{$bin}};$i++) {
	$j = int(rand()*@{$arr{$bin}});
	$x = $shuffle{$arr{$bin}->[$i]};
	$shuffle{$arr{$bin}->[$i]} = $shuffle{$arr{$bin}->[$j]};
	$shuffle{$arr{$bin}->[$j]} = $x;
    }
}

open FILE, $ARGV[0];
while($line=<FILE>) {
    chomp $line;
    ($x, $y) = split /\t/, $line;
    next unless($gene{$x} && $gene{$y});
    print "$x\t$shuffle{$y}\n";
    $k++ if($gene{$x} ne $gene{$shuffle{$y}});
    $n++;
}

print STDERR "[remaining = ", $k/$n,"]\n";

