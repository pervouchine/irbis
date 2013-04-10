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

# These scripts are to be documented ASAP
#
#
    $H1 = 2;        #half exon height
    $H2 = 5;        #arrow level
    $H3 = 9;        #top of the splicing arc
    $H5 = 25;       #distance between diagrams
    $TEXWIDTH = 500; # width of TEX picture environment
    $MAXSEGM  = 500; # maximal length of segment in bp before truncation
    @COLORS = ( 'red', 'green', 'blue', 'magenta', 'brown', 'olive', 'orange', 'cyan', 'purple', 'violet','BurntOrange','RoyalPurple','RoyalBlue','BrickRed'); 
		#'Plum', 'MidnightBlue','BurntOrange','Fuchsia','RoyalPurple','RoyalBlue','BrickRed','BlueViolet');

return(1);

sub drawgene {
    my $gene  = shift(@_);
    my $level = shift(@_);
    my $label = shift(@_);
    my %pos = @_;
    my @sites = sort {$a<=>$b} @{$SITES{$gene}};
    $label=~s/\_/-/g;

    foreach $site(@sites) {
        $pos{$site} = $POSITION{$site};
        $str = $STRAND{$site};
    }
    my @array = sort{$pos{$a}<=>$pos{$b}} keys(%pos);
    @array = reverse @array if($str<0);

    my %map = ();
    my $s=0;
    for($i=0;$i<@array-1;$i++) {
        $map{$array[$i]} = $s;
        $l = abs($pos{$array[$i+1]}-$pos{$array[$i]});
        $s+= (($l>$MAXSEGM) ? $MAXSEGM : $l);
    }
    $map{$array[$i]} = $s;

    my $output=undef;
    for($i=0;$i<@sites;$i++) {
        $x = round($map{$sites[$i]}*$TEXWIDTH/$s);
        $y = round($map{$sites[$i+1]}*$TEXWIDTH/$s);
        if($SEGM{$sites[$i]}=~/^[IX].$/) {
            $output.="\\draw ($x,$level) -- ($y, $level);\n";
        }
        if($SEGM{$sites[$i]}=~/^[EA].$/) {
            $output.= ($SEGM{$sites[$i]}=~/^.C$/ ? "\\draw[fill=gray] " : "\\draw ")." ($x,".($level-$H1).") rectangle ($y,".($level+$H1).");\n";
        }
        $output.= "\\draw[->] ($x,$level) |- ++(3,$H2);\n"          if($TYPE{$sites[$i]} eq "S");
        $output.= "\\draw ($x,$level) -- ++(0,".($H2-1).") ++(0,1) circle (1);\n" if($TYPE{$sites[$i]} eq "E");
    }

    foreach $intron(@{$INTRONS{$gene}}) {
        next unless($map{$intron->[0]} && $map{$intron->[1]});
        $x = round($map{$intron->[0]}*$TEXWIDTH/$s);
        $y = round($map{$intron->[1]}*$TEXWIDTH/$s);
        $z = round(($x + $y)/2);
        $h = round($H3 + log(1+abs($y-$x)));
        $output.="\\draw[thin] ($x,".($level+$H1).") .. controls ($z,".($level+$h).") .. ($y,".($level+$H1).");\n";
    }

    foreach $key(@_) {
        if($key=~/([lr][be]\d+)/) {
            $x = round($map{$1}*$TEXWIDTH/$s);
            @{$box_coordinates{$1}} = ($x, $level);
        }
    }
    $output.="\\draw($TEXWIDTH,$level) node[anchor=west]{\\small $label};\n";
    return($output);
}

sub drawboxes {
    my $out=undef;
    my $pair = @_[0];
    foreach $key(keys(%box_coordinates)) {
	if($key=~/lb(\d+)/) {
	    $b1 = ${$box_coordinates{"lb$1"}}[0];
            $e1 = ${$box_coordinates{"le$1"}}[0];
	    $h1 = ${$box_coordinates{"le$1"}}[1];
            $b2 = ${$box_coordinates{"rb$1"}}[0];
            $e2 = ${$box_coordinates{"re$1"}}[0];
            $h2 = ${$box_coordinates{"re$1"}}[1];
	    next unless($h1 ne undef && $h2 ne undef);
            $c  = $COLORS[$1 % @COLORS];
            $z1 = round(($b1 + $e1)/2);
            $z2 = round(($b2 + $e2)/2);
	    $z  = round(($z1 + $z2)/2);
	    $h = round(($h1 + $h2)/2);
	    $h-= round($H3 + log(1+abs($z1-$z2))) if($h1 eq $h2);
	    $head = "\\color{$c!95}\\draw[ultra thick] ($b1,$h1) -- ($e1, $h1);\\draw[ultra thick] ($b2,$h2) -- ($e2, $h2);";
            $out.= ($h1 eq $h2) ? "{$head \\draw ($z1,$h1) .. controls ($z,$h) .. ($z2,$h2);}" : "{$head \\draw ($z1,$h1) -- ($z2,$h2);}\n" if($1 eq $col{$pair->[0]}{$pair->[1]});
	    $out.= "{\\color{Yellow}\\draw[thick] ($z1,$h1) circle (3) ($z2,$h2) circle (4); }" if($1 eq $col{$pair->[0]}{$pair->[1]});
	}
    }
    return($out);
}

sub round {
    return(int(100*@_[0])/100);
}

###########################################################################################################


sub add_plots {
    my @a = split /\t/, $DATA{@_[0]};
    my @b = split /\t/, $DATA{@_[1]};
    my $out=undef;
    foreach $line(@PLOTS) {
        my @arr = split /\s+/, $line;
        my $color = 'blue';
        my $mark  = '*';
        my $scale = 0.71;
        my $xlab  = "x";
        my $ylab  = "y";
        my $ttl   = "title";
        my $x = $y = $xcol = $ycol = $labels = undef;
        foreach my $tag(@arr) {
            $x=$1       if($tag =~ /x=(.*)/);      # identifier of the x-dataset (for instance, GENE.rpkm)
            $y=$1       if($tag =~ /y=(.*)/);      # identifier of the y-dataset
            $xcol=$1-1  if($tag =~ /xcol=(\d*)/);  # column number for the x data
            $ycol=$1-1  if($tag =~ /ycol=(\d*)/);  # column number for the y data
            $xlab=$1    if($tag =~ /xlab=(.*)/);   # x-axis label 
            $ylab=$1    if($tag =~ /ylab=(.*)/);   # y-axis label
            $ttl=$1     if($tag =~ /main=(.*)/);   # graph title
            $color=$1   if($tag =~ /color=(.*)/);  # color of the points
            $mark=$1    if($tag =~ /mark=(.*)/);   # mark symbol
            $scale=$1   if($tag =~ /scale=(.*)/);  # scale of the graph
	    $labels=$1  if($tag =~ /labels=(.*)/); # labels yes/no
        }
        next unless($x && $y && $xcol && $ycol);
        %X = %{$VALUE{$x}{$a[$xcol]}};
        %Y = %{$VALUE{$y}{$b[$ycol]}};
        my $res = undef;
        foreach my $key(keys(%X)) {
            next unless($X{$key}=~/\d/ && $Y{$key}=~/\d/);
            $res.="($X{$key},$Y{$key})";
	    $res.=($labels eq "yes") ? "[$key]" : "[~]";
        }
        if($res) {
            $out.= "\\begin{tikzpicture}[scale=$scale]\n\\begin{axis}[nodes near coords, enlargelimits=0.2, grid=major, xlabel=$xlab, ylabel=$ylab, title=$ttl, ymin=0, ymax=1]\n";
            $out.= "\\addplot+[only marks,point meta=explicit symbolic,color=$color,mark=$mark,mark options={fill=$color}] coordinates {$res};\n\\end{axis}\n\\end{tikzpicture}\n";
        }
    }
    return($out);
}
