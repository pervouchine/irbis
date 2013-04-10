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
use Perl::annot;

$spacer = "\\\\\\vspace\{0.5cm\}";

if(@ARGV==0) {
    print STDERR "No args\n";
    exit(1);
}

for(my $i=0;$i<@ARGV;$i++) {
    $title = 1 		 if($ARGV[$i] eq "-title");
    $left_file     = $ARGV[++$i] if($ARGV[$i] eq "-l");
    $right_file    = $ARGV[++$i] if($ARGV[$i] eq "-r");
    push @tab_files, $ARGV[++$i] if($ARGV[$i] eq "-i");
    push @maf_files, $ARGV[++$i] if($ARGV[$i] eq "-m");
    $tex_file      = $ARGV[++$i] if($ARGV[$i] eq "-o");
    $expr_file	   = $ARGV[++$i] if($ARGV[$i] eq "-e");
    $cons_thr      = $ARGV[++$i] if($ARGV[$i] eq "-C");
    $binaryrel	   = $ARGV[++$i] if($ARGV[$i] eq "-b");
    $sources       = 1 if($ARGV[$i] eq "-s");

}

die("Output file not specified, exiting\n") unless($tex_file);
foreach $maf_file(@maf_files) {
    $flag = 1 if(-r $maf_file);
}
unless($flag) {
    print STDERR "Nothing to do, leaving\n";
    exit(0);
}

read_configuration($left_file);
read_annotation();
read_introns();
read_signatures();

read_configuration($right_file);
read_annotation();
read_introns();
read_signatures();

read_expression($expr_file);

if(-e $binaryrel) {
    print STDERR "[$binaryrel";
    open FILE, $binaryrel;
    while($line=<FILE>) {
        chomp($line);
	($l, $r, $label) = split /\t/, $line;
	$label =~ s/\_/\\\_/g;
	$BR{$l}{$r} = $label;
    }
    print STDERR "]";
}

##############################################################################################################

foreach $tab_file(@tab_files) {
    print STDERR "[Tab-input:\t$tab_file";
    open FILE, $tab_file || die(" can't open 'tab_file'\n");
    while($line=<FILE>) {
        chomp($line);
        ($specie, $left_site, $right_site, $len, $beg1, $end1, $beg2, $end2, $s1, $s2) = split /\t/, $line;
	$cons_score{$left_site}{$right_site}+=$s1+$s2;
    }
    seek FILE,0,0;
    print STDERR ", pass2";
    while($line=<FILE>) {
    	chomp($line);
    	($specie, $left_site, $right_site, $len, $beg1, $end1, $beg2, $end2, $s1, $s2) = split /\t/, $line;
        next if($cons_score{$left_site}{$right_site}<$cons_thr && $cons_thr);
    	push @{$site_pairs{$left_site}{$right_site}}, $line;
    }
    close FILE;
    print STDERR "]\n";
}

foreach $maf_file(@maf_files) {
    print STDERR "[Alignment:\t$maf_file ";
    open FILE, $maf_file || die(" - can't open '$maf_file'\n");
    my @array = split /\n\n\n/, `cat $maf_file`;
    foreach $pair(@array) {
    	($left, $right) = split /\n\n/, $pair;
    	($left_body,  $left_site)  = process($left);
    	($right_body, $right_site) = process($right);
        ($beg, $z, $z, $end) = sort{$a<=>$b} ($POSITION{$left_site}, $POSITION{$left_site+1}, $POSITION{$right_site}, $POSITION{$right_site+1});
        $coord = "$CHR{$left_site}:$beg-$end";
        $link = $CHR{$left_site} eq $CHR{$right_site} ? "\\href{http://genome.ucsc.edu/cgi-bin/hgTracks?db=$UCSC{$species[0]}&position=$coord}{VIEW}" : undef;
	$alignment{$left_site}{$right_site}="{\\footnotesize\\tt $left_body $spacer id=$left_site $DESC{$left_site} $spacer $right_body $spacer id=$right_site $DESC{$right_site} $link}\n";
    }
    print STDERR "]\n";
}

#die("No lines in input\n") unless(keys(%site_pairs)>0);

$number = 0;
foreach $left_site(keys(%site_pairs)) {
    foreach $right_site(keys(%{$site_pairs{$left_site}})) {
        $left_gene  = $GENE{$left_site};
        $right_gene = $GENE{$right_site};
        die("gene not found: ($left_gene $right_gene)") unless($left_gene && $right_gene);
        push @{$gene_pairs{$left_gene}{$right_gene}}, [$left_site, $right_site];
        $number++;
    }
}
print STDERR "Found $number site pairs\n";

$n=0;
$number = 0;
foreach $left_gene(keys(%gene_pairs)) {
    foreach $right_gene(keys(%{$gene_pairs{$left_gene}})) {
	print STDERR "Multiple match: $NAME{$left_gene} $NAME{$right_gene}\n" if(@{$gene_pairs{$left_gene}{$right_gene}}>1);
	foreach $pair(@{$gene_pairs{$left_gene}{$right_gene}}) {
	    foreach $line(@{$site_pairs{$pair->[0]}{$pair->[1]}}) {
		($specie, $left_site, $right_site, $len, $beg1, $end1, $beg2, $end2, $s1, $s2) = split /\t/, $line;
		$box{$left_gene}{"lb$n"}  = $POSITION{$left_site}  + $STRAND{$left_site}  * $beg1;
		$box{$left_gene}{"le$n"}  = $POSITION{$left_site}  + $STRAND{$left_site}  * $end1;
		$box{$right_gene}{"rb$n"} = $POSITION{$right_site} + $STRAND{$right_site} * $beg2;
		$box{$right_gene}{"re$n"} = $POSITION{$right_site} + $STRAND{$right_site} * $end2;
		$col{$left_site}{$right_site} = $n++;
		$ref{$left_gene}{$right_gene}.="(\\ref{$left_site:$right_site})" if($alignment{$left_site}{$right_site});
		last; 
	    }
	}
	$number++;
    }
}
print STDERR "Found $n box pairs\n";
print STDERR "Found $n gene pairs\n";
print STDERR "[WARNING: pgfplots suppressed]\n" unless($PGFPLOTS);

###############################################################################################################
my @document = ();

foreach $left_gene(sort{keys(%{$gene_pairs{$b}}) <=> keys(%{$gene_pairs{$a}})} keys(%gene_pairs)) {
    # draw title page
    $m = keys(%{$gene_pairs{$left_gene}});
    $m = 12 if($m>12);
    %box_coordinates=();
    $level = $H5*$m;
    $page = drawgene($left_gene, $level, $NAME{$left_gene}, %{$box{$left_gene}});
    foreach $right_gene(sort{@{$gene_pairs{$left_gene}{$b}} <=> @{$gene_pairs{$left_gene}{$a}}} keys(%{$gene_pairs{$left_gene}})) {
	next if($right_gene eq $left_gene);
	$level-=$H5;
	$page.=drawgene($right_gene, $level, $NAME{$right_gene}.$ref{$left_gene}{$right_gene}, %{$box{$right_gene}});
	last if($level<0);
    }
    $page.=drawboxes();
    $page ="\\begin{tikzpicture}[scale=0.05]\n$page\n\\end{tikzpicture}\n";
    push @document, "\\begin{figure}\n$page\\end{figure}\n\\clearpage\n\n" if($title);

    # draw additional pages
    foreach $right_gene(sort{@{$gene_pairs{$left_gene}{$b}} <=> @{$gene_pairs{$left_gene}{$a}}} keys(%{$gene_pairs{$left_gene}})) {
	%box_coordinates=();
	foreach $pair(@{$gene_pairs{$left_gene}{$right_gene}}) {
	    $page = drawgene($left_gene, $H5, $NAME{$left_gene}, %{$box{$left_gene}});
	    $page.= drawgene($right_gene, 0, $NAME{$right_gene}, %{$box{$right_gene}}) unless($right_gene eq $left_gene);
	    $page.= drawboxes($pair);
	    $page = "\\begin{tikzpicture}[scale=0.05]\n$page\n\\end{tikzpicture}\n";	
	    $text = $alignment{$pair->[0]}{$pair->[1]};
	    $text.= $spacer if($text);
	    $text.= add_plots($pair->[0],$pair->[1]) if($PGFPLOTS);
	    @cap = ();
	    push @cap, 'CS='.$cons_score{$pair->[0]}{$pair->[1]} if($cons_score{$pair->[0]}{$pair->[1]});
 	    push @cap, 'rel='.$BR{$pair->[0]}{$pair->[1]} if($BR{$pair->[0]}{$pair->[1]});
            $cap = join(",", @cap);
	    push @document, "\\begin{figure}\\centering $page$spacer$text\\caption{$cap}\\label{$pair->[0]:$pair->[1]}\\end{figure}\n\\clearpage\n\n" if($text);
        }
    }
}

$texfile = replace_extention($tex_file,".tex");
$texfile =~ s/\.tex$//i;
compile_tex($texfile,tex_document(join("",@document)),($title ? 2 : 1), $sources);
exit(0);



###############################################################################################################

sub process {
#this function works out ONE maf alignment block
    my @array = split /\n/, @_[0];
    my @data = @title = @link = ();
    my $id = undef;
    my $f = undef;
    my $g = undef;
    for $line(@array) {
	$id = $1 if($line =~ /^a id=(\d+)/);
	next unless($line =~ /^s/);
	($s, $chr, $pos, $len, $strand, $totlen, $seq) = split /\s+/, $line;
	$beg = ($strand eq "+") ? $pos : $totlen-$pos-$len;
	$end = ($strand eq "+") ? $pos+$len : $totlen-$pos;
	($name, $chr) = split /\./, $chr;
	$coord = "$chr:$beg-$end";
	$link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=$UCSC{$name}&position=$coord&Track=refGene";
	push @data, $seq;
	push @title, $sources ? $ABBR{$name} : "\\href\{$link\}\{$ABBR{$name}\}";
        ($f, $g) = features($seq,$id) if($name eq $species[0]); # splice sites as breaks
    }
    push @title,"";
    push @data, asterisks(@data);
    if($f) {
	my @tmp = split /\s/, $f;
	push @data, $g;
        for(my $i=0;$i<@data;$i++) {
	    for(my $s=0, my $j=0; $j<@tmp-1; $j++) {
	        $s+=length($tmp[$j]);
		die("$i $data[$i] $s $j $f\n") if($s+$j>length($data[$i]));
		substr($data[$i], $s+$j, 0)=" ";
	    }
    	}
	$g = pop @data;
    }
    return(do_tabular(\@title,\@data,$g), $id);
}

sub do_tabular {
# This function takes an array of sepcies's names, an array of (uncut) sequences, and calls tex_tabular to
# create a tex teable. It cuts sequences with the lentgh exceeds $TEXWINDOW
# Input:  1) ptr to an array of names, 2) ptr to an array with sequences
# Output: tex tabular environment 
    my @title = @{@_[0]};
    my @data  = @{@_[1]};
    my $bord  = @_[2];
    my $res = undef;
    my $window = int(length($data[0]) / (int(length($data[0])/$TEXWINDOW) + 1)) + 1;

    for(my $n = length($data[0]);$n>0;$n-=$window) {
	my @arr = ();
	for($i=0;$i<@title;$i++) {
	    my $seq = substr($data[$i],0,$window);
	    $arr[$i]="$title[$i] $seq";
	    substr($data[$i],0,$window)="";
	}
	my $b = " ".substr($bord,0,$window);
	substr($bord,0,$window)="";
	$res.=tex_tabular(\@arr, $b);
    }
    return($res);
}


sub features {
# This function takes sequence with gaps (0) and looks for signatures in it given the id number (1)
# Returns spaces in the places where features were residing
# Input: Seq, ID
# output Seq (with spaces), 
    my $seq = @_[0];
    my $id  = @_[1];
    $seq =~s/\-//g;
    $seq =~tr/[a-z]/[A-Z]/;
    my @arr=();
    foreach my $z(@{$SITES{$GENE{$id}}}) {
	my $l = index($seq,$LSGN{$z}.substr($RSGN{$z},0,1));
	my $r = index($seq,substr($LSGN{$z},-1,1).$RSGN{$z});
	$r++ if($r>=0);
	my $k = undef;
	$k = $r if($l>=0 && $r>=0 && $l+length($LSGN{$z})==$r);
	$k = $r if($r>=0 && $r<length($LSGN{$z}));
	$k = $l + length($LSGN{$z}) if($l>=0 && $l + length($LSGN{$z}) + length($RSGN{$z}) >= length($seq));
	if($k) {
	    substr($seq,$k,0) = " ";
	    push @arr, $z;
        }
    }
    print STDERR "$seq\n" if(@_[0] =~/cgcgtaataaat/);
    unshift @arr, ($arr[0]-1);
#back from $seq (without gaps) to @_[0] (with gaps)
    my @a = split //, $seq;
    my @b = split //, @_[0];
    for(my $i=0, my $j=0; $i<@a; $i++) {
	for(;$j<@b && $b[$j] eq "-";$j++) {;}
	if($a[$i] eq " ") {
	    $b[$j] = " $b[$j]";
	}
	else {
	    $j++;
	}
    }
    my $f = join("",@b);
    my @c = split /\s/, $f;
    my @d = ();
    for($i=0;$i<@c;$i++) {
	$s = ($SEGM{$arr[$i]}=~/^[EA].$/) ? "-" : ".";
	push @d, ($s x length($c[$i]));
    }
    $g = join("",@d);
    return($f, $g);
}
