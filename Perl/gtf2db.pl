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
    print STDERR "This script transforms .gtf (gencode or ensembl) into .cps format; details can be found in README file\n";
    print STDERR "Usage: $0 -i gtf input -db database file name (without .cps) -test number of lines to read in test regime";
    print STDERR " -level max transcript level -full don't remove the dot part of the gene name -only require that transcript type be the same as gene biotype\n";
    exit(0);
}

$maxl = 65535;

for(my $i=0;$i<@ARGV;$i++) {
    $gtf  = $ARGV[++$i] if($ARGV[$i] eq "-i");
    $ref  = $ARGV[++$i] if($ARGV[$i] eq "-db");
    $test = $ARGV[++$i] if($ARGV[$i] eq "-test");
    $maxl = $ARGV[++$i] if($ARGV[$i] eq "-level");
    $glc  = $ARGV[++$i] if($ARGV[$i] eq "-g");
    $full = 1 if($ARGV[$i] eq "-full");
    $only = 1 if($ARGV[$i] eq "-only");
}

die("gtf file cannot be opened, exiting\n") unless(-e $gtf);
die("database name not specified, exiting\n") unless($ref);

open SING, ">$ref.cps";
open INTR, ">$ref.int";
open EXON, ">$ref.exn";
open SEGM, ">$ref.sgm";

if($glc) {
    print STDERR "[Definition of genomic loci: $glc";
    %relation = split /[\n\t]/, `cat $glc`;
    print STDERR "]\n";
}

print STDERR "[Input set to: $gtf]\n";
open FILE, $gtf;
seek(FILE,0,2);
$endpos = tell(FILE);
seek(FILE,0,0);
while($line=<FILE>) {
    $currpos = tell(FILE);
    progressbar(int($currpos/1000), int($endpos/1000)+1,"Reading gtf\t");
    next if($line=~/^\#/);
    chomp($line);
    ($chr, $source, $tag, $beg, $end, $dot, $str, $dot, $rest) = split /\s*\t\s*/, $line;
    next unless($tag eq "exon" || $tag eq "CDS");
    $chr = "chr$chr" unless($chr=~/^chr/);
    $str = ($str eq "+") ? 1 : -1;
    %feature = ();
    foreach $field(split /\s*;\s*/, $rest) {
	($key, $value) = split /\s*\"\s*/, $field;
	$feature{$key}=$value if($key=~/\w/);
    }
    $feature{'gene_type'} = $feature{'gene_biotype'} unless($feature{'gene_type'});
    $feature{'gene_type'} = $source unless($feature{'gene_type'});
    next unless($feature{'gene_id'} && $feature{'gene_type'} && $feature{'transcript_id'});
    $gene_id = $feature{'gene_id'};
    $transcript_id = $feature{'transcript_id'};
    unless($full) {
        ($gene_id) = split /\./, $gene_id;
        ($transcript_id) = split /\./, $transcript_id;
    }

    $gene_id = $relation{$transcript_id} if($relation{$transcript_id});

    next unless($gene_id);
    if($tag eq "exon") {
	assure($chromosome{$gene_id}, $chr);
	assure($strand{$gene_id}, $str);
	$gene_type{$gene_id}  = $feature{'gene_type'};
        $genename{$gene_id}   = $feature{'gene_name'} ? $feature{'gene_name'} : "NA";
	$chromosome{$gene_id} = $chr;
	$strand{$gene_id}     = $str;

        push @{$transcripts{$gene_id}}, $transcript_id unless($exon_beg{$transcript_id});
	$transcript_type{$transcript_id} = $feature{'transcript_type'};
	$level{$transcript_id}   = $feature{'level'};
        push @{$exon_beg{$transcript_id}}, $beg;
        push @{$exon_end{$transcript_id}}, $end;
    }
    if($tag eq "CDS") {
        push @{$cds_beg{$transcript_id}}, $beg;
        push @{$cds_end{$transcript_id}}, $end;
    }
    $n++;
    last if($n>$test && $test>0);
}
progressbar(100,100,"Reading gtf\t");

foreach $gene_id(keys(%transcripts)) {
    progressbar(++$num, 0+keys(%transcripts),"Processing\t");
    $chr = $chromosome{$gene_id};
    $str = $strand{$gene_id};

    # scan through all transcripts of the given gene (duplicate transcripts excluded)
    # position means exon's 1st nt or intron's 1st nt (not exon's last nt!)
    # exons, introns, splice sites, transcript start and ends are protocoled
    # %site keeps the number of transcripts in which the given position was used as a SS
    # %exon stores the other end of the SHORTEST exon ending at the given site

    %site = %type = %min = %max = %exon = %intr = ();
    foreach $transcript_id(@{$transcripts{$gene_id}}) {
        next if($transcript_type{$transcript_id} ne $gene_type{$gene_id} && $only);
	next unless($level{$transcript_id}<=$maxl);

        my @beg = sort {$a<=>$b} @{$exon_beg{$transcript_id}};
        my @end = sort {$a<=>$b} @{$exon_end{$transcript_id}};
        if($str<0) {
            @tmp = @beg; @beg = reverse @end; @end = reverse @tmp;
        }
        for(my $i=0;$i<@beg;$i++) {
            $end[$i]+=$str;
            $site{$beg[$i]}++;
            $site{$end[$i]}++;
            $type{$beg[$i]} = ($i>0 ? "A" : "S")        unless($type{$beg[$i]} eq "A");
            $type{$end[$i]} = ($i<@beg-1 ? "D" : "E")   unless($type{$end[$i]} eq "D");
	    $exon{$beg[$i]."\t".$end[$i]}{$transcript_id}++;
            $intr{$end[$i-1]."\t".$beg[$i]}{$transcript_id}++ if($i>0);
        }
    }
    @array = sort{$a<=>$b} keys(%site);
    @array = reverse @array if($str<0);

    %site_tot = %segm_tot = %segm_exn = %segm_cds = %segm_tp = %segm_lft = %segm_rgh = %coverage = ();

    foreach $transcript_id(@{$transcripts{$gene_id}}) {
        next if($transcript_type{$transcript_id} ne $gene_type{$gene_id} && $only);
        next unless($level{$transcript_id}<=$maxl);

	# This part concerns segments' inclusion and exclusion into pre-mRNAs
	# and with the property of a segment to be coding or non-coding
	# %site_tot keeps the total number of transcripts COVERING the site
	# %segm_tot keeps the total number of transcripts COVERING the segment [site, site+1]
	# %segm_exn keeps the number of transcripts in which [site, site+1] is exonic
	# %segm_cds keeps the number of CDSs that OVERLAP with the segment [site, site+1]
	# %segm_lft is the number of transcripts in which the segment [site, site+1] is in the left  UTR (5' for +, 3' for -)
	# %segm_rgh is the number of transcripts in which the segment [site, site+1] is in the right UTR (3' for +, 5' for -)

        my @beg  = sort {$a<=>$b} @{$exon_beg{$transcript_id}};
        my @end  = sort {$a<=>$b} @{$exon_end{$transcript_id}};
        my @cbeg = sort {$a<=>$b} @{$cds_beg{$transcript_id}};
        my @cend = sort {$a<=>$b} @{$cds_end{$transcript_id}};

	if($str<0) {
            @tmp = @beg;  @beg  = reverse @end;  @end  = reverse @tmp;
            @tmp = @cbeg; @cbeg = reverse @cend; @cend = reverse @tmp;
        }

        for(my $i=0;$i<@beg;$i++) {
            $end[$i]+= $str;
        }
        for(my $i=0;$i<@cbeg;$i++) {
	    $cend[$i]+=$str;
	}
        ($min,  $max)  = sort{$a<=>$b} ($beg[0],  $end[@end-1]);
	($cmin, $cmax) = sort{$a<=>$b} ($cbeg[0], $cend[@cend-1]);
        for(my $j=0;$j<@array;$j++) {
            $site_tot{$array[$j]}++ if($min<=$array[$j] && $array[$j]<=$max);
        }
        for(my $j=0;$j<@array-1;$j++) {
            ($x, $y) = sort {$a<=>$b} ($array[$j], $array[$j+1]);
            $segm_tot{$array[$j]}++ if($min<=$x && $y<=$max);
	    $segm_lft{$array[$j]}++ if($min<=$x && $y<$cmin && $cmin>0);
	    $segm_rgh{$array[$j]}++ if($cmax<$x && $y<=$max && $cmax>0);
            for(my $i=0;$i<@beg;$i++) {
                ($p, $q) = sort {$a<=>$b} ($beg[$i], $end[$i]);
                $segm_exn{$array[$j]}++ if($p<=$x && $y<=$q);
            }
	    for(my $i=0;$i<@cbeg;$i++) {
		($p, $q) = sort {$a<=>$b} ($cbeg[$i], $cend[$i]);
		$segm_cds{$array[$j]}++ unless($y<=$p || $q<=$x);
	    }
	    $coverage{$array[$j]}{$array[$j+1]}{$transcript_id}++ if($min<=$x && $y<=$max);
        }
    }

    for(my $j=0;$j<@array-1;$j++) {
	print STDERR "!" if($segm_exn{$array[$j]}>$segm_tot{$array[$j]});
	$segm_tp{$array[$j]} = ($segm_tot{$array[$j]} > 0) ? ( ($segm_exn{$array[$j]} > 0) ? ( ($segm_exn{$array[$j]} < $segm_tot{$array[$j]}) ? "A" : "E"  ) : "I" ) : "X";
	$z = ($segm_cds{$array[$j]} > 0) ? "C" : "N";
	$z = ($str>0 ? "5" : "3") if($segm_tot{$array[$j]}==$segm_lft{$array[$j]});
	$z = ($str>0 ? "3" : "5") if($segm_tot{$array[$j]}==$segm_rgh{$array[$j]});
	$segm_tp{$array[$j]}.= $z; #($segm_cds{$array[$j]} > 0) ? "C" : "N";
    }

    ++$gn;
    %sidx = ();
    for(my $j=0;$j<@array;$j++) {
	++$sn;
	print SING "$chr\t$array[$j]\t$str\t$gn\t$sn\t$type{$array[$j]}\t",($j>0?abs($array[$j]-$array[$j-1]):0),"\t",($j+1<@array?abs($array[$j+1]-$array[$j]):0),"\t";
	print SING "$gene_id\t$gene_type{$gene_id}\t$genename{$gene_id}\t",$site{$array[$j]}+0,"\t",$site_tot{$array[$j]}+0,"\t";
	print SING ($j>0 ? $segm_tp{$array[$j-1]}:"NA"),"\t",($j+1<@array ? $segm_tp{$array[$j]}:"NA"),"\t";
        print SING ($j>0 ? segm_repr($chr,$array[$j-1],$array[$j],$str):"NA"), "\t", ($j+1<@array ? segm_repr($chr,$array[$j],$array[$j+1],$str) : "NA"),"\n";
	$sidx{$array[$j]} = $sn;
	if($j>0) {
	    $segm = segm_repr($chr,$array[$j-1],$array[$j],$str);
	    next if($been{$segm});
	    $been{$segm}=1;
	    print SEGM "$segm\t$gene_id\n";
	}
    }
    foreach my $i(keys(%intr)) {
	(my $d, my $a) = split /\t/, $i;
	print INTR "$sidx{$d}\t$sidx{$a}\t",join(",",keys(%{$intr{$i}})),"\t",join(",",keys(%{$coverage{$d}{$a}})),"\n";
    }
    foreach my $e(keys(%exon)) {
        (my $a, my $d) = split /\t/, $e;
        print EXON "$sidx{$a}\t$sidx{$d}\t",join(",",keys(%{$exon{$e}})),"\t",join(",",keys(%{$coverage{$a}{$d}})),"\n";
    }
}
close SEGM;
close SING;
close INTR;
close EXON;

sub segm_repr {
    (my $chr, my $x, my $y, my $strand) = @_;
    $y-= $strand;
    ($x, $y) = sort {$a<=>$b} ($x, $y);
    return("$chr\_$x\_$y\_$strand");
}

