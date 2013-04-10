#!/usr/bin/perl
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

use Perl::utils;
use Perl::binrel;

if(@ARGV==0) {
    print STDERR "This script transforms .gtf (gencode or ensembl) into .glc format; details can be found in README file\n";
    print STDERR "Usage: $0 -i gtf -full don't remove the dot part of the gene name\n"; 
    exit(0);
}

for(my $i=0;$i<@ARGV;$i++) {
    $gtf  = $ARGV[++$i] if($ARGV[$i] eq "-in");
    $full = 1 if($ARGV[$i] eq "-full");
}

die("gtf file cannot be opened, exiting\n") unless(-e $gtf);

print STDERR "[Input set to: $gtf]\n";
open FILE, $gtf;
seek(FILE,0,2);
$endpos = tell(FILE);
seek(FILE,0,0);

while($line=<FILE>) {
    $currpos = tell(FILE);
    progressbar(int($currpos/1000), int($endpos/1000)+1,"Reading exons\t");
    next if($line=~/^\#/);
    chomp($line);

    ($chr, $source, $tag, $beg, $end, $dot, $str, $dot, $rest) = split /\s*\t\s*/, $line;
    next unless($tag eq "exon");

    %feature = ();
    foreach $field(split /\s*;\s*/, $rest) {
        ($key, $value) = split /\s*\"\s*/, $field;
        $feature{$key}=$value if($key=~/\w/);
    }

    $gene_id = $feature{'gene_id'};
    $transcript_id = $feature{'transcript_id'};

    unless($full) {
        ($gene_id) = split /\./, $gene_id;
        ($transcript_id) = split /\./, $transcript_id;
    }

    $chr = "chr$chr" unless($chr=~/^chr/);
    assure($chromosome{$transcript_id}, $chr),
    $chromosome{$transcript_id} = $chr;

    $str = ($str eq "+") ? 1 : -1;
    assure($strand{$transcript_id}, $str);
    $strand{$transcript_id} = $str;

    ($beg, $end) = reverse ($beg, $end) if($str<0);
    push @{$exons{$transcript_id}}, [$beg,$end] if($tag eq "exon");

    $transcript_gene{$transcript_id}{$gene_id} = $gene_transcript{$gene_id}{$transcript_id} = 1;
}

progressbar(100,100,"Reading exons\t");

$n = 0+keys(%exons);
foreach $transcript_id(keys(%exons)) {
    progressbar(++$j, $n,"Processing\t");
    @arr = sort{$a->[0]<=>$b->[0]} @{$exons{$transcript_id}};
    $chr = $chromosome{$transcript_id};
    $str = $strand{$transcript_id};
    @arr = reverse @arr if($str<0);

    for($i=0;$i<@arr;$i++) {
	$transcript_acc{$transcript_id}{"$chr\t$arr[$i]->[0]\t$str"}=$acc_transcript{"$chr\t$arr[$i]->[0]\t$str"}{$transcript_id}=1 if($i>0);
        $transcript_don{$transcript_id}{"$chr\t$arr[$i]->[1]\t$str"}=$don_transcript{"$chr\t$arr[$i]->[1]\t$str"}{$transcript_id}=1 if($i<@arr-1);
    }
}

print STDERR "[Transcript-donor and donor-transcript relations";
%rel1 = product(\%transcript_don,\%don_transcript);
print STDERR "]\n";

print STDERR "[Transcript-acceptor and acceptor-transcript relations";
%rel2 = product(\%transcript_acc,\%acc_transcript);
adjoint(\%rel1,\%rel2);
print STDERR "]\n";

print STDERR "[Transcript-gene and gene-transcript relations";
%rel3 = product(\%transcript_gene,\%gene_transcript);
adjoint(\%rel2,\%rel3);
print STDERR "]\n";

print STDERR "[Transitive closure of transcript-transcript relation";
%res = cliques(\%rel3);
print STDERR "]\n";

foreach $transcript_id(keys(%res)) {
    print "$transcript_id\t$res{$transcript_id}\n";
}

