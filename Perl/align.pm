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

###########################################################################################################
# This package contains various subroutines for intelligent alignments of sequenes with boxes
###########################################################################################################

    $pdflatex="pdflatex";		#name of the pdflatex program
    $muscle="./muscle";			#name of the muscle program;
    $tmpfilename = "tmp";		#temporary file + PID will be used by muscle and pdflatex
    $alnparams = "-maxiters 1 -diags";	#parameters passed to the muscle alignmer
    $INFTY = 65535;			#plus infinity
    $MAXALNLEN = 65;			#sequences longer than this limit are cut in halves, only ALNWINDOW from the left and right are left
    $ALNWINDOW = 20;			#then the left and the right parts are alignmed separately and the resulting alignments are merged
    $MAXLENMUSCLE = 800;		#max sequence length tolerated by muscle
    $ASTERISK_SCALE = 1;		#grey gradations of asterisks
    $ASTERISK_POWER = 1;		#=alpha, where black_intensity=relative_entropy^alpha
    $TEXWINDOW = 200;			#width of the alignment to be cut per 1 line

    return(1);
###########################################################################################################

sub hash2fasta {
# this procedure saves a hash into a FASTA file
# input:  (hash, filename)
# output: number of non-empty sequences
    my $filename = pop(@_);
    my %hash = @_;
    my $n=undef;
    my $flag=undef;
    open XFILE, ">$filename" || die("can't open $filename");
    foreach $key(keys(%hash)) {
	if($hash{$key}=~/(\w)/) {
	    if(length($hash{$key})>$MAXLENMUSCLE) {
		$hash{$key} = substr($hash{$key},0,$MAXLENMUSCLE);
		$flag=1;
	    }
	    print XFILE ">$key\n$hash{$key}\n";
	    $n++;
	}
    }
    close XFILE;
    $nwarnings++ if($flag);
    return($n);
}

sub fasta2hash {
# this procedure reads FASTA into a hash
# input:  filename
# output: hash
    my $data = undef;
    my %hash = ();
    open XFILE, @_[0] || die("can't open @_[0]");

    while($line=<XFILE>) {
	$data.=$line;
    }
    close XFILE;
    my @array = split /\>/, $data;
    shift(@array);
    foreach $key(@array) {
	($id, $rest) = split /\n/, $key, 2;
	$id   =~ s/\s//g;
	$rest =~ s/\n//g;
	$hash{$id} = $rest;
    }
    return(%hash);
}

sub matchcase {
# this procedure reads line in @_[0] (gaps allowed) and transforms cases according to @_[1] (no gaps)
# input 1: string, i.e. AGGGTACGATCTAGTCGATCGAG
# input 2: string, i.e. AAAAAAAAA......AAAAAAAA
# output: 		AGGGTACGAtctagtCGATCGAG
    my @a = split //, @_[0];
    my @b = split //, @_[1];
    my $i=0;
    my $j=0;
    while($i<@a) {
	while($a[$i] eq "-" && $i<@a) {$i++;}
	if($i<@a) {
	    $a[$i] =~ tr/[a-z]/[A-Z]/ if($b[$j]=~/[A-Z]/);
	    $a[$i] =~ tr/[A-Z]/[a-z]/ if($b[$j]=~/[a-z.]/);
	}
	$i++;
	$j++;
    }
    return(join("",@a));
}

sub align {
# this procedure passes a hash (NO BLOCKS!) to the alignment program and returns the aligned hash
# input:  hash sequences not aligned
# output: hash sequences aligned
    my %hash = @_;
    my %res = ();
    if(hash2fasta(%hash,"$tmpfilename$$.fa")) {
    	system("$muscle -in $tmpfilename$$.fa -out $tmpfilename$$.aln -quiet $alnparams");
    	%res = fasta2hash("$tmpfilename$$.aln");
    	foreach my $key(keys(%res)) {
	    $res{$key} = matchcase($res{$key},$hash{$key})
    	}
    	system("rm -f $tmpfilename$$.aln $tmpfilename$$.fa");
        return(%res);
    }
    else {
	return(@_);
    }
}


sub reduce {
#this subroutine checkes if the aligned block is not too wide, and if it is then it removes a piece in the middle
# in theory a better way would be to detect bad alignments i.e. ones which picked up a wrong seed
    my %hash = @_;
    my $n = 0;
    my $MAXALNLEN1 = $MAXALNLEN*2;
    my $ALNWINDOW1 = $ALNWINDOW*2;
    my %res=();
    foreach my $key(keys(%hash)) {
        my $l = length($hash{$key});
        $n = $l if($l>$n);
    }
    if($n<$MAXALNLEN1) {
        return(%hash);
    }
    foreach my $key(keys(%hash)) {
	$res{$key} = substr($hash{$key},0,$ALNWINDOW1)."::::::".substr($hash{$key},-$ALNWINDOW1,$ALNWINDOW1);
    }
    return(%res);
}

sub malign {
# same as align if the shortest of the sequences in the hash is less that $MAXALNLEN
# if not, $ALNWINDOW from the left of each sequnce are cut and aligned, same from the right, the results are docked with ...
# input:  hash sequences not aligned
# output: hash sequences aligned
    my %hash = @_;
    my $n = $INFTY;
    foreach my $key(keys(%hash)) {
	my $l = length($hash{$key});
	$n = $l if($l<$n);
    }
    if($n<$MAXALNLEN) {
	return(align(%hash));
    }
    my %left=();
    my %right=();
    my %insert=();
    my $ni=0;
    foreach my $key(keys(%hash)) {
	$left{$key}  = substr($hash{$key},0,$ALNWINDOW);
	$right{$key} = substr($hash{$key},-$ALNWINDOW,$ALNWINDOW);
	$insert{$key}= length($hash{$key}) - 2*$ALNWINDOW;
	$ni = length($insert{$key}) if(length($insert{$key})>$ni);
    }
    %left  = align(%left);
    %right = align(%right);
    my %res=();
    foreach my $key(keys(%hash)) {
	$res{$key} = "$left{$key}...".("." x ($ni-length($insert{$key})))."$insert{$key}...$right{$key}";
    }
    return(%res);
}

sub zalign {
# yet to be done better: block alignment with RE-alignment of sequences where no seedds have been found
    my %hash1 = %{@_[0]};
    my %hash2 = %{@_[1]};
    my %hash = ();
    my %res = ();
    %hash = qalign(%hash1);
    if(hash2fasta(%hash,"$tmpfilename$$.profile")) {
        return(%hash) unless(hash2fasta(%hash2,"$tmpfilename$$.fa"));
        system("muscle -profile -in1 $tmpfilename$$.profile -in2 $tmpfilename$$.fa -out $tmpfilename$$.aln -quiet $params");
        %res = fasta2hash("$tmpfilename$$.aln");
        foreach my $key(keys(%res)) {
            $res{$key} = matchcase($res{$key},$hash1{$key}.$hash2{$key});
        }
        system("rm -f $tmpfilename$$.profile $tmpfilename$$.aln $tmpfilename$$.fa");
        return(%res);
    }
}

######################################################################################################

sub qsplit {
# this procedure takes sequences stored in a hash, where sequecnes have breaks (spaces) and splits it
# into an array of sequne hashes without brakes. terminates if block sizes are different
# input:  sequence hash
# output: 2D mixed array/hash; 1st key - index, second key - seq_id
    my %hash = @_;
    my @data = ();
    my $n=undef;
    foreach my $key(keys(%hash)) {
        my @array = split /\s+/, $hash{$key};
        for(my $i=0;$i<@array;$i++) {
            $data[$i]{$key} = $array[$i];
        }
        $n = 0 + @array unless($n);
        #die("different number of blocks for qalign") unless($n==0+@array);
    }
    return(@data);
}

sub qjoin {
# this procedure joins blocks in 2D mixed array/hash; 1st key - index, second key - seq_id
# into a single hash by ids in the first block (with index = 0)
# inverse to qsplit
# input:  2D mixed array/hash; 1st key - index, second key - seq_id
# output: sequence hash
    my @data = @_;
    my %res = ();
    my @array = keys(%{$data[0]});
    foreach my $key(@array) {
	for(my $i=0;$i<@data;$i++) {
	    $res{$key}.=$data[$i]{$key};
	}
    }
    return(%res);
}

sub qdock {
# this procedure docks two alignments given by blocks references by hash1 and hash2
# by removing dashes at the right end of hash1 and at the left end of hash2
# input:  (hash REFERENCE to the first alingment, hash REFERENCE to the second alignment)
# output: number of dashes removed
    my %hash1 = %{@_[0]};
    my %hash2 = %{@_[1]};
    my $n = $INFTY;
    foreach my $key(keys(%hash1)) {
	my $x = ($hash1{$key} =~ /(\-+)$/) ? length($1) : 0;
	my $y = ($hash2{$key} =~ /^(\-+)/) ? length($1) : 0;
	$n = $x + $y if($n > $x + $y);
    }
    foreach my $key(keys(%hash1)) {
	my $x = ($hash1{$key} =~ /(\-+)$/) ? length($1) : 0;
	my $y = ($hash2{$key} =~ /^(\-+)/) ? length($1) : 0;
	my $k = ($x<$n ? $x : $n);
	my $l = $n - $k;
	substr(${@_[0]}{$key},-$k) = "" if($k>0);
	substr(${@_[1]}{$key},0,$l) = "" if($l>0);
    }
    return($n);
} 

sub qalign {
# does block alignment, docking, and union
# input:  hash
# output: hash
    my %hash = @_;
    my %res = ();
    my @data = qsplit(%hash);
    for(my $i=0;$i<@data;$i++) {
	%{$data[$i]} = malign(%{$data[$i]});
	%{$data[$i]} = reduce(%{$data[$i]});
    }
    for(my $i=0;$i<@data-1;$i++) {
        qdock(\%{$data[$i]},\%{$data[$i+1]});
    }
    %res = qjoin(@data);
    return(%res);
}

#############################################################################################################################

sub compile_tex {
# runs TEX compiler on the body supplied in the first argument
# input:  (filename, tex_text, number of iterations)
    my $filename = @_[0];
    my $tex_text = @_[1];
    my $n_passes = @_[2];
    my $line;
    $n_passes = 1 unless($n_passes);

    #$tex_text =~ s/\_/\\\_/g;
    open  XFILE, ">$filename.tex";
    print XFILE $tex_text;
    close XFILE;

    my @aux = split /\//, $filename;
    pop(@aux);
    my $param = @aux > 0 ? "-output-directory=".join("/",@aux) : undef;
    for(my $i=1;$i<=$n_passes;$i++) {
        print STDERR "Latex pass $i ";
	open XFILE, "$pdflatex $param $filename |";
	while($line=<XFILE>) {
	    if($line=~/^\!/ && $line!~/^\! Package pgfplots Warning/) {
		print STDERR " encountered problem (see $filename.aux), exiting\n$line";
		exit(1);
	    }
	    print STDERR "." if($line=~/\[(\d+)/);
	}
	close XFILE;
	print STDERR "\n";
    }
    system("rm -f $filename.tex $filename.aux $filename.log $filename.out") unless(@_[3]);
}


sub tex_document {
# Produces TEX document
# input:  (body, preamble)
# output: tex_document
#\\textwidth=21cm
#\\textheight=27cm
#\\hoffset=-4cm
#\\voffset=-3cm
(my $body, my $head) = @_;
$head.="\\usepackage{pgfplots}\n" if($PGFPLOTS);
return <<END
\\NeedsTeXFormat{LaTeX2e}
\\documentclass[english,a4paper,onecolumn,10pt,rotating,landscape]{article}
\\usepackage{longtable}
\\usepackage[usenames,dvipsnames]{xcolor}
\\usepackage{soul}
\\usepackage{lscape}
\\usepackage{rotating}
\\usepackage{hyperref}
\\usepackage{tikz}
\\renewcommand{\\baselinestretch}{0.9}
\\textwidth=27cm
\\textheight=20cm
\\hoffset=-8cm
\\voffset=-4cm
$head
\\begin{document}
$body
\\end\{document\}
END
}

sub tex_table {
# Produces TEX table
# input: array of space-separated columns 
# output: table_text 
my $caption = @_[1] ? "\\caption\{@_[1]\}" : undef;
return <<END
\\begin\{table\}
\\begin\{center\}
\\tt\\footnotesize
@_[0]
\\end\{center\}
\\footnotesize
$caption
\\end\{table\}\n
END
}

sub highlight {
    my $z = @_[0];
    $z =~ s/([a-z\&\-\:\.])([A-Z])/$1\\hl\{$2/g;
    $z =~ s/^([A-Z])/\\hl\{$1/;
    $z =~ s/([A-Z])([a-z\&\-\:\.])/$1\}$2/g;
    $z =~ s/([A-Z])$/$1\}/;
    return($z);
}


sub tex_tabular {
# Produces TEX table
# input: array of space-separated columns 
# output: table_text 
    my @table = @{@_[0]};
    my @res = ();
    my $n=undef;
 
    foreach my $line(@table) {
	my @a = split /\s/, $line;
	$n = 0 + @a unless($n);
	die("Dimension error") unless($n == @a);
	my $start = shift(@a);
        my $rest = highlight(join("\&",@a));
	$rest =~ s/\&/\\,\&\\,/g;
	$rest = replace_asterisks($rest) unless($rest=~/[A-Z]/i);
	push @res, "$start\&$rest";
    }

    my $header = "l". ("c\@{\\,}|\@{\\,}" x ($n-2)). "l";
    my $header = "l". ("c\@{}|\@{}" x ($n-2)). "l";

    my @a = split /\s/, @_[1];
    print STDERR "!" unless(@a == $n);
    my $borders = undef;
    for($i=1;$i<=@a;$i++) {
	$borders.= "\\cline{$i-$i}" if($a[$i-1]=~/-/);
    }
    my $top = $borders. ("\&" x ($n-1)). "\\\\\n" if($borders);
    my $bot = $borders if($borders);
    #$top=$bot=undef unless(@a == $n);


    my $output = join("\\\\\n",@res)."\\\\";
    return("\\begin\{tabular\}\{$header\}\n$top\n$output\n$bot\\end\{tabular\}\n");
}

sub landscape {
    return("\\begin\{landscape\}\n@_[0]\n\\end\{landscape\}");
}

sub enthropy {
    my %hash=@_;
    my $s=0;
    my $e=0;
    my $m=0;
    my $c=undef;
    foreach $val(keys(%hash)) {
	$e+= ($hash{$val}>0) ? $hash{$val} * log($hash{$val}) : 0;
	$s+= $hash{$val};
	if($hash{$val}>$m) {
	    $m=$hash{$val};
	    $c=$val;
	}
    }
    return(0) unless($s>0);
    $e=int($ASTERISK_SCALE*($e/$s/log($s)));
    return($c =~ /[a-z]/ ? $e : 0);
}

sub asterisks {
    my @table = @_;
    my $n=0;
    my %frequency;
    my $res="";

    for(my $i=0;$i<@table;$i++) {
        $table[$i] =~ tr/[A-Z]/[a-z]/;
        $f[$i] = [ split //, $table[$i] ];
        $n = length($table[$i]) if($n<length($table[$i]));
    }

    for(my $j=0;$j<$n;$j++) {
        %frequency=();
        for(my $i=0;$i<@table;$i++) {
            $frequency{$f[$i][$j]}++;
        }
        substr($res,$j,1) = enthropy(%frequency);
    }
    return($res);
}

sub replace_asterisks {
    my @arr = split //, @_[0];
    for(my $i=0;$i<@arr;$i++) {
	$arr[$i]="{\\color{black!".int((($arr[$i]/$ASTERISK_SCALE)**$ASTERISK_POWER)*100)."}*}" if($arr[$i] =~/\d/);
    }
    return(join("",@arr));
}


sub replace_extention {
    my $s = @_[0];
    $s =~ s/\.\w+$//;
    return($s.@_[1]);
}
