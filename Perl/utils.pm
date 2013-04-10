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


use Perl::setup;

##########################################################################################
# This package contains subroutines
# 1) load_global_configuration_file
# 2) read_configuration
# 3) read_annotation
# 4) read_introns
# 5) read_signatures
# 6) read_sequences
# 7) read_expression
###########################################################################################

    %ABBR = %UCSC = %DATATYPE = %REFDB = %CLADE = %NICKNAME = %FULLNAME = ();

    load_global_configuration_file();
    return(1);

##########################################################################################

sub progressbar {
    my ($current, $last, $message) = @_;
    my $width = 2**int(log($WCHAR)/log(2));
    my $i;
    if(int(($width*($current-1))/$last) < int(($width*$current)/$last)) {
        my $k = int(($width*$current)/$last);
        print STDERR "\r$message\[";
        for($i=0;$i<$k;$i++) {print STDERR "=";}
        print STDERR ">" if($k<$width);
        for($i++;$i<$width;$i++) { print STDERR " ";}
        print STDERR "\] ",($current/$last < 0.1 ? " " : ""), int(100*$current/$last),"%";
    }
    print STDERR "\n" if($current==$last);
}

sub assure {
    (my $a, my $b) = @_;
    die("Conflict: $a vs $b\n") if($a && $a ne $b);
}

sub revcomp {
    my $res = join("", reverse split //,@_[0]);
    $res =~ tr/[acgt]/[tgca]/;
    $res =~ tr/[ACGT]/[TGCA]/;
    return($res);
}

sub randseq {
    my $res=undef;
    my @arr = split //, "acgt";
    for(my $i=0;$i<@_[0];$i++) {
        $res.=$arr[int(rand()*4)];
    }
    return($res)
}

sub make {
    my %param=@_;
    my @input =();
    my @commandline = ();
    my @output = ();
    my %mkdir = ();
    push @commandline, $param{'before'};
    push @input, $param{'depends'}; 
    foreach $key(keys(%{$param{'input'}})) {
        push @commandline, ($key, $param{'input'}{$key}) if($param{'input'}{$key} =~/\w/);
        push @input, $param{'input'}{$key};
    }
    push  @commandline, $param{'after'};
    foreach $key(keys(%{$param{'output'}})) {
        push @commandline, ($key, $param{'output'}{$key});
        push @output, $param{'output'}{$key};
	my @arr = split /\//, $param{'output'}{$key};
	pop(@arr);
	$mkdir{join("/",@arr,'')}++;
    }
    #push @input, $param{'script'} unless($param{'script_not_required'});
    $param{'script'}= "perl $param{'script'}" if($param{'script'}=~/\.pl$/);
    print join(" ",@output)," : ",join(" ",@input), "\n";
    print "\tmkdir -p ",join(" ", keys(%mkdir)), "\n";
    print "\t$param{'script'} ",join(" ",@commandline)," \n";
    print "$param{'group'} :: ", join(" ",@output), "\n" if($param{'group'});
}

##########################################################################################

sub load_global_configuration_file {
    my $cfg_file = "config.dat";
    open FILE, $cfg_file || die("Config file ($cfg_file) does not exist, exiting\n");
    while(my $line=<FILE>) {
	next if($line=~/^#/);
	chomp $line;
	(my $id) = split /\s+/, $line;
	(my $id, $ABBR{$id}, $UCSC{$id}, $DATATYPE{$id}, $REFDB{$id}, $CLADE{$id}, $NICKNAME{$id}, $FULLNAME{$id}) = split /\s+/, $line;
	push @{$CLADELIST{$CLADE{$id}}}, $id;
	$BASESPECIES{$CLADE{$id}} = $id unless($BASESPECIES{$CLADE{$id}});
    }
    close FILE;
}

###########################################################################################################
sub read_configuration {
# Read configuration file specified in $config_file_name
# Populate variables $annotation_file_name and $introns_file_name and $path
    $config_file_name = @_[0] if(@_[0]);
    return unless($config_file_name);
    print STDERR "[Configuration:\t$config_file_name";
    unless($filewasread{$config_file_name}) {
	$path = $annotation_file_name = $introns_file_name = $signatures_file_name = $extention = undef;
	%species = %index = ();
    	$nsp = 0;
    	open FILE, $config_file_name || die(" - can't open '$config_file_name'\n");
    	while($line=<FILE>) {
	    chomp $line;
            ($term, $val, $weight) = split /\s+/, $line;
	    $halfsize		  = $val if($term eq "halfsize");
	    $gapsize		  = $val if($term eq "gap");
	    if($term eq "path") {
		$path = $val;
		$path = "$ENV{HOME}$1" if($path=~/^\~(.*)$/);
	    }
            $annotation_file_name = $path.$val if($term eq "annotation");
            $introns_file_name    = $path.$val if($term eq "introns");
	    $signatures_file_name = $path.$val if($term eq "signatures");
            $extention		  = $val if($term eq "extention");
            $index{$val}	  = $nsp if($term eq "species");
	    $weigh{$val}	  = $weight if($term eq "species");
            $species[$nsp++]	  = $val if($term eq "species");
	    $sumweigh             += $weight if($term eq "species");
	}
	close FILE;
	$filewasread{$config_file_name}=1;
    }
    print STDERR "]\n";
}

sub read_annotation {
    return unless($annotation_file_name);
    print STDERR "[Annotation:\t$annotation_file_name";
    unless($filewasread{$annotation_file_name}) {
    	open FILE, $annotation_file_name || die("- can't open '$annotation_file_name'\n");
    	while($line=<FILE>) {
            chomp $line;
	    ($chr, $pos, $str, $gene, $site, $type, $ups, $dws, $ensg, $biotype, $name, $use, $tot, $tup, $tdw, $sup, $sdw) = split /\t/, $line;
	    $DATA{$site}     = $line;
            $GENE{$site}     = $gene;
	    $SEGM{$site}     = $tdw;
	    $CHR{$site}	     = $chr;
	    $POSITION{$site} = $pos;
	    $STRAND{$site}   = $str;
	    $TYPE{$site}     = $type;
	    $DESC{$site}     = "gene=$ensg name=$name segment=$sdw type=$tdw";
	    $DESC{$site} =~ s/\_/\\\_/g; 
	    push @{$SITES{$gene}}, $site;
	    $NAME{$gene}     = $name;
        }
        close FILE;
	print STDERR " ",0+keys(%DATA);
	$filewasread{$annotation_file_name} = 1;
    }
    print STDERR "]\n";
}

sub getbox {
    my $site = @_[0];
    my $x = $POSITION{$site} + $STRAND{$site}*@_[1];
    my $y = $POSITION{$site} + $STRAND{$site}*@_[2];
    ($x, $y) = sort {$a<=>$b} ($x, $y);
    return("$CHR{$site}\_$x\_$y\_$STRAND{$site}");
}

sub read_introns {
    return unless($introns_file_name);
    print STDERR "[Introns:\t$introns_file_name";
    unless($filewasread{$introns_file_name}) {
    	open FILE, $introns_file_name || die(" - can't open '$introns_file_name'\n");
    	while($line=<FILE>) {
            chomp $line;
	    ($left, $right) = split /\t/, $line;
	    my $gene = $GENE{$left};
   	    next unless($gene);
	    push @{$INTRONS{$gene}}, [$left, $right];
        }
        close FILE;
	print STDERR " ",0+keys(%INTRONS);
	$filewasread{$introns_file_name}=1;
    }	
    print STDERR "]\n";
}

sub read_signatures {
    return unless($signatures_file_name);
    print STDERR "[Signatures:\t$signatures_file_name";
    unless($filewasread{$signatures_file_name}) {
    	open FILE, $signatures_file_name || die(" can't open '$signatures_file_name'\n");
    	while($line=<FILE>) {
            chomp $line;
	    ($site, $chr, $pos, $len, $str, $full, $seq) = split /\t/, $line;
            $seq =~ tr/[a-z]/[A-Z]/;
            $LSGN{$site} = substr($seq, 0,10);
            $RSGN{$site} = substr($seq,10,10);
        }
        close FILE;
	print STDERR " ",0+keys(%LSGN);
	$filewasread{$signatures_file_name}=1;
    }
    print STDERR "]\n";
}

sub read_sequences {
    my %selected = @_;
    if($maf_file_name) {
  	print STDERR "[MAF input $maf_file_name";
  	open FILE, $maf_file_name || die(" can't open '$maf_file_name'\n");
  	while($line=<FILE>) {
    	    if($line=~/^a id=(\d+)/) {
	    	$id = $1;
	    	next unless($selected{$id});
	    	while($line=<FILE>) {
	    	    ($s, $chr, $pos, $len, $str, $full, $seq) = split /\s+/, $line;
	    	    last unless($s eq "s");
	    	    ($name) = split /\./,$chr;
	    	    $org = $index{$name};
	    	    next if($org eq undef);
            	    $seq =~ tr/[A-Z]/[a-z]/;
	    	    $attr[$org][$id] = "$chr\t$pos\t$len\t$str\t$full";
		    $seq = revcomp($seq) if($rc);
	    	    $data[$org][$id] = $seq;
	    	}
    	    }
  	}
  	close FILE;
	print STDERR "(rc)" if($rc);
	print STDERR "]\n";
    }
    else {
  	for($i=0;$i<$nsp;$i++) {
    	    $filename = "$path$species[$i]$extention";
    	    print STDERR "[Reading $filename";
    	    open FILE, "$filename" || die(" can't open '$filename'\n");
	    if(-e "$filename.ind") {
		%ind = split /[\t\n]/, `cat $filename.ind`;
		print STDERR "(ind)";
		foreach $id(keys(%selected)) {
		    seek FILE, $ind{$id}, 0;
		    $line=<FILE>;
		    chomp($line);
		    ($id, $chr, $pos, $len, $str, $full, $seq) = split /\t/, $line;
		    $seq =~ tr/[A-Z]/[a-z]/;
		    $attr[$i][$id] = "$species[$i].$chr\t$pos\t$len\t$str\t$full";
                    $seq = revcomp($seq) if($rc);
                    $data[$i][$id] = $seq;
		}
	    }
	    else {
    	    	while($line=<FILE>) {
 		    chomp($line);
		    ($id, $chr, $pos, $len, $str, $full, $seq) = split /\t/, $line;
		    if($selected{$id}) {
            	    	$seq =~ tr/[A-Z]/[a-z]/;
	    	    	$attr[$i][$id] = "$species[$i].$chr\t$pos\t$len\t$str\t$full";
		    	$seq = revcomp($seq) if($rc);
	    	    	$data[$i][$id] = $seq;
		    }
    	        }
	    }
	    close FILE;
	    print STDERR "(rc)" if($rc);
    	    print STDERR "]\n";
  	}
    }
}

####################################################################################################################################################

sub read_expression {
    my $expression_file_name = @_[0];
    return unless($expression_file_name);
    open FILE, $expression_file_name || die("Can't open '$expression_file_name'\n");
    while($line=<FILE>) {
    	chomp $line;
	my @arr = split /\s+/, $line;
	my $key = shift(@arr);
	if($key eq "data") {
	    my @files=();
	    foreach my $tag(@arr) {
		$id = $1 if($tag =~ /id=(.+)/);
		push @files, $1 if($tag =~ /file=(.+)/);
	    }
	    next unless($id);
	    foreach $file(@files) {
		read_datafile($id, $file);
	    }
	}
	if($key eq "plot") {
	    push @PLOTS, $line 
	}
    }
    close FILE;
}

sub read_datafile {
    my $id = @_[0];
    my $file = @_[1];
    $file = "$ENV{HOME}$1" if($file=~/^\~(.*)$/);
    print STDERR "[Data: $file -> $id ";
    open INPUT, $file || die("Can't open '$file'\n");
    my @header = ();
    while($line=<INPUT>) {
        chomp $line;
        my @arr = split /\s+/, $line;
        if(@header) {
            my $key = shift(@arr);
            next unless($key);
            for($i=0;$i<@header;$i++) {
            	$VALUE{$id}{$key}{$header[$i]} = $arr[$i];
            }
	} 
        else {
            @header = @arr;
        }
    }
    close INPUT;
    print STDERR 0+keys(%{$VALUE{$id}}),"]\n";
}
