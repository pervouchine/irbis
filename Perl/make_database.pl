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
#

use Perl::utils;
print "include config.mk\n";
print "include download.mk\n";
print "include special.mk\n";
print "SEGPAR  = -margin 8 -limit 5000\nSGNPAR  = -we 10 -wi 10 -cis -all\n\n";

%biotypes = (   'vertebrate'=>[],
                'insect'=>['ncRNA'],
                'nematode'=>['ncRNA']
            );
%genes    = (   'vertebrate'=>["SF1", "HNRNPL", "PTPRC", "DST", "TTN", "KCNMA1","FOXP1","MAP3K4","FGFR3","FGFR1","SNORD115-1","SNORD116-4","RP11-439A17.4","HIST"],
                'insect'=>["Nmnat", "slo", "Ca-alpha1D", "Dscam","14-3-3zeta","Mhc","GluClalpha", "MRP"],
                'nematode'=>["lin-49","slo-1","mrp-1","clp-1"]
            );

foreach my $key(keys(%UCSC)) {
    $clade = $CLADE{$key};
    push @{$sus{$clade}}, "\${METADATA}$clade/$key.sus";
    next unless($key eq $BASESPECIES{$clade});
    foreach $biotype("protein_coding", "snoRNA", "snRNA", @{$biotypes{$clade}}) {
	make(before=>"awk '\$\$10==\"$biotype\"'",input=>{''=>"\${METADATA}$clade/$REFDB{$key}.cps"}, output=>{'>'=>"\${METADATA}$clade/$biotype.cps"}, group=>"$clade");
    }
    make(before=>"awk '\$\$15~/[EAI][5N3]/'", input=>{''=>"\${METADATA}$clade/protein\_coding.cps"}, output=>{'>'=>"\${METADATA}$clade/ncpcg.cps"}, group=>"$clade");
    make(before=>"awk '\$\$15~/[EAI]C/'",     input=>{''=>"\${METADATA}$clade/protein\_coding.cps"}, output=>{'>'=>"\${METADATA}$clade/cspcg.cps"}, group=>"$clade");
}

foreach my $key(keys(%UCSC)) {
    $clade = $CLADE{$key};
    next unless($key eq $BASESPECIES{$clade});
    $files = join(" ",@{$sus{$clade}});
    foreach $gene(@{$genes{$clade}}) {
	make(before=>"awk '\$\$11==\"$gene\"'", input=>{''=>"\${METADATA}$clade/$REFDB{$key}.cps"}, output=>{'>'=>"\${METADATA}$clade/$gene.cps"}, group=>"$clade");
	make(script=>"\${CDIR}getmuf", input=>{-files=>"$files",-in=>"\${METADATA}$clade/$gene.cps"}, output=>{-out=>"\${METADATA}$clade/$gene.mus"}, group=>"$clade");
    }
}

foreach my $name(keys(%UCSC)) {
    $clade = $CLADE{$name};
    next unless($name eq $BASESPECIES{$clade});
    make(before=>"sort -k2,2n", input=>{''=>"\${METADATA}$clade/$REFDB{$name}.cps"}, output=>{'>'=>"\${METADATA}$clade/$REFDB{$name}.cps.srt"});
    make(script=>"\${PDIR}all2all.pl", input=>{-in=>"\${METADATA}$clade/$REFDB{$name}.cps"}, output=>{-out=>"\${METADATA}$clade/$REFDB{$name}.a2a"}, after=>"-bin");
    foreach my $key(@{$CLADELIST{$clade}}) {
	$ref = $REFDB{$key};
	make(script=>"\${CDIR}transf", input=>{-dir=>"\${DOWNLOAD}$key/md5sum.txt"}, output=>{-idx=>"\${SEQUENCE}$clade/$key.idx", -dbx=>"\${SEQUENCE}$clade/$key.dbx"});
        my $map = ($key eq $BASESPECIES{$clade}) ? "\${METADATA}$clade/$ref.cps" : "\${METADATA}$clade/$key.map";
	if($key eq $BASESPECIES{$clade}) {
	    make(script=>"\${CDIR}getwind", input=>{-in=>"\${METADATA}$clade/$ref.cps",-idx=>"\${SEQUENCE}$clade/$key.idx",-dbx=>"\${SEQUENCE}$clade/$key.dbx"},
					    output=>{-out=>"\${METADATA}$clade/$ref.sgn"}, after=>"\${SGNPAR}", script_not_required=>yes);
	}
	else {
	    make(script=>"\${CDIR}map_single", input=>{-in=>"\${METADATA}$clade/$ref.cps.srt",-chain=>"\${NETCHAIN}$ref.$key.all.chain"}, output=>{-out=>"\${METADATA}$clade/$key.aln"});
	    make(script=>"\${CDIR}best_match", input=>{-in=>"\${METADATA}$clade/$ref.cps",-aln=>"\${METADATA}$clade/$key.aln"}, output=>{-out=>"\${METADATA}$clade/$key.map"});
	}

	make(script=>"\${PDIR}out2bed.pl", input=>{-in=>"$map"},output=>{-out=>"\${METADATA}$clade/$key.bed"}, after=>($key eq $BASESPECIES{$clade} ? "-cis" : undef));
	make(script=>"\${CDIR}getsegm", input=>{-in=>"\${METADATA}$clade/$key.bed",-idx=>"\${SEQUENCE}$clade/$key.idx", -dbx=>"\${SEQUENCE}$clade/$key.dbx"},
					output=>{'>'=>"\$(METADATA)$clade/$key.sus"}, after=>"\$(SEGPAR) | sort -k1,1n");
	make(script=>"\${CDIR}indexing", input=>{-in=>"\${METADATA}$clade/$key.sus"}, output=>{-out=>"\${METADATA}$clade/$key.sus.ind"}, group=>"$clade");
    }
}

print "all :: ",join(" ",keys(%BASESPECIES)),"\n";
