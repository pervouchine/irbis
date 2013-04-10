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

print "include database.mk
.PHONY: all clean\n\n";

do_cmd('insect',	  'insect',  ["snoRNA", "snRNA", "ncRNA", "cspcg", "ncpcg"],  {"cspcg"=>"ncpcg"}, "-t .75 -L 12 -g 1 ");
do_cmd('nematode',  'nematode',["snoRNA", "snRNA", "ncRNA", "cspcg", "ncpcg"],  {"cspcg"=>"ncpcg"}, "-t .60 -L 12 -g 1 ");
do_cmd('vertebrate','mammal',  ["snoRNA", "snRNA", "lncRNA", "cspcg", "ncpcg"], {"cspcg"=>"ncpcg"}, "-t .80 -L 12 -g 1 ");


sub do_cmd {
    my $domain = @_[0];
    my $clade  = @_[1];
    my @biotypes = @{@_[2]};
    my %donottouch = %{@_[3]};
    my $params   = @_[4];

    my $db = $REFDB{$BASESPECIES{$domain}};

    foreach $biotype(@biotypes) {
	make(script=>"\${XDIR}trim", input=>{-i=>"$clade.cfg",-r=>"\${METADATA}$domain/$biotype.cps"}, output=>{-o=>"\${OUTDIR}$clade/$biotype.met"}, group=>"$clade");
	make(script=>"\${XDIR}trim", input=>{-i=>"$clade.cfg",-r=>"\${METADATA}$domain/$biotype.cps"}, output=>{-o=>"\${OUTDIR}$clade/$biotype\_rc.met"}, after=>'-rc');
    }

    for($i=0;$i<@biotypes;$i++) {
    	for($j=$i+1;$j<@biotypes;$j++) {
	    next if($donottouch{$biotypes[$i]} eq $biotypes[$j]);
	    make(script=>"\${XDIR}irbis", input=>{-l=>"\${OUTDIR}$clade/$biotypes[$i].met", -r=>"\${OUTDIR}$clade/$biotypes[$j].met"},
					  output=>{-o=>"\${OUTDIR}$clade/$biotypes[$i]\_$biotypes[$j].tab"}, after=>"$params -v", group=>"$clade");
	    make(script=>"\${XDIR}irbis", input=>{-l=>"\${OUTDIR}$clade/$biotypes[$i].met", -r=>"\${OUTDIR}$clade/$biotypes[$j]\_rc.met"},
					  output=>{-o=>"\${OUTDIR}$clade/$biotypes[$i]\_$biotypes[$j]\_rc.tab"}, after=>"$params -v", group=>"$clade");
    	}
    }

    foreach $biotype(@biotypes) {
	make(script=>"\${XDIR}irbis",	   input=>{-l=>"\${OUTDIR}$clade/$biotype.met",-r=>"\${OUTDIR}$clade/$biotype.met",-B=>"\${METADATA}$domain/$db.a2a"},
					   output=>{-o=>"\${OUTDIR}$clade/$biotype\_SS.tab"}, after=>"$params -u 1");
	make(script=>"\${PDIR}tab2maf.pl", input=>{-l=>"$clade.cfg",-r=>"$clade.cfg",-i=>"\${OUTDIR}$clade/$biotype\_SS.tab"},
					   output=>{-o=>"\${OUTDIR}$clade/$biotype\_SS.maf"}, after=>"$params");
	make(script=>"\${PDIR}tab2pdf.pl", input=>{-l=>"$clade.cfg",-r=>"$clade.cfg",-i=>"\${OUTDIR}$clade/$biotype\_SS.tab",-m=>"\${OUTDIR}$clade/$biotype\_SS.maf"},
					   output=>{-o=>"\${OUTDIR}$clade/$biotype\_SS.pdf"}, depends=>"\${METADATA}$domain/$db.sgn", group=>"$clade-ss");
    }
    print "all :: $clade\nss :: $clade-ss\n";
}
