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
print "GOLDENPATH=http://hgdownload.soe.ucsc.edu/goldenPath/\n";

foreach my $key(keys(%UCSC)) {
    my $name = $UCSC{$key};	# ucsc name
    my $Name = $UCSC{$key};	# ucsc name first character capital (for vsMm9, for instance)
    substr($Name,0,1) =~ tr/[a-z]/[A-Z]/;

    next unless($name);

    #Genome download
    print "\$(DOWNLOAD)$key/md5sum.txt : \n";
    print "\tmkdir -p \$(DOWNLOAD)$key\n";
    if($DATATYPE{$key} eq "fa") {
	print "\twget \${GOLDENPATH}$name/bigZips/$name.fa.gz -O \$(DOWNLOAD)$key/$name.fa.gz\n";
	print "\tgunzip -f \$(DOWNLOAD)$key/$name.fa.gz\n";
	print "\twget \${GOLDENPATH}$name/bigZips/$name.fa.masked.gz -O \$(DOWNLOAD)$key/$name.fa.masked.gz\n";
	print "\tgunzip -f \$(DOWNLOAD)$key/$name.fa.masked.gz\n";
    }
    if($DATATYPE{$key} eq "scaffold") {
        print "\twget \${GOLDENPATH}$name/bigZips/scaffoldFa.gz -O \$(DOWNLOAD)$key/scaffold.fa.gz\n";
	print "\tgunzip -f \$(DOWNLOAD)$key/scaffold.fa.gz\n";
        print "\twget \${GOLDENPATH}$name/bigZips/scaffoldFaMasked.gz -O \$(DOWNLOAD)$key/scaffold.fa.masked.gz\n";
        print "\tgunzip -f \$(DOWNLOAD)$key/scaffold.fa.masked.gz\n";
    }
    if($DATATYPE{$key} eq "scaffold-zip") {
        print "\twget \${GOLDENPATH}$name/bigZips/scaffoldFa.zip -O \$(DOWNLOAD)$key/scaffold.fa.zip\n";
        print "\tunzip -f \$(DOWNLOAD)$key/scaffold.fa.zip -d\$(DOWNLOAD)$key/\n";
        print "\twget \${GOLDENPATH}$name/bigZips/scaffoldFaMasked.zip -O \$(DOWNLOAD)$key/scaffold.fa.masked.zip\n";
        print "\tunzip -f \$(DOWNLOAD)$key/scaffold.fa.masked.zip -d \$(DOWNLOAD)$key/\n";
    }
    if($DATATYPE{$key} eq "chrom-zip") {
	print "\twget \${GOLDENPATH}$name/bigZips/chromFa.zip -O \$(DOWNLOAD)$key/chromFa.zip\n";
	print "\tunzip \$(DOWNLOAD)$key/chromFa.zip -d \$(DOWNLOAD)$key/\n";
	print "\twget \${GOLDENPATH}$name/bigZips/chromFaMasked.zip -O \$(DOWNLOAD)$key/chromFaMasked.zip\n";
	print "\tunzip \$(DOWNLOAD)$key/chromFaMasked.zip -d \$(DOWNLOAD)$key/\n";
    }
    if($DATATYPE{$key} eq "chrom" || $DATATYPE{$key} eq "chrom-sub" || $DATATYPE{$key} eq "chrom-soft") {
	print "\twget \${GOLDENPATH}$name/bigZips/chromFa.tar.gz -O \$(DOWNLOAD)$key/chromFa.tar.gz\n";
        print "\tgunzip -f \$(DOWNLOAD)$key/chromFa.tar.gz\n";
	print "\ttar -xf \$(DOWNLOAD)$key/chromFa.tar -C \$(DOWNLOAD)$key/\n";
        print "\twget \${GOLDENPATH}$name/bigZips/chromFaMasked.tar.gz -O \$(DOWNLOAD)$key/chromFaMasked.tar.gz\n";
        print "\tgunzip -f \$(DOWNLOAD)$key/chromFaMasked.tar.gz\n";
        print "\ttar -xf \$(DOWNLOAD)$key/chromFaMasked.tar -C \$(DOWNLOAD)$key/\n";
	if($DATATYPE{$key} eq "chrom-sub") {
            print "\tmv \$(DOWNLOAD)$key/?/*  \$(DOWNLOAD)$key/\n";
            print "\tmv \$(DOWNLOAD)$key/??/* \$(DOWNLOAD)$key/\n";
            print "\trm -f -r \$(DOWNLOAD)$key/? \$(DOWNLOAD)$key/??\n";
	}
        if($DATATYPE{$key} eq "chrom-soft") {
            print "\tmv \$(DOWNLOAD)$key/softMask/* \$(DOWNLOAD)$key/\n";
            print "\tmv \$(DOWNLOAD)$key/hardMask/* \$(DOWNLOAD)$key/\n";
            print "\trm -f -r \$(DOWNLOAD)$key/softMask/ \$(DOWNLOAD)$key/hardMask/\n";
        }

	print "\trm -f \$(DOWNLOAD)$key/chromFa.tar \$(DOWNLOAD)$key/chromFaMasked.tar\n";
    }
    print "\twget \${GOLDENPATH}$name/bigZips/md5sum.txt -O \$(DOWNLOAD)$key/md5sum.txt\n";
    print "genome :: \$(DOWNLOAD)$key/md5sum.txt\n\n";

    unless($key eq $basespecies{$CLADE{$key}}) {	
	print "\$(NETCHAIN)$REFDB{$key}.$key.all.chain:\n";
	print "\tmkdir -p \$(NETCHAIN)\n";
	print "\twget \${GOLDENPATH}$REFDB{$key}/vs$Name/$REFDB{$key}.$name.all.chain.gz -O \$(NETCHAIN)$REFDB{$key}.$key.all.chain.gz\n";
	print "\tgunzip -f \$(NETCHAIN)$REFDB{$key}.$key.all.chain.gz\n";
	print "chain :: \$(NETCHAIN)$REFDB{$key}.$key.net\n\n";

        print "\$(NETCHAIN)$REFDB{$key}.$key.net:\n";
        print "\tmkdir -p \$(NETCHAIN)\n";
        print "\twget \${GOLDENPATH}$REFDB{$key}/vs$Name/$REFDB{$key}.$name.net.gz -O \$(NETCHAIN)$REFDB{$key}.$key.net.gz\n";
        print "\tgunzip -f \$(NETCHAIN)$REFDB{$key}.$key.net.gz\n\n";
	print "net :: \$(NETCHAIN)$REFDB{$key}.$key.net\n\n";
    }
}
