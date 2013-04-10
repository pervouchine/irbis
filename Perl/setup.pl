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


    $READKEY  = 'yes' if(`perl -e 'use Term::ReadKey;print 1;' 2>>/dev/null`);
    if($READKEY) {
        print "use Term::ReadKey;\n\$READKEY=1;\n";
    	print "(\$WCHAR) = GetTerminalSize();\n";
	print STDERR "Found:\tTerm::ReadKey\n";
    }
    else {
        print "\$WCHAR = 192;\n"; 
    }

    `echo \\\\documentclass{article}\\\\usepackage{pgfplots}\\\\begin{document}1\\\\end{document}>tmp$$.tex\npdflatex -interaction nonstopmode tmp$$ 2>/dev/null`;
    if(-e "tmp$$.pdf") {
	print "\$PGFPLOTS = 1;\n";
	print STDERR "Found:\tPGFPLOTS\n";
    }
    `rm -f tmp$$.*`;

    `muscle 2>>tmp$$`;
    if(`cat tmp$$`=~/Edgar/) {
	print STDERR "Found:\tMuscle\n";
    }
    `rm -f tmp$$`;


    print "return(1);\n";
