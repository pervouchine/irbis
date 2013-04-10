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

include config.mk

CC = g++ -O3
PROGRAMS=${CDIR}map_agnostic ${CDIR}map_single ${CDIR}net_filter ${CDIR}transf ${CDIR}getsegm ${CDIR}getwind ${CDIR}getmuf \
	 ${CDIR}genutils.o ${CDIR}progressbar.o ${CDIR}indexing ${CDIR}syntenic_filter ${CDIR}best_match \
	 ${XDIR}trim ${XDIR}irbis ${XDIR}extend 

.PHONY:	all clean

all ::	${PROGRAMS} database.mk download.mk metacalc.mk

clean ::
	rm -f ${PROGRAMS} database.mk download.mk metacalc.mk ${PDIR}setup.pm

##################################################################################################################

${PDIR}setup.pm:	${PDIR}setup.pl
		perl ${PDIR}setup.pl > ${PDIR}setup.pm

download.mk : config.dat ${PDIR}make_download.pl ${PDIR}utils.pm ${PDIR}setup.pm
	perl ${PDIR}make_download.pl > download.mk

database.mk : config.dat ${PDIR}make_database.pl ${PDIR}utils.pm ${PDIR}setup.pm
	perl ${PDIR}make_database.pl > database.mk

metacalc.mk : config.dat ${PDIR}make_metacalc.pl ${PDIR}utils.pm ${PDIR}setup.pm
	perl ${PDIR}make_metacalc.pl > metacalc.mk

##################################################################################################################

${CDIR}genutils.o : ${CDIR}genutils.c ${CDIR}genutils.h
	${CC} -c ${CDIR}genutils.c -o ${CDIR}genutils.o

${CDIR}map_agnostic : ${CDIR}map_agnostic.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}map_agnostic ${CDIR}map_agnostic.c ${CDIR}genutils.o ${CDIR}progressbar.o

${CDIR}map_single : ${CDIR}map_single.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}map_single ${CDIR}map_single.c ${CDIR}genutils.o ${CDIR}progressbar.o

${CDIR}net_filter : ${CDIR}net_filter.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}net_filter ${CDIR}net_filter.c ${CDIR}genutils.o ${CDIR}progressbar.o

${CDIR}syntenic_filter : ${CDIR}syntenic_filter.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}syntenic_filter ${CDIR}syntenic_filter.c ${CDIR}genutils.o ${CDIR}progressbar.o 

${CDIR}best_match : ${CDIR}best_match.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}best_match ${CDIR}best_match.c ${CDIR}genutils.o ${CDIR}progressbar.o

${CDIR}progressbar.o:  ${CDIR}progressbar.c ${CDIR}progressbar.h
	${CC} -c ${CDIR}progressbar.c -o ${CDIR}progressbar.o

${CDIR}getwind : ${CDIR}getwind.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}getwind ${CDIR}getwind.c ${CDIR}genutils.o ${CDIR}progressbar.o

${CDIR}transf : ${CDIR}transf.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}transf ${CDIR}transf.c ${CDIR}genutils.o ${CDIR}progressbar.o

${CDIR}getsegm : ${CDIR}getsegm.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}getsegm ${CDIR}getsegm.c ${CDIR}genutils.o ${CDIR}progressbar.o

${CDIR}getmuf : ${CDIR}getmuf.c ${CDIR}genutils.o ${CDIR}progressbar.o
	${CC} -o ${CDIR}getmuf ${CDIR}getmuf.c ${CDIR}genutils.o ${CDIR}progressbar.o

${CDIR}indexing : ${CDIR}indexing.c ${CDIR}genutils.o
	${CC} -o ${CDIR}indexing ${CDIR}indexing.c ${CDIR}genutils.o

##################################################################################################################


${XDIR}trim : ${XDIR}trim.c ${XDIR}dictionary.h ${XDIR}subset.h ${XDIR}orderedset.h ${XDIR}conservation.h ${CDIR}genutils.h ${CDIR}genutils.o
	${CC} -o ${XDIR}trim ${CDIR}genutils.o ${XDIR}trim.c -I ${CDIR} -I ${XDIR}

${XDIR}irbis : ${XDIR}irbis.c ${XDIR}dictionary.h ${XDIR}subset.h ${XDIR}orderedset.h ${XDIR}conservation.h ${CDIR}genutils.h ${CDIR}genutils.o ${XDIR}relation.h
	${CC} -o ${XDIR}irbis ${CDIR}genutils.o ${XDIR}irbis.c -I ${CDIR} -I ${XDIR}

${XDIR}extend : ${XDIR}extend.c 
	${CC} -o ${XDIR}extend ${XDIR}extend.c

muscle:
		wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
		tar -xf muscle3.8.31_i86linux64.tar.gz
		rm muscle3.8.31_i86linux64.tar.gz
		mv muscle3.8.31_i86linux64 muscle
