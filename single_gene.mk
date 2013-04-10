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

include metacalc.mk

${OUTDIR}${SPECIES}/${NAME}.met:        ${METADATA}${DOMAIN}/${NAME}.cps
	${XDIR}trim -i ${SPECIES}.cfg -r ${METADATA}${DOMAIN}/${NAME}.cps -o ${OUTDIR}${SPECIES}/${NAME}.met -maf ${METADATA}${DOMAIN}/${NAME}.mus

${OUTDIR}${SPECIES}/${NAME}.tab: ${OUTDIR}${SPECIES}/${NAME}.met ${OUTDIR}${SPECIES}/${NAME}.met
	${XDIR}irbis -l ${OUTDIR}${SPECIES}/${NAME}.met -r ${OUTDIR}${SPECIES}/${NAME}.met -o ${OUTDIR}${SPECIES}/${NAME}.tab ${PARAMS}

${OUTDIR}${SPECIES}/${NAME}.bed: ${OUTDIR}${SPECIES}/${NAME}.tab ${SPECIES}.cfg
	perl ${PDIR}tab2bed.pl -i ${OUTDIR}${SPECIES}/${NAME}.tab -s ${SPECIES}.cfg -o ${OUTDIR}${SPECIES}/${NAME}.bed

${OUTDIR}${SPECIES}/${NAME}.maf: ${SPECIES}.cfg ${OUTDIR}${SPECIES}/${NAME}.tab
	perl ${PDIR}tab2maf.pl -l ${SPECIES}.cfg -r ${SPECIES}.cfg -i ${OUTDIR}${SPECIES}/${NAME}.tab -o ${OUTDIR}${SPECIES}/${NAME}.maf -maf ${METADATA}${DOMAIN}/${NAME}.mus ${PARAMS}

${OUTDIR}${SPECIES}/${NAME}.pdf: ${SPECIES}.cfg ${OUTDIR}${SPECIES}/${NAME}.tab ${OUTDIR}${SPECIES}/${NAME}.maf
	perl ${PDIR}tab2pdf.pl -l ${SPECIES}.cfg -r ${SPECIES}.cfg -i ${OUTDIR}${SPECIES}/${NAME}.tab -m ${OUTDIR}${SPECIES}/${NAME}.maf -o ${OUTDIR}${SPECIES}/${NAME}.pdf ${PARAMS}

tab:    ${OUTDIR}${SPECIES}/${NAME}.tab

bed:	${OUTDIR}${SPECIES}/${NAME}.bed

maf:    ${OUTDIR}${SPECIES}/${NAME}.maf

pdf:    ${OUTDIR}${SPECIES}/${NAME}.pdf

pdf-:   
	rm -f ${OUTDIR}${SPECIES}/${NAME}.pdf ${OUTDIR}${SPECIES}/${NAME}.aux ${OUTDIR}${SPECIES}/${NAME}.tex ${OUTDIR}${SPECIES}/${NAME}.log ${OUTDIR}${SPECIES}/${NAME}.out
maf-:   
	rm -f ${OUTDIR}${SPECIES}/${NAME}.maf
tab-:
	rm -f ${OUTDIR}${SPECIES}/${NAME}.tab
