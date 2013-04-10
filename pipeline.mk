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

include metacalc
#OUTDIR comes from metacalc
#METADATA comes from metacalc
#SPECIES <command line>
#DOMAIN  <command line>
#LEFT    <command line>
#RIGHT	 <command line>
#OUT	 <command line>

${OUTDIR}${SPECIES}/${LEFT}.met:	${METADATA}${DOMAIN}/${LEFT}.cps
	./trim -i ${SPECIES}.cfg -r ${METADATA}${DOMAIN}/${LEFT}.cps -o ${OUTDIR}${SPECIES}/${LEFT}.met

${OUTDIR}${SPECIES}/${RIGHT}.met:	${METADATA}${DOMAIN}/${RIGHT}.cps
	./trim -i ${SPECIES}.cfg -r ${METADATA}${DOMAIN}/${RIGHT}.cps -o ${OUTDIR}${SPECIES}/${RIGHT}.met

${OUTDIR}${SPECIES}/${OUT}.tab:	${OUTDIR}${SPECIES}/${LEFT}.met ${OUTDIR}${SPECIES}/${RIGHT}.met
	./irbis -l ${OUTDIR}${SPECIES}/${LEFT}.met -r ${OUTDIR}${SPECIES}/${RIGHT}.met -o ${OUTDIR}${SPECIES}/${OUT}.tab ${PARAMS}

${OUTDIR}${SPECIES}/${OUT}.bed: ${OUTDIR}${SPECIES}/${OUT}.tab ${SPECIES}.cfg
	./_tab2bed -s ${SPECIES}.cfg -i ${OUTDIR}${SPECIES}/${OUT}.tab -o ${OUTDIR}${SPECIES}/${OUT}.bed

${OUTDIR}${SPECIES}/${OUT}.maf: ${SPECIES}.cfg ${OUTDIR}${SPECIES}/${OUT}.tab
	./_tab2maf -l ${SPECIES}.cfg -r ${SPECIES}.cfg -i ${OUTDIR}${SPECIES}/${OUT}.tab -o ${OUTDIR}${SPECIES}/${OUT}.maf ${PARAMS}

${OUTDIR}${SPECIES}/${OUT}.pdf: ${SPECIES}.cfg ${OUTDIR}${SPECIES}/${OUT}.tab ${OUTDIR}${SPECIES}/${OUT}.maf
	./_tab2pdf -l ${SPECIES}.cfg -r ${SPECIES}.cfg -i ${OUTDIR}${SPECIES}/${OUT}.tab -m ${OUTDIR}${SPECIES}/${OUT}.maf -o ${OUTDIR}${SPECIES}/${OUT}.pdf -notitle ${PARAMS}

view:	${OUTDIR}${SPECIES}/${OUT}.pdf
	acroread ${OUTDIR}${SPECIES}/${OUT}.pdf

tab:	${OUTDIR}${SPECIES}/${OUT}.tab

bed:	${OUTDIR}${SPECIES}/${OUT}.bed

maf:	${OUTDIR}${SPECIES}/${OUT}.maf

pdf:	${OUTDIR}${SPECIES}/${OUT}.pdf

pdf-:	
	rm -f ${OUTDIR}${SPECIES}/${OUT}.pdf ${OUTDIR}${SPECIES}/${OUT}.aux ${OUTDIR}${SPECIES}/${OUT}.tex ${OUTDIR}${SPECIES}/${OUT}.log ${OUTDIR}${SPECIES}/${OUT}.out
maf-:	
	rm -f ${OUTDIR}${SPECIES}/${OUT}.maf
tab-:
	rm -f ${OUTDIR}${SPECIES}/${OUT}.tab
