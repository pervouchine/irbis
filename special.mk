include config.mk
VERTEBRATE = ${METADATA}vertebrate/
INSECT	   = ${METADATA}insect/
NEMATODE   = ${METADATA}nematode/

${ANNOTATION}Drosophila_melanogaster.BDGP5.25.64.gtf:
	mkdir -p ${ANNOTATION}
	wget ftp://ftp.ensembl.org/pub/release-64/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.25.64.gtf.gz -O ${ANNOTATION}Drosophila_melanogaster.BDGP5.25.64.gtf.gz
	gunzip ${ANNOTATION}Drosophila_melanogaster.BDGP5.25.64.gtf.gz

${ANNOTATION}gencode.v7.annotation.gtf:
	mkdir -p ${ANNOTATION}
	wget ftp://ftp.sanger.ac.uk/pub/gencode/release_7/gencode.v7.annotation.gtf.gz -O ${ANNOTATION}gencode.v7.annotation.gtf.gz
	gunzip ${ANNOTATION}gencode.v7.annotation.gtf.gz

${ANNOTATION}Caenorhabditis_elegans.WS220.65.gtf:
	mkdir -p ${ANNOTATION}
	wget ftp://ftp.ensembl.org/pub/release-65/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WS220.65.gtf.gz -O ${ANNOTATION}Caenorhabditis_elegans.WS220.65.gtf.gz
	gunzip ${ANNOTATION}Caenorhabditis_elegans.WS220.65.gtf.gz

${ANNOTATION}Gencode_lncRNAsv7_summaryTable.txt:
	mkdir -p ${ANNOTATION}
	wget http://genome.crg.es/~rjohnson/lncrna_webpage_data/Gencode_lncRNAsv7_summaryTable_02_14_2012.txt -O ${ANNOTATION}Gencode_lncRNAsv7_summaryTable.txt
	#gunzip ${ANNOTATION}Gencode_lncRNAsv7_summaryTable_10-01-2011.txt.gz

#######################################################################################################################################

${VERTEBRATE}hg19.cps : ${ANNOTATION}gencode.v7.annotation.gtf ${PDIR}gen_loci.pl ${PDIR}gtf2db.pl
	mkdir -p ${VERTEBRATE}
	perl ${PDIR}gen_loci.pl -in ${ANNOTATION}gencode.v7.annotation.gtf > ${VERTEBRATE}hg19.glc
	perl ${PDIR}gtf2db.pl -db ${VERTEBRATE}hg19 -i ${ANNOTATION}gencode.v7.annotation.gtf

${INSECT}dm3.cps : ${ANNOTATION}Drosophila_melanogaster.BDGP5.25.64.gtf
	mkdir -p ${INSECT}
	perl ${PDIR}gtf2db.pl -db ${INSECT}dm3 -i ${ANNOTATION}Drosophila_melanogaster.BDGP5.25.64.gtf

${NEMATODE}ce6.cps:	${ANNOTATION}Caenorhabditis_elegans.WS220.65.gtf
	mkdir -p ${NEMATODE}
	perl ${PDIR}gtf2db.pl -db ${NEMATODE}ce6 -i ${ANNOTATION}Caenorhabditis_elegans.WS220.65.gtf -full

#######################################################################################################################################

${ANNOTATION}Gencode_lncRNAsv7_summaryTable_intergenic.txt: ${ANNOTATION}Gencode_lncRNAsv7_summaryTable.txt
	awk '$$7=="intergenic"' ${ANNOTATION}Gencode_lncRNAsv7_summaryTable.txt > ${ANNOTATION}Gencode_lncRNAsv7_summaryTable_intergenic.txt

${VERTEBRATE}lncRNA.cps: ${ANNOTATION}Gencode_lncRNAsv7_summaryTable_intergenic.txt ${VERTEBRATE}hg19.cps ${PDIR}intergenic_lncrna.pl
	perl ${PDIR}intergenic_lncrna.pl -s ${VERTEBRATE}hg19.cps -i ${ANNOTATION}Gencode_lncRNAsv7_summaryTable_intergenic.txt | sort -u -k 5 > ${VERTEBRATE}lncRNA.cps

