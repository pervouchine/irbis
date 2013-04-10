all :: benchmark/test3_ce6.pdf benchmark/test3_dm3.pdf benchmark/test3_hg19.pdf

benchmark/test1.res : benchmark/test1.pl
	perl benchmark/test1.pl > benchmark/test1.res

benchmark/test2.res : benchmark/test2.pl
	perl benchmark/test2.pl > benchmark/test2.res

benchmark/test3_ce6.res : ~/db/metadata/nematode/ce6.cps ~/db/output/nematode/ncpcg.met 
	perl benchmark/test3.pl  ~/db/metadata/nematode/ ~/db/output/nematode/ ce6 nematode.cfg > benchmark/test3_ce6.res

benchmark/test3_ce6.pdf : benchmark/test3_ce6.res
	R --no-save --args benchmark/test3_ce6.res benchmark/test3_ce6.pdf < benchmark/plot1.r

benchmark/test3_dm3.res : ~/db/metadata/insect/dm3.cps ~/db/output/insect/ncpcg.met
	perl benchmark/test3.pl  ~/db/metadata/insect/ ~/db/output/insect/ dm3 insect.cfg > benchmark/test3_dm3.res

benchmark/test3_dm3.pdf : benchmark/test3_dm3.res
	R --no-save --args benchmark/test3_dm3.res benchmark/test3_dm3.pdf < benchmark/plot1.r

benchmark/test3_hg19.res : ~/db/metadata/vertebrate/hg19.cps ~/db/output/mammal/ncpcg.met
	perl benchmark/test3.pl  ~/db/metadata/vertebrate/ ~/db/output/mammal/ hg19 mammal.cfg > benchmark/test3_hg19.res

benchmark/test3_hg19.pdf : benchmark/test3_hg19.res
	R --no-save --args benchmark/test3_hg19.res benchmark/test3_hg19.pdf < benchmark/plot1.r



