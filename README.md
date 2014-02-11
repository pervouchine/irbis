=============================================================================================================================
Copyright 2011,2012 Dmitri Pervouchine

This file is part of the IRBIS package.

IRBIS package is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

IRBIS package is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with IRBIS package.  If not, see <http://www.gnu.org/licenses/>.

=============================================================================================================================

INSTALLATION INSTRUCTIONS

1) Since you are reading this file, you must have unzipped the archive already
2) Make sure that you have installed
 - g++ compiler
 - perl
 - muscle 
 - latex + tikz + pgfplots (pgfplots is optional)
3) make sure that you don`t have the directory called 'db' in your home directory; it will be replaced
4) To install and download data, type
    make all
    make -f database.mk all
    make -f metacalc.mk all
5) 'make -f database.mk all' downloads all genomes specified in config.dat from UCSC website
   If you already have them, you can save some time by providing files directly
   You may find it useful to look at database and download makefiles
6) 'make -f metacalc.mk all' creates hash tables and does some matching of conserved complementary motifs (see metacalc.mk)

Should you have any questions, please email Dmitri Pervouchine dp@crg.eu
A short versioning track is provided in VERSION file

=============================================================================================================================

GENOMIC COORDINATE CONVENTION

Splice site coordinates refer to to the phosphodisester bond at the 5` position of the nucleotide, i.e.,
    for D/E first nt of the intron (G of GT) or after the end of transcript (D=donor, E=end)
    for A/S fist  nr of the exon (one after AG) or the first nt of the gene (A=acceptor, S=start)

=============================================================================================================================

FILE TYPES

Besides, BED and GTF formats, a few other tab-delimited file types are used.

1) .cps = chromosome-position-strand
 0 chromosome
 1 position
 2 strand

 3 gene index
 4 site index
 5 type (D=Donor, A=Acceptor, S=Start, E=End)

 6 number of nucleotides to the nearest site upstream
 7 number of nucleotides to the nearest site downstream

 8 db gene_id (Fbgn or ENSG)
 9 type or bio_type
10 db gene_name

11 number of transcripts in which the site is used
12 number of transcripts in covering the site 

13 symbolic notation for the type of the segment [site - 1, site] (see below)
14 symbolic notation for the type of the segment [site, site + 1] (see below)

15 Shortest upstream exon   relative to [site, site + 1] segment
16 Shortest downstream exon relative to [site, site + 1] segment

Segment type consists of two letters. The first can be "E","I","A" or "X". The second can be "C" or "N".
E = always exonic
A = alternative (can be exonic or intronic)
I = always intronic
X = enigmatic (intergenic within a gene) ==> to be resolved with Gencode people

N = doesn`t overlap with coding parts
C = overlaps with at least one CDS

.aln = chromosome-position-strand-chromosome-position-strand
 0 chromosome	FROM
 1 position	FROM
 2 strand	FROM

 3 chromosome	TO
 4 position	TO
 5 strand	TO

 6 gene index
 7 site index
 8 type		(D=Donor, A=Acceptor, S=Start, E=End)

 9 number of nucleotides upstream before the nearest site in the same gene (0 if S)
10 number of nucleotides downstream before the nearest site in the same gene (0 if E)

11 chain alignment score

out = same as aln
However, note that .aln maps FROM=>possibly many TOs, while .out maps FROM=>exactly one TO
See best_match.c for details on how the one TO is selected

.bed = begin-end (same as in bed-tols, UCSC)
0 = chromosome
1 = begin
2 = end
3 = score
4 = name => used for the indexing, smaller of the two indices is taken for each segment when bed files are generated from out files
5 = strand

.tsq = tabulated sequence segment file
0 = index
1 = sequence
[deprecated after 2.4]

.tsw = tabulated sequence window file
0 = index
1 = sequence
[deprecated after 2.4]

.sus = single unaligned segment; .suw = single unaligned window ; .suf = sigle unaligned file
0 = index
1 = chomosome
2 = start
3 = size
4 = strand
5 = srcSize (as in MAF)
6 = sequence

.mus = multiple unaligned segment; .muw = multiple unaligned window ; .muf = multiple unaligned file
LINE-based format same as MAF
a id = index
s specie.chromosome start size strand srcSize ACGATCGGAGCAGCA
ends with \n\n

.int = intron file
0 = gene id
1 = donor id
2 = acceptor id
3 = comma-separated list of transcript ids in which the intron is used

.exn = exon file
0 = gene id
1 = acceptor id
2 = donor id
3 = comma-separated list of transcript ids in which the exon is present
4 = comma-separated list of transcript ids that cover given exon 

===============================================================================================================================================

BRIEF DESCRIPTION OF THE ALGORITHM (might be outdated)

1. There are two sets of packs orthologous sequences called LEFT and RIGHT sets. Packs of orthologous sequences are, for instance, orthologous exons
   or introns, such that no MAS is produced for these sets, but it is known a priori that they are orthologous.
2. For each sequence in the LEFT set in each species i, construct a dictionary (gapped version)
        D(i): gapped n-mer |-> list of (id, pos, gap) or any other instance of the class ordered_set where it occurs
        instances of ordered_set must go in ascending order wrt the "<" operation defined for ordered_set
   Similarly, for each sequence in the RIGHT set in each specie, construct a dictionary
        D*(i): gapped n-mer |-> list of (id, pos, gap) or any other instance of the class ordered_set where occurs its reverse complement

    At that, the convention is that for the LEFT set the positions start from the leftl for the RIGHT set the positions start from the right.
        XXXXggXXXX
                 ^
                 i <- each position in the left SET points to the last nt of the n-mer;

        YYYYggYYYY
        ^
        j <- for the RIGHT set it`s symmetirc
   That is, each dictionary D(i) appoints to a word the set of locations, where the word occurs. The locations are of the form (id, pos, gap) and are listed
   in ascending order for each word. For the D*(i) everything is symmetric, i.e., reverse complements are considered for each word and triplets are listed in
   the reverse order.

3. These dictionaries are TRIMMED as follows:
        a) for each n-mer look up all its ordered_set instances in species i
        b) they have to be sorted in acsending order!
        c) set the array of pointers P(i) at the first elements of the list of ordered_set in species i
        d) select the min element of {P(i)}
        e) look at how many P(i) the min is achieved, where ordered_set instances are equivalent by .equiv() relationship defined in ordered_set
           basically .equiv() means that the ids are the same and positions are not too different
           not too different is defined withon ordered_set class by SAMEWINDOW parameter
        f) if the min is achieved at sufficiently many elements (or more intricate way is to sum weights w(i)) then those ordered_set instances are kept
           sufficiently many is defined by THRESHOLD
        e) increment pointers at which min was achieved and repeat
    By the end of the day, the dictionaries will contain only gapped n-mers which are conserved in at least THRESHOLD species

4.  The next step is taking a direct product. For each word in dictionary D(i) and the same word in dictionary D*(i) we create a list of ordered pairs
    of elements of D(i) and D*(i), i.e. (id1, pos1, gap1, id2, pos2, gap2) sorted in lexicographic order. This direct product of two dictionaries is
    subject to intersection procedure as in 3 (THRESHOLD might have rised).

5.  Next, there is a customary folding procedure, in which each of the conserved complementary words is split into base pairs and those are aligned by
    using max bp approach.


===============================================================================================================================================
