# evaluate against perfect?
QUAKE=/Users/t/Documents/papers/2014-streaming/Quake

all: ref-errors.details.txt corr-errors.details.txt

display: compare.txt refcompare.txt 

ecoli-subset.fq: ecoli-mapped.fq.gz
	sample-reads-randomly.py -R 1 -N 10000 ecoli-mapped.fq.gz -o ecoli-subset.fq

ecoli-subset.100k.fq: ecoli-mapped.fq.gz
	sample-reads-randomly.py -R 2 -N 100000 ecoli-mapped.fq.gz -o ecoli-subset.100k.fq

ecoli-subset.100k.cor.fq.gz: ecoli-subset.100k.fq
	echo ecoli-subset.100k.fq > ecoli_quake_list.txt
	#ln -fs /usr/local/bin/jellyfish .
	#cat ecoli-mapped.fq.gz.keep | $(QUAKE)/bin/count-qmers -q 33 -k 14 > ecoli_dn_counts.out
	#$(QUAKE)/bin/cov_model.py ecoli_dn_counts.out > ecoli_dn_counts.cov
	$(QUAKE)/bin/correct -f ecoli_quake_list.txt -p 4 -k 14 -q 33 -c 7.94 -z -m ecoli_dn_counts.out
	#mv ecoli-mapped.cor.fq.gz ecoli-mapped.dn.cor.fq.gz

ecoli-subset.100k.quake.fq: ecoli-subset.100k.cor.fq.gz ecoli-subset.100k.fq
	extract-original-reads-from-quake-cor.py ecoli-subset.100k.fq ecoli-subset.100k.cor.fq.gz ecoli-subset.100k.quake.fq
	

original.sam: ecoli-subset.fq
	cat ecoli-subset.fq | bowtie2 -p 4 -x ecoli -U - -S original.sam

original.sam.pos: original.sam
	./sam-scan.py ecoliMG1655.fa original.sam -o original.sam.pos

###

ecoli.dn.k13.kh: ecoli-mapped.fq.gz.keep.gz
	load-into-counting.py -k 13 -x 4e7 ecoli.dn.k13.kh ecoli-mapped.fq.gz.keep.gz

corr.k13.C17.fq: ecoli.dn.k13.kh ecoli-subset.fq
	../sandbox/error-correct-pass2.py ecoli.dn.k13.kh ecoli-subset.fq -o corr.k13.C17.fq --trusted-cov 17

corr.k13.C17.sam: corr.k13.C17.fq
	cat corr.k13.C17.fq | bowtie2 -p 4 -x ecoli -U - -S corr.k13.C17.sam

corr.k13.C17.sam.pos: corr.k13.C17.sam
	./sam-scan.py ecoliMG1655.fa corr.k13.C17.sam -o corr.k13.C17.sam.pos

###

ecoli.dn.k15.kh: ecoli-mapped.fq.gz.keep.gz
	load-into-counting.py -k 15 -x 4e7 ecoli.dn.k15.kh ecoli-mapped.fq.gz.keep.gz

corr.k15.C15.fq: ecoli.dn.k15.kh ecoli-subset.fq
	../sandbox/error-correct-pass2.py ecoli.dn.k15.kh ecoli-subset.fq -o corr.k15.C15.fq --trusted-cov 15

corr.k15.C15.sam: corr.k15.C15.fq
	cat corr.k15.C15.fq | bowtie2 -p 4 -x ecoli -U - -S corr.k15.C15.sam

corr.k15.C15.sam.pos: corr.k15.C15.sam
	./sam-scan.py ecoliMG1655.fa corr.k15.C15.sam -o corr.k15.C15.sam.pos

###

ecoli.dn.k17.kh: ecoli-mapped.fq.gz.keep.gz
	load-into-counting.py -k 17 -x 4e7 ecoli.dn.k17.kh ecoli-mapped.fq.gz.keep.gz

corr.k17.C15.fq: ecoli.dn.k17.kh ecoli-subset.fq
	../sandbox/error-correct-pass2.py ecoli.dn.k17.kh ecoli-subset.fq -o corr.k17.C15.fq --trusted-cov 15

corr.k17.C15.sam: corr.k17.C15.fq
	cat corr.k17.C15.fq | bowtie2 -p 4 -x ecoli -U - -S corr.k17.C15.sam

corr.k17.C15.sam.pos: corr.k17.C15.sam
	./sam-scan.py ecoliMG1655.fa corr.k17.C15.sam -o corr.k17.C15.sam.pos

###

ecoli.dn.k21.kh: ecoli-mapped.fq.gz.keep.gz
	load-into-counting.py -k 21 -x 8e7 ecoli.dn.k21.kh ecoli-mapped.fq.gz.keep.gz

ecoli.dn.k23.kh: ecoli-mapped.fq.gz.keep.gz
	load-into-counting.py -k 23 -x 8e7 ecoli.dn.k23.kh ecoli-mapped.fq.gz.keep.gz

ecoli.dn.k31.kh: ecoli-mapped.fq.gz.keep.gz
	load-into-counting.py -k 31 -x 8e7 ecoli.dn.k31.kh ecoli-mapped.fq.gz.keep.gz

corr.k21.C5.fq: ecoli.dn.k21.kh ecoli-subset.fq
	../sandbox/error-correct-pass2.py ecoli.dn.k21.kh ecoli-subset.fq -o corr.k21.C5.fq --trusted-cov 5

corr.k21.C5.sam: corr.k21.C5.fq
	cat corr.k21.C5.fq | bowtie2 -p 4 -x ecoli -U - -S corr.k21.C5.sam

corr.k21.C5.sam.pos: corr.k21.C5.sam
	./sam-scan.py ecoliMG1655.fa corr.k21.C5.sam -o corr.k21.C5.sam.pos

ecoli-ref.dn.k21.kh: ecoliMG1655.fa
	load-into-counting.py -k 21 -x 1e7 ecoli-ref.dn.k21.kh ecoliMG1655.fa

refcorr.k21.fq: ecoli-ref.dn.k21.kh
	../sandbox/error-correct-pass2.py ecoli-ref.dn.k21.kh ecoli-subset.fq -o refcorr.k21.fq --trusted-cov=1

refcorr.k21.sam: refcorr.k21.fq
	cat refcorr.k21.fq | bowtie2 -p 4 -x ecoli -U - -S refcorr.k21.sam

refcorr.k21.sam.pos: refcorr.k21.sam
	./sam-scan.py ecoliMG1655.fa refcorr.k21.sam -o refcorr.k21.sam.pos

###

ref-errors.sam: ref-errors.fq
	bowtie2 -p 4 -x ecoli -U ref-errors.fq -S ref-errors.sam

ref-errors.details.txt: ref-errors.sam
	./sam-scan-details.py ecoliMG1655.fa ref-errors.sam ecoli-ref.dn.k21.kh -o ref-errors.details.txt -C 1

corr-errors.sam: corr-errors.fq
	bowtie2 -p 4 -x ecoli -U corr-errors.fq -S corr-errors.sam

corr-errors.details.txt: corr-errors.sam
	./sam-scan-details.py ecoliMG1655.fa corr-errors.sam ecoli.dn.k21.kh -o corr-errors.details.txt -C 5

###

test-set.sam: test-set.fq
	bowtie2 -p 4 -x ecoli -U test-set.fq -S test-set.sam

test-set.sam.pos: test-set.sam
	./sam-scan.py ecoliMG1655.fa test-set.sam -o test-set.sam.pos

test-set-report: test-set.sam.pos
	./summarize-pos-file.py test-set.sam.pos test-set.fq

###


compare.txt: original.sam.pos corr.k13.C17.sam.pos corr.k15.C15.sam.pos \
	corr.k17.C15.sam.pos corr.k21.C5.sam.pos
	./summarize-pos-file.py original.sam.pos ecoli-subset.fq
	./summarize-pos-file.py corr.k13.C17.sam.pos corr.k13.C17.fq
	./summarize-pos-file.py corr.k15.C15.sam.pos corr.k15.C15.fq
	./summarize-pos-file.py corr.k17.C15.sam.pos corr.k17.C15.fq
	./summarize-pos-file.py corr.k21.C5.sam.pos corr.k21.C5.fq --save-erroneous-to=corr-errors.fq

refcompare.txt: refcorr.k21.sam.pos
	./summarize-pos-file.py original.sam.pos ecoli-subset.fq
	./summarize-pos-file.py refcorr.k21.sam.pos refcorr.k21.fq --save-erroneous-to=ref-errors.fq

ecoli.1.bt2: ecoliMG1655.fa
	bowtie2-build ecoliMG1655.fa ecoli
	samtools faidx ecoliMG1655.fa

#####

NULLGRAPH=~/dev/nullgraph

simple-haplo-reads.fa: simple-haplo.fa
	$(NULLGRAPH)/make-reads.py -C 100 -S 1 simple-haplo.fa > simple-haplo-reads.fa

simple-reads.fa: simple.fa
	$(NULLGRAPH)/make-biased-reads.py -S 1 -C 100 simple.fa > simple-reads.fa

simple-reads.ct: simple-reads.fa
	load-into-counting.py -k 20 -x 1.1e6 simple-reads.ct simple-reads.fa

simple-haplo-reads.dn.ct: simple-haplo-reads.fa
	normalize-by-median.py -k 20 -C 20 simple-haplo-reads.fa -s simple-haplo-reads.dn.ct -x 1.1e6

ex1: simple-haplo-reads.dn.ct
	find-variant-by-align.py simple-haplo-reads.dn.ct simple-alt.fa

ex2: simple-reads.ct
	count-by-align.py simple-reads.ct simple-haplo.fa > simple-haplo.counts
	count-by-align.py simple-reads.ct simple-alt.fa > simple-alt.counts

ex3:simple-haplo-reads.dn.ct
	find-variant-by-align.py simple-haplo-reads.dn.ct simple-long.fa