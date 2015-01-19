# evaluate against perfect?

all: compare.txt

ecoli-subset.fq:
	sample-reads-randomly.py -R 1 -N 10000 ecoli-mapped.fq.gz.keep.gz -o ecoli-subset.fq

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
	load-into-counting.py -k 21 -x 4e7 ecoli.dn.k21.kh ecoli-mapped.fq.gz.keep.gz

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
