#! /usr/bin/env python
import sys
import argparse
import screed
import khmer
from graphAlignment import align_long, AlignmentIndex


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('table')
    parser.add_argument('ref')
    args = parser.parse_args()

    ct = khmer.new_counting_hash(20, 1.1e6, 4)
    #ct = khmer.load_counting_hash('simple-haplo-reads.dn.ct')
    ct.consume_fasta('simple-haplo-reads.fa.keep')
    #ct = khmer.load_counting_hash(args.table)
    aligner = khmer.new_readaligner(ct, 5, 1.0)

    for record in screed.open(args.ref):
        seq = record.sequence
        seq = seq.replace('N', 'A')

        score, alignment = align_long(ct, aligner, seq)

        g = alignment.g
        r = alignment.r

        m, n = alignment.compare()
        print record.name, m, n, n - m, "%.3f%%" % (float(m)/ n * 100)
        for start in range(0, len(alignment), 60):
            print alignment[start:start+60]

        if 1:
            print len(seq), alignment.refseqlen()

            gidx = AlignmentIndex(alignment)

            for i in range(len(alignment)):
                print i, alignment[i], gidx.get_ri(i), seq[gidx.get_ri(i)]
                if alignment[i][1] in 'ACGTagct':
                    assert alignment[i][1].upper() == seq[gidx.get_ri(i)]

            for j in range(len(seq)):
                gi = gidx.get_gi(j)
                print j, seq[j], gi, alignment[gi], len(seq)
                if alignment[gi][1] in 'ACGTacgt' or 1:
                    assert seq[j] == alignment[gi][1].upper()


if __name__ == '__main__':
    main()
