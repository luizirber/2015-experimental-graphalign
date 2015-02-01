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

    ct = khmer.load_counting_hash(args.table)
    aligner = khmer.new_readaligner(ct, 5, 1.0)

    for record in screed.open(args.ref):
        seq = record.sequence
        seq = seq.replace('N', 'A')

        score, alignment = align_long(ct, aligner, seq)

        g = alignment.g
        r = alignment.r

        m, n = alignment.compare()
        print record.name, m, n, n - m, "%.3f%%" % (float(m)/ n * 100)
        #for start in range(0, len(alignment), 60):
        #    print alignment[start:start+60]

        if 0:
            gidx = AlignmentIndex(alignment)
            for gi, a, b in alignment.mismatches():
                print gi, a, b, gidx.get_ri(gi)

        if 0:
            print len(seq), alignment.refseqlen()
            gidx._sanityCheck(seq)

if __name__ == '__main__':
    main()
