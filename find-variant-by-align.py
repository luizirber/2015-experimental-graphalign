#! /usr/bin/env python
import sys
import argparse
import screed
import khmer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('table')
    parser.add_argument('ref')
    args = parser.parse_args()

    ct = khmer.load_counting_hash(args.table)
    aligner = khmer.ReadAligner(ct, 5, 1.0)

    for record in screed.open(args.ref):
        s = record.sequence
        s = s.replace('N', 'A')

        score, graph_alignment, read_alignment, truncated = \
               aligner.align(s)

        #assert not truncated

        g = graph_alignment #.replace('-', '')
        r = read_alignment  #.replace('-', '')

        line1 = []
        line2 = []
        line3 = []
        for n, (a, b) in enumerate(zip(g, r)):
            line1.append(a)
            line3.append(b)
            if a != b:
                line2.append(' ')
            else:
                line2.append('|')

        print '::', record.name, score, truncated
        for start in range(0, len(line1), 60):
            print "".join(line1[start:start+60])
            print "".join(line2[start:start+60])
            print "".join(line3[start:start+60])
            print '--'

        #print record.name, ct.find_spectral_error_positions(s, 10)

        #print covs


if __name__ == '__main__':
    main()
