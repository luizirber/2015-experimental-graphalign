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

        assert not truncated

        g = graph_alignment.replace('-', '')
        r = read_alignment.replace('-', '')

        print record.name
        for kstart in range(0, len(g) - ct.ksize() + 1):
            kmer = g[kstart:kstart+ct.ksize()]
            print kstart, ct.get(kmer)


if __name__ == '__main__':
    main()
