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

    for record in screed.open(args.ref):
        s = record.sequence
        s = s.replace('N', 'A')

        for i in range(0, len(s) - ct.ksize() + 1):
            kmer = s[i:i+ct.ksize()]
            print i, ct.get(kmer)
        
        print record.name, ct.find_spectral_error_positions(s, 10)

if __name__ == '__main__':
    main()
