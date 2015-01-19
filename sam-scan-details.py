#! /usr/bin/env python
import sys
import argparse
import screed
import math
import khmer

def ignore_at(iter):
    for item in iter:
        if item.startswith('@'):
            continue
        yield item

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genome')
    parser.add_argument('samfile')
    parser.add_argument('counts_table')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-C', '--cutoff', type=int, default=2)

    args = parser.parse_args()

    kh = khmer.load_counting_hash(args.counts_table)

    genome_dict = dict([ (record.name, record.sequence) for record in \
                        screed.open(args.genome) ])

    n = 0
    n_skipped = 0
    n_rev = n_fwd = 0

    for samline in ignore_at(open(args.samfile)):
        n += 1
        if n % 100000 == 0:
            print >>sys.stderr, '...', n

        readname, flags, refname, refpos, _, _, _, _, _, seq = \
                  samline.split('\t')[:10]
        if refname == '*' or refpos == '*':
            # (don't count these as skipped)
            continue
        
        refpos = int(refpos)
        try:
            ref = genome_dict[refname][refpos-1:refpos+len(seq) - 1]
        except KeyError:
            print >>sys.stderr, "unknown refname: %s; ignoring (read %s)" % (refname, readname)
            n_skipped += 1
            continue

        errors = []

        out1 = []
        out2 = []
        out3 = []
        out4 = []
        out5 = []
        
        for pos, (a, b) in enumerate(zip(ref, seq)):
            if a != b:
                # see http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output - '16' is the revcomp flag.
                if int(flags) & 16:
                    pos = len(seq) - pos - 1
                    n_rev += 1
                else:
                    n_fwd += 1
                errors.append(pos)

            out1.append(a)
            out2.append(b)
            if a != b:
                out3.append('x')
            else:
                out3.append(' ')

            kmer = ref[pos:pos+kh.ksize()]
            if len(kmer) == kh.ksize() and kh.get(kmer) < args.cutoff:
                out4.append('|')
            else:
                out4.append(' ')

            kmer = seq[pos:pos+kh.ksize()]
            if len(kmer) == kh.ksize() and kh.get(kmer) < args.cutoff:
                if kh.get(kmer) == 0:
                    out5.append('0')
                else:
                    out5.append('<')
            else:
                out5.append(' ')

        for start in range(0, len(out1), 60):
            if start == 0:
                print >>args.outfile, '%s positions %s:%s (len %d)' % \
                      (readname, start, min(start+60, len(out1)), len(out1))
            else:
                print >>args.outfile, '%s:%s (len %d)' % \
                      (start, min(start+60, len(out1)), len(out1))
                
            print >>args.outfile, '      ref:', "".join(out1[start:start+60])
            print >>args.outfile, '     read:', "".join(out2[start:start+60])
            print >>args.outfile, ' mismatch:', "".join(out3[start:start+60])
            print >>args.outfile, ' ref kmer:', "".join(out4[start:start+60])
            print >>args.outfile, 'read kmer:', "".join(out5[start:start+60])
            print >>args.outfile, ''
        print >>args.outfile, '-------'


if __name__ == '__main__':
    main()
