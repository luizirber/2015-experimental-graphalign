#! /usr/bin/env python
import sys
import argparse
import screed


def output_single(read):
    name = read.name
    sequence = read.sequence

    accuracy = None
    if hasattr(read, 'accuracy'):
        accuracy = read.accuracy
        assert len(sequence) == len(accuracy), (sequence, accuracy)
        return "@%s\n%s\n+\n%s\n" % (name, sequence, accuracy)
    else:
        return ">%s\n%s\n" % (name, sequence)


def read_pos_file(filename):
    for line in open(filename):
        line = line.strip()
        try:
            read, posns = line.split(' ', 2)
            posns = map(int, posns.split(','))
        except ValueError:
            read = line
            posns = []
            
        yield read, posns

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('posfile')
    parser.add_argument('reads')
    parser.add_argument('--limitreads', default=None)
    parser.add_argument('--save-erroneous-to', dest='save_err_to',
                        help="output erroneous reads to this file",
                        type=argparse.FileType('w'), default=None)
    args = parser.parse_args()

    print 'reading files...', args.posfile, args.reads
    posdict = dict(read_pos_file(args.posfile))

    limitnames = None
    if args.limitreads:
        limitnames = set([ readname for readname, _ in \
                           read_pos_file(args.limitreads) ])
    
    all_reads = 0
    sum_bp = 0

    print 'reading sequences...'
    for n, record in enumerate(screed.open(args.reads)):
        if n % 100000 == 0:
            print >>sys.stderr, '...', n

        if args.limitreads and record.name not in limitnames:
            continue

        all_reads += 1
        sum_bp += len(record.sequence)

    print 'done!'

    n_reads = 0
    n = 0
    m = 0
    skipped = 0
    for k, v in posdict.iteritems():
        if args.limitreads and k not in limitnames:
            skipped += 1
            continue

        n_reads += 1

        if not v:
            continue

        n += 1
        m += len(v)

    if args.save_err_to:
        for nn, record in enumerate(screed.open(args.reads)):
            if nn % 100000 == 0:
                print >>sys.stderr, '... x save', nn
            if posdict.get(record.name):
                args.save_err_to.write(output_single(record))

    print 'XXX', all_reads, n_reads

    print 'posfile %s: %d mutated reads of %d; %d mutations total' % \
          (args.posfile, n, n_reads, m)
    print 'skipped:', skipped
    print '%d bp total' % (sum_bp,)
    print 'overall error rate: %f%%' % (100. * m / float(sum_bp))


if __name__ == '__main__':
    main()


