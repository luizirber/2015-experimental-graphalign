#! /usr/bin/env python
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


# read in list of error positions per read
def read_pos_file(filename, ignore_set=set()):
    for line in open(filename):
        line = line.strip()
        try:
            read, posns = line.split(' ', 2)
        except ValueError:
            read = line
            posns = []

        if read in ignore_set:
            posns = []
        elif posns:
            if posns is 'V':            # ignore variable coverage foo
                ignore_set.add(read)
                posns = []
            else:
                posns = list(sorted(map(int, posns.split(','))))

        yield read, posns


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('posfile')
    parser.add_argument('readsfile')
    args = parser.parse_args()

    posdict = dict(read_pos_file(args.posfile))

    out_1err = open(args.readsfile + '.err1', 'w')
    out_2err = open(args.readsfile + '.err2', 'w')
    out_3err = open(args.readsfile + '.err3', 'w')
    out_4err = open(args.readsfile + '.err4', 'w')

    for n, read in enumerate(screed.open(args.readsfile)):
        errors = posdict.get(read.name, [])
        if errors:
            if len(errors) == 1:
                out_1err.write(output_single(read))
            elif len(errors) == 2:
                out_2err.write(output_single(read))
            elif len(errors) == 3:
                out_3err.write(output_single(read))
            elif len(errors) >= 4:
                out_4err.write(output_single(read))

if __name__ == '__main__':
    main()
