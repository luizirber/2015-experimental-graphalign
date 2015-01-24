#! /usr/bin/env python
import screed
import argparse


def output_single(read):
    if hasattr(read, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.accuracy)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('original')
    parser.add_argument('quake')
    parser.add_argument('output')
    args = parser.parse_args()

    read_info = {}
    for record in screed.open(args.quake):
        read_info[record.name] = len(record.sequence)

    fp = open(args.output, 'w')

    for record in screed.open(args.original):
        if record.name in read_info:
            length = read_info[record.name]
            record.sequence = record.sequence[:length]
            record.accuracy = record.accuracy[:length]
            fp.write(output_single(record))

if __name__ == '__main__':
    main()
