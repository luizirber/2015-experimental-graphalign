#! /usr/bin/env python
import argparse
import khmer
import screed
import sys
import graphAlignment

REGIONSIZE=100


def turn_locations_into_regions(refposns):
    d = {}
    for ref, pos in refposns:
        x = d.get(ref, set())
        x.add(pos)
        d[ref] = x

    regions = []
    for k in d.keys():
        posns = list(d[k])
        posns.sort()

        start = posns[0]
        i = 0
        while 1:
            while i < len(posns) - 1 and \
                  (posns[i+1] - posns[i]) < REGIONSIZE:
                i += 1

            # done with region - record.
            end = posns[i]
            regions.append((k, start, end))

            # done over all? exit.
            if i == len(posns) - 1:
                break

            # move on to next region.
            i += 1
            start = posns[i]

    return regions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference')
    parser.add_argument('readfile')
    args = parser.parse_args()

    ct = khmer.new_counting_hash(21, 1e7, 4)

    tags_to_positions = {}
    references = {}
    for record in screed.open(args.reference):
        # store for later retrieval - in memory, for now.
        references[record.name] = record.sequence

        # load into graph & tag
        ct.consume_and_tag(record.sequence)

        # track positions in reference by tag
        tagposns = ct.get_tags_and_positions(record.sequence)
        for pos, tag in tagposns:
            x = tags_to_positions.get(tag, [])
            x.append((record.name, pos))
            tags_to_positions[tag] = x

    # now, walk through the reads and map to graph
    aligner = khmer.ReadAligner(ct, 0, 1.0)
    for read in screed.open(args.readfile):

        # align to graph, where possible
        readseq = read.sequence.replace('N', 'A')
        score, g, r, truncated = aligner.align(readseq)
        if truncated:
            print >>sys.stderr, "IGNORING read", read.name
            continue

        # find locations in reference where read alignment overlaps a tag
        refseq = g.replace('-', '')
        ptags = ct.get_tags_and_positions(refseq)
        assert len(ptags)

        refposns = []
        for pos, tag in ptags:
            refposns.extend(tags_to_positions[tag])

        # extract the larger region, remap read to get exact positions
        regions = turn_locations_into_regions(refposns)
        for (ref, start, end) in regions:

            # pull out reference region
            referenceseq = references[ref]
            start = max(start - REGIONSIZE/2, 0)
            end = min(end + REGIONSIZE/2, len(referenceseq))
            regionseq = referenceseq[start:end]

            # align region back to read
            nct = khmer.new_counting_hash(21, 1e5, 4)
            nct.consume(readseq)
            naligner = khmer.ReadAligner(nct, 1, 1.0)
            score, galign = graphAlignment.align_long(nct, naligner, regionseq)

            for n, (a, b) in enumerate(galign):
                if a != '=':
                    break

            o = len(galign)
            while 1:
                (a, b) = galign[o - 1]
                if a != '=':
                    break
                o -= 1

            if '=' in galign[n:o].g:
                assert 0

            gidx = graphAlignment.AlignmentIndex(galign)

            print 'Read %s aligns to %s[%s:%s]' % (read.name,
                                                   ref, start + n, start + o)
            print galign[n:o]
            #print referenceseq[start + n:start + o]

if __name__ == '__main__':
    main()
