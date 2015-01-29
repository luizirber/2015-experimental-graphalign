#! /usr/bin/env python
import sys
import argparse
import screed
import khmer
import string
import array

REGIONSIZE = 50

__complementTranslation = string.maketrans('ACTGactg-Nn', 'TGACtgac-Nn')


def reverse_complement(s):
    """
    Build reverse complement of 's', alignment-aware.
    """
    r = array.array('c', s)
    r.reverse()
    r = string.join(r, '')
    c = string.translate(r, __complementTranslation)

    return c


def align_long(ct, aligner, sta):
    K = ct.ksize()
    
    # first, pick seeds for each chunk
    seeds = []
    for start in range(0, len(sta), REGIONSIZE):
        region = sta[start:start + REGIONSIZE + K]
        seed_pos = find_highest_abund_kmer(ct, region)
        seeds.append(start + seed_pos)

    print >>sys.stderr, 'picked %d seeds' % len(seeds)

    # then, break between-seed intervals down into regions, starting at
    # first seed.
    region_coords = []
    for i in range(len(seeds) - 1):
        seed_pos = seeds[i]
        end_seed = seeds[i + 1] - 1
        region_coords.append((seed_pos, end_seed + K))

    region_coords.append((seeds[-1], len(sta)))

    # start building piecewise alignments
    g_alignments = []
    r_alignments = []
    scores = []

    n = 0
    for (start, end) in region_coords:
        #print 'aligning', n, start, end, (end - start - K)
        score, g, r, trunc = aligner.align(sta[start:end])

        if trunc:
            g = '-' * (end-start)
            r = sta[start:end]
            score = 0

        scores.append(score)
        if end != len(sta):
            g_alignments.append(g[:-K + 1])
            r_alignments.append(r[:-K + 1])
        else:
            g_alignments.append(g)
            r_alignments.append(r)

        n += 1

    # deal with beginning, too.
    leftend = sta[0:seeds[0] + K + 1]
    leftend_rc = reverse_complement(leftend)
    score, g, r, trunc = aligner.align(leftend_rc)
    
    if trunc:
        score = 0
        g = '-' * len(leftend)
        r = leftend_rc

    g = reverse_complement(g[K + 1:])
    r = reverse_complement(r[K + 1:])

    scores.insert(0, score)
    g_alignments.insert(0, g)
    r_alignments.insert(0, r)

    final_g = "".join(g_alignments)
    final_r = "".join(r_alignments)

    return sum(scores), final_g, final_r


def find_highest_abund_kmer(ct, r):

    pos = 0
    abund = ct.get(r[:ct.ksize()])
    for kstart in range(1, len(r) - ct.ksize() + 1):
        kmer = r[kstart:kstart + ct.ksize()]
        count = ct.get(kmer)

        if count > abund:
            pos = kstart

    return pos


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('table')
    parser.add_argument('ref')
    args = parser.parse_args()

    ct = khmer.load_counting_hash(args.table)
    aligner = khmer.new_readaligner(ct, 5, 1.0)

    for record in screed.open(args.ref):
        s = record.sequence
        s = s.replace('N', 'A')

        score, graph_alignment, read_alignment = align_long(ct, aligner, s)

        g = graph_alignment
        r = read_alignment

        line1 = []
        line2 = []
        line3 = []
        line4 = []
        for n, (a, b, c) in enumerate(zip(g, r, s)):
            line1.append(a)
            line3.append(b)
            if a != b:
                line2.append(' ')
            else:
                line2.append('|')
            line4.append(c)

        print '::', record.name, score, len(s), len(line1)
        for start in range(0, len(line1), 60):
            print "".join(line1[start:start+60])
            print "".join(line2[start:start+60])
            print "".join(line3[start:start+60])
            print "".join(line4[start:start+60])
            print ''


if __name__ == '__main__':
    main()
