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


def stitch(gs, rs, K):
    ga = []
    ra = []
    for (g, r) in zip(gs[:-1], rs[:-1]):
        ga.append(g[:-K + 1])
        ra.append(r[:-K + 1])

    ga.append(gs[-1])
    ra.append(rs[-1])

    return "".join(ga), "".join(ra)
    

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
        g_alignments.append(g)
        r_alignments.append(r)

        n += 1

    # deal with beginning, too.
    leftend = sta[0:seeds[0] + K]
    leftend_rc = reverse_complement(leftend)
    score, g, r, trunc = aligner.align(leftend_rc)
    
    if trunc:
        score = 0
        g = '-' * len(leftend)
        r = leftend_rc

    g = g[1:]
    r = r[1:]
    g = reverse_complement(g)      # trim off 1-base overlap
    r = reverse_complement(r)      # trim off 1-base overlap

    scores.insert(0, score)
    g_alignments.insert(0, g)
    r_alignments.insert(0, r)

    final_g, final_r = stitch(g_alignments, r_alignments, K)

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
        m = 0
        for n, (a, b, c) in enumerate(zip(g, r, s)):
            line1.append(a)
            line3.append(b)
            if a == '-' or b == '-' or a != b:
                line2.append(' ')
            else:
                line2.append('|')
                m += 1

            if b == '-':
                line4.append('-')
            else:
                line4.append(c)

        ident = int(float(m) / n * 100)

        print '::', record.name, score, len(s), len(line1), "%s%% ident" % ident
        for start in range(0, len(line1), 60):
            print "".join(line1[start:start+60])
            print "".join(line2[start:start+60])
            print "".join(line3[start:start+60])
#            print "".join(line4[start:start+60])
            print ''


if __name__ == '__main__':
    main()
