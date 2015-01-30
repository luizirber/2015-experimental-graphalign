#! /usr/bin/env python
import sys
import argparse
import screed
import khmer
from graphAlignment import GraphAlignment, make_gap, reverse_complement

REGIONSIZE = 50

def len_nogap(s):
    return len(s.replace('-', ''))


def stitch(galign, K):
    ga = []
    ra = []

    n = 0
    for gg in galign[:-1]:
        g = gg.g
        r = gg.r
        
        #assert len(g) == len(r), n
        minlen = min(len(g), len(r))
        g = g[:minlen] # @CTB
        r = r[:minlen]
        
        ga.append(g[:-K + 1])
        ra.append(r[:-K + 1])
        n += 1

    ga.append(galign[-1].g)
    ra.append(galign[-1].r)

    return GraphAlignment("".join(ga), "".join(ra))


def align_segment_right(aligner, seq, next_ch=None):
    assert len(seq) >= 21
    score, g, r, truncated = aligner.align(seq)
    galign = GraphAlignment(g, r)

    if truncated:
        aligned_length = galign.refseqlen()
        unaligned_len = len(seq) - aligned_length

        if next_ch:                     # try aligning backwards from next seed
            unaligned_seq = seq[-unaligned_len:]
            sr, ggr = align_segment_left(aligner, unaligned_seq + next_ch)

            middle_len = len(seq) - galign.refseqlen() - ggr.refseqlen()
            galign += make_gap(unaligned_seq[:middle_len])[:-1]
            score += sr
        else:
            galign += make_gap(seq[-unaligned_len:])

    return score, galign


def align_segment_left(aligner, seq):
    assert len(seq) > 21
    seq_rc = reverse_complement(seq)
    score, g, r, truncated = aligner.align(seq_rc)
    galign = GraphAlignment(g, r)

    if truncated:
        aligned_length = min(len_nogap(g), len_nogap(r))
        unaligned_len = len(seq) - aligned_length
        galign += make_gap(seq_rc[-unaligned_len:])

    return score, galign.reverse_complement()


def align_long(ct, aligner, sta):
    K = ct.ksize()
    
    # first, pick seeds for each chunk
    seeds = []
    for start in range(0, len(sta), REGIONSIZE):
        region = sta[start:start + REGIONSIZE + K - 1]
        seed_pos = find_highest_abund_kmer(ct, region)
        seeds.append(start + seed_pos)

    print >>sys.stderr, 'picked %d seeds' % len(seeds)
    assert len(seeds) == len(set(seeds))

    # then, break between-seed intervals down into regions, starting at
    # first seed.
    region_coords = []
    for i in range(len(seeds) - 1):
        seed_pos = seeds[i]
        end_seed = seeds[i + 1] - 1
        region_coords.append((seed_pos, end_seed + K))

    if len(sta) - seeds[-1] > K:
        region_coords.append((seeds[-1], len(sta)))

    # start building piecewise alignments
    alignments = []
    scores = []

    n = 0
    for (start, end) in region_coords[:-1]:
        #print 'aligning', n, start, end, (end - start - K)
        score, galign = align_segment_right(aligner, sta[start:end],
                                          next_ch=sta[end])
        
        scores.append(score)
        alignments.append(galign)
        
        n += 1

    # deal with end:
    (start, end) = region_coords[-1]
    score, galign = align_segment_right(aligner, sta[start:end])
    
    alignments.append(galign)
    scores.append(score)
    
    # deal with beginning, too: reverse align from first seed.
    leftend = sta[0:seeds[0] + K]
    score, galign = align_segment_left(aligner, leftend)

    galign = galign[:-1]                # trim off seed k-mer
    alignments.insert(0, galign)
    scores.insert(0, score)

    # stitch all the alignments together
    final = stitch(alignments, K)

    return sum(scores), final


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

        score, alignment = align_long(ct, aligner, s)

        g = alignment.g
        r = alignment.r

        m, n = alignment.compare()
        print record.name, m, n, n - m, "%.3f%%" % (float(m)/ n * 100)
        for start in range(0, len(alignment), 60):
            print alignment[start:start+60]

        print len(s), alignment.refseqlen()


if __name__ == '__main__':
    main()
