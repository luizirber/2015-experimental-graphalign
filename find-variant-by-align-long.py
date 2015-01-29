#! /usr/bin/env python
import sys
import argparse
import screed
import khmer
import string
import array

REGIONSIZE = 50

__complementTranslation = string.maketrans('ACTGactg-Nn', 'TGACtgac-Nn')


def len_nogap(s):
    return len(s.replace('-', ''))


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

    n = 0
    for (g, r) in zip(gs[:-1], rs[:-1]):
        #assert len(g) == len(r), n
        minlen = min(len(g), len(r))
        g = g[:minlen] # @CTB
        r = r[:minlen]
        
        ga.append(g[:-K + 1])
        ra.append(r[:-K + 1])
        n += 1

    ga.append(gs[-1])
    ra.append(rs[-1])

    return "".join(ga), "".join(ra)


def align_segment_right(aligner, seq, next_ch=None):
    assert len(seq) >= 21
    score, g, r, truncated = aligner.align(seq)

    if truncated:
        aligned_length = min(len_nogap(g), len_nogap(r))
        unaligned_len = len(seq) - aligned_length

        if next_ch:                     # try aligning backwards from next seed
            unaligned_seq = seq[-unaligned_len:]
            sr, gr, rr, truncr = align_segment_left(aligner,
                                                    unaligned_seq + next_ch)
            #assert truncr               # _should_ be truncated!

            middle_len = min(len_nogap(g + gr), len_nogap(r + rr))
            g += '-' * middle_len + gr[:-1]
            r += unaligned_seq[:middle_len] + rr[:-1]
            score += sr
        else:
            g += '-' * unaligned_len
            r += seq[-unaligned_len:]

    return score, g, r, truncated


def align_segment_left(aligner, seq):
    assert len(seq) > 21
    seq_rc = reverse_complement(seq)
    score, g, r, truncated = aligner.align(seq_rc)

    if truncated:
        aligned_length = min(len_nogap(g), len_nogap(r))
        unaligned_len = len(seq) - aligned_length

        g += '-' * unaligned_len
        r += seq_rc[-unaligned_len:]

    g = reverse_complement(g)
    r = reverse_complement(r)
    
    return score, g, r, truncated


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
    g_alignments = []
    r_alignments = []
    scores = []

    n = 0
    for (start, end) in region_coords[:-1]:
        #print 'aligning', n, start, end, (end - start - K)
        score, g, r, trunc = align_segment_right(aligner, sta[start:end],
                                                 next_ch=sta[end])

        scores.append(score)
        g_alignments.append(g)
        r_alignments.append(r)

        n += 1

    # deal with end:
    (start, end) = region_coords[-1]
    score, g, r, trunc = align_segment_right(aligner, sta[start:end])
    
    scores.append(score)
    g_alignments.append(g)
    r_alignments.append(r)

    # deal with beginning, too: reverse align from first seed.
    leftend = sta[0:seeds[0] + K]
    score, g, r, trunc = align_segment_left(aligner, leftend)
    
    g = g[:-1]                          # trim off seed match
    r = r[:-1]
    
    scores.insert(0, score)
    g_alignments.insert(0, g)
    r_alignments.insert(0, r)

    # stitch all the alignments together
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
        for n, (a, b) in enumerate(zip(g, r)):
            line1.append(a)
            line3.append(b)
            if a == '-' or b == '-' or a != b:
                line2.append(' ')
            else:
                line2.append('|')
                m += 1

        ident = float(m) / n * 100

        print '::', record.name, score, len(s), len(line1), "%.2f%% ident" % ident, m, n
        for start in range(0, len(line1), 60):
            print "".join(line1[start:start+60])
            print "".join(line2[start:start+60])
            print "".join(line3[start:start+60])
#            print "".join(line4[start:start+60])
            print ''


if __name__ == '__main__':
    main()
