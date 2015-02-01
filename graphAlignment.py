import sys
import string
import array
__complementTranslation = string.maketrans('ACTGactg-=Nn', 'TGACtgac-=Nn')


def len_nogap(s):
    return len(s.replace('-', ''))


__complementTranslation = string.maketrans('ACTGactg-=Nn', 'TGACtgac-=Nn')


def reverse_complement(s):
    """
    Build reverse complement of 's', alignment-aware.
    """
    r = array.array('c', s)
    r.reverse()
    r = string.join(r, '')
    c = string.translate(r, __complementTranslation)

    return c


class GraphAlignment(object):
    def __init__(self, g, r):
        assert len(g) == len(r), (g, len(g), r, len(r))
        self.g = g
        self.r = r

        counts = [0] * len(r)
        c = 0
        for ch in r:
            if ch == ch.upper():
                c += 1
            counts.append(c)
        self.truepos = counts

    def reverse_complement(self):
        return GraphAlignment(reverse_complement(self.g),
                              reverse_complement(self.r))

    rc = reverse_complement

    def __add__(self, other):
        return GraphAlignment(self.g + other.g, self.r + other.r)

    def __len__(self):
        return len(self.g)

    def refseqlen(self):
        return self.r.count('G') + self.r.count('C') + self.r.count('T') + \
               self.r.count('A')

    def __getitem__(self, i):
        if isinstance(i, slice):
            start, stop, step = i.indices(len(self.g))
            assert step == 1
            return GraphAlignment(self.g[start:stop], self.r[start:stop])

        return (self.g[i], self.r[i])

    def compare(self):
        n = len(self.g)
        m = 0
        for (a, b) in zip(self.g, self.r):
            if not (a in '-=' or b in '-=' or a != b):
                m += 1

        return m, n

    def __str__(self):
        line1 = []
        line2 = []
        line3 = []
        
        for (a, b) in zip(self.g, self.r):
            line1.append(a)
            line3.append(b)
            
            if a in '-=' or b in '-=' or a != b:
                line2.append(' ')
            else:
                line2.append('|')

        return "%s\n%s\n%s\n" % ("".join(line1), "".join(line2), "".join(line3))


def make_gap(r):
    return GraphAlignment('='*len(r), r)


def _index_alignment(galign, freq=100):
    g_to_r = {}
    r_to_g = {}

    si = 0
    for gi, (a, b) in enumerate(zip(galign.g, galign.r)):
        if gi % freq == 0:
            g_to_r[gi] = si
        if si % freq == 0:
            r_to_g[si] = gi

        if b in 'ACGT':
            si += 1

    return g_to_r, r_to_g


class AlignmentIndex(object):
    def __init__(self, galign, freq=100):
        self.g_to_r, self.r_to_g = _index_alignment(galign, int(freq))
        self.galign = galign
        self.freq = int(freq)

    def get_gi(self, ri):
        rpost = int(ri / self.freq) * self.freq
        gpost = self.r_to_g[rpost]

        diff = ri - rpost
        gi = gpost
        
        while diff > 0:
            (a, b) = self.galign[gi]
            if b in 'ACGT':
                diff -= 1
            gi += 1

        # make sure it's on a valid letter ;)
        while 1:
            (a, b) = self.galign[gi]
            if b in 'ACGT':
                break
            gi += 1
            
        return gi

    def get_ri(self, gi):
        gpost = int(gi / self.freq) * self.freq
        ri = self.g_to_r[gpost]

        diff = gi - gpost
        while diff > 0:
            (a, b) = self.galign[gpost]
            diff -= 1
            gpost += 1
            if b in 'ACGT':
                ri += 1

        return ri


REGIONSIZE = 50


def stitch(galign, K):
    """Stitch together multiple alignments, assuming all overlaps are
    off-by-one on the right end."""
    ga = []
    ra = []

    len_so_far = 0

    n = 0
    for gg in galign[:-1]:
        g = gg.g
        r = gg.r
        assert len(g) == len(r)
        
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

    # did it fail to align across entire segment?
    if truncated:
        aligned_length = galign.refseqlen()
        unaligned_len = len(seq) - aligned_length

        # if we are given next ch, try aligning backwards from next seed
        if next_ch:
            unaligned_seq = seq[-unaligned_len:]
            sr, ggr = align_segment_left(aligner, unaligned_seq + next_ch)
            galign += ggr[:-1]
            score += sr                 # need to adjust score down... @@
        else:
            # if not, just build a gap...
            galign += make_gap(seq[-unaligned_len:])

    return score, galign


def align_segment_left(aligner, seq):
    seq_rc = reverse_complement(seq)
    score, galign = align_segment_right(aligner, seq_rc)
    return score, galign.rc()


def align_long(ct, aligner, sta):
    K = ct.ksize()
    
    # first, pick seeds for each chunk
    seeds = []
    for start in range(0, len(sta), REGIONSIZE):
        region = sta[start:start + REGIONSIZE + K - 1]
        if len(region) < K:
            assert len(sta) - start < REGIONSIZE
            break
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

    # start building piecewise alignments, anchored by seeds.
    alignments = []
    scores = []

    n = 0
    for (start, end) in region_coords[:-1]:
        score, galign = align_segment_right(aligner, sta[start:end],
                                            next_ch=sta[end])
        
        
        scores.append(score)
        alignments.append(galign)
        
        n += 1

    # deal with end (no possibility of alignment from right)
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
