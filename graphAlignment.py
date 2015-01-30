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
               self.r.count('A') + self.r.count('N')

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
