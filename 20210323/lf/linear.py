from Bio.Align.substitution_matrices import Array
import numpy
from Bio import SeqIO
from time import time
from utils import get_random_string, load_score, save_result
import sys


def alignment(str1, str2, gap, cost, tags):
    n = len(str1)
    m = len(str2)

    Table = [[0 for i in range(len(str1)+1)] for j in range(len(str2)+1)]

    for i in range(1, m+1):
        Table[i][0] = i*gap
    for j in range(1, n+1):
        Table[0][j] = j*gap
    for i in range(1, m+1):
        for j in range(1, n+1):
            Table[i][j] = min(
                Table[i][j - 1] + gap,
                Table[i - 1][j] + gap,
                Table[i - 1][j - 1] +
                cost[tags[str2[i - 1]], tags[str1[j - 1]]]
            )

    return Table


def RecurBackTrackWrapper(B, A, T, cost, tags, g):
    aligs = []
    current_alig = []

    def RecurBackTrack(i, j):
        if (i > 0) and (j > 0) and T[i][j] == T[i-1][j-1] + cost[tags[A[i-1]], tags[B[j-1]]]:
            current_alig.append((A[i-1], B[j-1]))
            RecurBackTrack(i-1, j-1)
            current_alig.pop()
        if (i > 0) and (j >= 0) and T[i][j] == T[i-1][j] + g:
            current_alig.append((A[i-1], "-"))
            RecurBackTrack(i-1, j)
            current_alig.pop()
        if (i >= 0) and (j > 0) and T[i][j] == T[i][j-1] + g:
            current_alig.append(("-", B[j-1]))
            RecurBackTrack(i, j-1)
            current_alig.pop()
        if (i == 0) and (j == 0):
            alig = current_alig.copy()
            alig.reverse()
            aligs.append(numpy.transpose(alig))

    RecurBackTrack(len(A), len(B))
    return aligs


def IterativeBackTrack(B, A, T, cost, tags, g):
    i = len(A)
    j = len(B)
    alig = []
    while (i > 0 and j > 0):
        if (i > 0) and (j > 0) and T[i][j] == T[i-1][j-1] + cost[tags[A[i-1]], tags[B[j-1]]]:
            alig.append((A[i-1], B[j-1]))
            i = i-1
            j = j-1
        elif (i > 0) and (j >= 0) and T[i][j] == T[i-1][j] + g:
            alig.append((A[i-1], "-"))
            i = i-1
        elif (i >= 0) and (j > 0) and T[i][j] == T[i][j-1] + g:
            alig.append(("-", B[j-1]))
            j = j-1
    return numpy.array(alig).transpose()


def main(str1, str2, gap, cost, tags, save):
    str1 = str1.upper()
    str2 = str2.upper()

    begin = time()
    T = alignment(str1, str2, gap, cost, tags)
    after_align = time()
    align = IterativeBackTrack(str1, str2, T, cost, tags, gap)
    after_backtrack = time()

    optimal = [
        "".join(seq[::-1]) for seq in align
    ]

    if save:
        save_result(optimal)

    print("\nInput:", str1, str2, gap, sep="\n\n")
    print("\nSubstitution matrix:")
    print(cost)
    print("\nAlignment graph:")
    print(numpy.matrix(T))
    print("\nOptimal alignment:")
    print("\n".join(optimal))
    # print("\nNumber of alignments:")
    # print(len(alignments))
    print("\nOptimal cost:")
    print(T[-1][-1])
    print("\nAlignment time:")
    print(after_align - begin)
    print("\nBacktrack time:")
    print(after_backtrack - after_align)


def correctness(gap, cost, tags):
    seqs = list(SeqIO.parse("correctness.fasta", "fasta"))
    scores = numpy.zeros((len(seqs), len(seqs)))
    for i, seqi in enumerate(seqs):
        for j, seqj in enumerate(seqs):
            T = alignment(seqi.upper(), seqj.upper(), gap, cost, tags)
            scores[i, j] = T[-1][-1]
    print(f"{gap}*k", sep="\n")
    print(scores)


def eval_running_time(gap, cost, tags):
    import matplotlib.pyplot as pyplot

    xs = []
    ys = []
    for i in range(5):
        N = 256 * 2 ** i
        str1 = get_random_string(N)
        str2 = get_random_string(N)

        begin = time()
        T = alignment(str1, str2, gap, cost, tags)
        alig = IterativeBackTrack(str1, str2, T, cost, tags, gap)
        elapsed = time() - begin

        xs.append(N)
        ys.append(elapsed / N ** 2)
        print(N, elapsed)

    pyplot.scatter(xs, ys)
    pyplot.plot(xs, ys)
    pyplot.semilogx(base=2)
    pyplot.xlabel('Sequence length')
    pyplot.ylabel('Running time divided by N^2')
    pyplot.title('Running time complexity evaluation')
    pyplot.show()


if __name__ == '__main__':

    cost, tags = load_score()

    gap = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    save = len(sys.argv) > 2

    if (gap == 0):
        raise Exception("Wrong gap cost")

    str1 = list(SeqIO.parse("seq1.fasta", "fasta"))[0]
    str2 = list(SeqIO.parse("seq2.fasta", "fasta"))[0]

    # str1 = 'acgtgtcaacgt'
    # str2 = 'acgtcgtagcta'

    # str1 = 'AATAAT'
    # str2 = 'AAGG'

    # str1 = 'tccagaga'
    # str2 = 'tcgat'

    # str1 = get_random_string(4000)
    # str2 = get_random_string(4000)

    main(str1, str2, gap, cost, tags, save)
    # correctness(gap, cost, tags)
    # eval_running_time(gap, cost, tags)
