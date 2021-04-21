import random
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def get_random_string(length):
    letters = ["A", "T", "G", "C"]
    return ''.join(random.choice(letters) for i in range(length))


def load_score():
    with open("score.txt", "r") as f:
        lines = [line.split() for line in f.readlines()]
        tags = [line[0] for line in lines]
        tags = {name: i for i, name in enumerate(tags)}
        if (set(tags) != {"A", "C", "G", "T"}):
            raise Exception("Wrong letters in the alphabet")

        cost = [line[1:] for line in lines]
        cost = np.matrix(cost, int, )
        if (cost.shape != (4, 4)):
            raise Exception("Wrong shape of the score matrix")

    return cost, tags


def save_result(optimal):
    frags = []
    for i, x in enumerate(optimal):
        record = SeqRecord(Seq(x), f"seq{i+1}", "", "")
        frags.append(record)
    SeqIO.write(frags, "output.fasta", "fasta")


if __name__ == '__main__':
    print(get_random_string(50))
