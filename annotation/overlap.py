import pandas as pd
import pickle
from joblib import Parallel, delayed
import time


def timer(func):
    def wrapper(*args, **kwargs):
        t_start = time.time()
        result = func(*args, **kwargs)
        t_total = time.time() - t_start
        print('{} took {}s'.format(func.__name__, t_total))
        return result
    return wrapper


old = pd.read_csv('old.gff3', sep='\t', header=None, comment='#')
kai = pd.read_csv('kai.gff3', sep='\t', header=None, comment='#')

old.columns = ['Col_' + str(x) for x in range(9)]
kai.columns = old.columns

old8 = old.Col_8.apply(lambda x: x.split(';')[0])
kai8 = kai.Col_8.apply(lambda x: x.split(';')[0])


@timer
def position(a):
    start = dict()
    end = dict()
    for i in a:
        d = old8.loc[i] == kai8
        d = d.index[d]
        if len(d) == 1:
            x = old.Col_3.loc[i] != kai.Col_3.loc[d]
            y = old.Col_4.loc[i] != kai.Col_4.loc[d]
            if x.bool():
                start[i] = kai.Col_3.loc[d]
            if y.bool():
                end[i] = kai.Col_4.loc[d]
    return start, end


result = Parallel(n_jobs=32)(delayed(position)(i) for i in [old.index])

start_ = result[0][0]
end_ = result[0][1]

pickle.dump(start_, open('start.pickle', 'wb'))
pickle.dump(end_, open('end.pickle', 'wb'))
