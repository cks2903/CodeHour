import pandas as pd

start = pd.read_pickle('start.pickle')
end = pd.read_pickle('end.pickle')

old = pd.read_csv('old.gff3', sep='\t', header=None, comment='#')
old.columns = ['Col_' + str(x) for x in range(9)]
# old.Col_8 = old.Col_8.apply(lambda x: x.replace(' ', ''))

newcol3 = []
for i in old.index:
    if i in start.keys():
        cond = old.Col_3.iloc[i] != start[i]
        if cond.bool():
            newcol3.append(int(start[i]))
        else:
            newcol3.append(old.Col_3.iloc[i])
    else:
        newcol3.append(old.Col_3.iloc[i])

newcol4 = []
for i in old.index:
    if i in end.keys():
        cond = old.Col_4.iloc[i] != end[i]
        if cond.bool():
            newcol4.append(int(end[i]))
        else:
            newcol4.append(old.Col_4.iloc[i])
    else:
        newcol4.append(old.Col_4.iloc[i])

old.iloc[:, 3] = newcol3
old.iloc[:, 4] = newcol4

old.to_csv('annotation_new_version.gff3', sep='\t', index=None, header=None)
