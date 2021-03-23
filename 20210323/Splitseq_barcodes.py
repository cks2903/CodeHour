
import itertools

barcodes1 = open("Barcodes1.txt",  "r")
barcodes2 = open("Barcodes2.txt",  "r")
barcodes3 = open("Barcodes3.txt", "r")


# Make sure the are no empty line in the input files
list1 = barcodes1.read().split("\n")
list2 = barcodes2.read().split("\n")
list3 = barcodes3.read().split("\n")

# initializing list of list
all_list = [list1, list2, list3]

# printing lists
#print("The original lists are : " + str(all_list))

res1 = itertools.product(*all_list)

#print("All possible permutations are : " + str(res1))


file1 = open("Splitseq_barcodes.txt", "w")
for x in res1:
    barcodes = (' '.join(x))
    barcodes1 = str(barcodes.split("\n"))
    # print(barcodes1)
    barcodes2 = barcodes1.replace(" ", "")[2:26]
    # print(barcodes2)
    file1.write("%s\n" % str(barcodes2))
file1.close()
