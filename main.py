import os
import re
import matplotlib.pyplot as plt
import itertools
import os
from typing import List, Any
import numpy as np
import statistics
import matplotlib.pyplot as plt
import pandas as pd

scriptdir = os.path.dirname(os.path.realpath(__file__))


def drawplot(x, y):
    plt.bar(x, y)
    plt.show()
    plt.gcf().savefig('yourAAs.png')


# def file_lines(directory,x_list): #x_list is the list that is created after getting user input to be counted in files
#     x_list = input("enter AAs with space: ").upper().split()
#     str_count = 0
#     for file in os.listdir(directory + '/files'):  # gets all the files that are in this directory
#         with open(directory + '/files/' + file, 'r') as f:
#             for line in f:
#                 print(i)
#                 str_count += line.count(y)
#                     print(aa_num_y)
#             num_y_list.append(aa_num_y)
#             aa_num_y = 0
# """def compare_files(directory):
#     aa_names_and_numbers = {}
#     file_names = []  # the file names from our folder(directory) will go in this list #later could be used in showing trhe results of each file
#     # used for comparison witrh the file names to be more clear
#     aa_name_x = input(
#         "Please enter the symbol of your amino acid(s):,\n (if more than one, separate them with a space)").upper()
#     print("\n")  # previously referred as name_x
#     x_list = aa_name_x.split()
#     num_y_list = []
#     aa_num_y = 0
#     # aa_list = [] #TODO make it later to put each line as a parameter for furthur usage and manipulation of the file
#     for file in os.listdir(directory + '/files'):  # gets all the files that are in this directory
#
#         file_names.append(file)
#         with open(directory + '/files/' + file, 'r') as f:
#             # lines = f.readlines()
#             # aa_num_y = 0
#             # print(x_list)
#             # below is the original which does not count the next aa
#
#             for i in x_list:
#                 print(i)
#                 for line in f:  # write if to prevent counting the line with ">"
#                     print(line)
#                     aa_num_y += line.count(i)
#                     print(aa_num_y)
#             num_y_list.append(aa_num_y)
#             aa_num_y = 0
#             for line in f:
#                 for i in x_list:
#                     print(i)
#                     for line in f: #write if to prevent counting the line with ">"
#                         print(line)
#                         aa_num_y += line.count(i)
#                         print(aa_num_y)
#                 num_y_list.append(aa_num_y)
#                 aa_num_y = 0
#
#         for i in x_list:  # assiging keys from x_list for the dictionary and the list of each aa symbol number in all files as the value
#             aa_names_and_numbers[
#                 i] = None  # this is for setting dict keys from a list when we don't have values yet, in this case we can also comment it cuz the next line sets keys as well as values
#             aa_names_and_numbers[i] = num_y_list
#
#     print(aa_names_and_numbers)
#     print(file_names)
#     print(num_y_list)
#     # print(aa_lists)
#
#     # for i in x_list:
#     # drawplot(i, int(aa_names_and_numbers[x_list])) #not working
#     # drawplot(x_list, aa_names_and_numbers[i]) #fix plot to show several aa charts in onr
#
#
# # compare_files(scriptdir)
# """


# found here: https://stackoverflow.com/questions/29805642/learning-to-parse-a-fasta-file-with-python/29805905
"""Use a list to accumulate lines until you reach a new id. Then join the lines together and store them
 with the id in a dictionary. The following function takes an open file and yields each pair of (id, sequence)."""


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


# calling the read_fasta() function here, I had to put the whole directory address(as scriptdir variable)
# cuz it is not entered as parameter here
"""for file in os.listdir(scriptdir + '/files'):
    with open(scriptdir + '/files/' + file, 'r') as fp:
        for name, seq in read_fasta(fp):
            print(name, seq)"""


# now we can use this in our compare files function as another function (and used read_fasta method here too):
def name_seq_dict(directory):  # puts names as keys and seq as values of the dictionary
    name_seq_dict = {}
    for file in os.listdir(directory + '/files'):
        with open(directory + '/files/' + file, 'r') as fp:
            for name, seq in read_fasta(fp):
                name_seq_dict[name] = seq
                # print(name, seq)
    # print(name_seq_dict)
    return name_seq_dict
    # return name, seq does not show several files when returning, but surprisingly it prints them all


# print(name_seq_dict(scriptdir))


def compare_files(directory):  # this method asks user their desired aminoacids and it counts them and gives them plots
    aa_names_and_numbers = {}
    file_names = []  # the file names from our folder(directory) will go in this list #later could be used in showing
    # the results of each file
    # used for comparison witrh the file names to be more clear
    x_list = input(
        "Please enter the symbol of your amino acid(s):,\n (if more than one, separate them with a space)") \
        .upper().split()
    print("\n")  # previously referred as name_x
    num_y_list = []
    # aa_num_y = 0
    name_seq = name_seq_dict(directory)
    file_number = len(name_seq)
    a = 0
    b = file_number
    c = 0
    for i in x_list:
        aa_names_and_numbers[i] = None
        aa_num_y = 0
        # num_y_list.clear()
        for seq in name_seq.values():  # iterate through dict values
            aa_num_y = seq.count(i)
            num_y_list.append(aa_num_y)
            if len(num_y_list) % file_number == 0:
                aa_names_and_numbers[i] = num_y_list[a:b]
                a += file_number  # adding file number makes sublists in range of number of each AA symbol in all the
                # files, it separates them from the num_y_list because if we give R and L as our aa symbols,
                # our num_y_list will be like [1,2,3,4,5,6] which the first 3 numbers are R nums in our 3 files(in
                # this case) and second 3 numbers are L nums in our 3 files
                b += file_number
            else:
                continue

    print(aa_names_and_numbers)
    y_list = []
    plot_num = file_number  # number of output plots, changed name just to be more clear that we want the same number
    # plots as the files we have
    m = 0
    n = len(x_list)
    while plot_num != 0:
        # y_list = []

        for i in x_list:
            # print(aa_names_and_numbers[i][c])
            # y_list = [] # we are going to take corresponding y values from dict list
            # values, the 0 index is the R of the first file,the 1st index is
            # the R num of second file ,this goes the same for the second dict key and value
            y_list.append([aa_names_and_numbers[i][c]])
            # print(y_list)
        y_list_merged = list(itertools.chain.from_iterable(
            y_list))  # to remove the brackets and make one flat list out of list of lists , found here: https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
        drawplot(x_list, y_list_merged[m:n])  # fix plot to show several aa charts in onr #does not work
        c += 1
        plot_num -= 1
        m += len(x_list)
        n += len(x_list)
    return aa_names_and_numbers

    # TODO: it does not separate each sequence's X AA counts, e.g: it sums all Rs from all seqs together=> put the
    #  following
    """# aa_names_and_numbers dict code in the upper loop
    for i in x_list:  # assigning keys from x_list for the dictionary and the list of each aa symbol number in all files
        # as the value
        aa_names_and_numbers[
            i] = None  # this is for setting dict keys from a list when we don't have values yet, in this case we can
        # also comment it cuz the next line sets keys as well as values
        aa_names_and_numbers[i] = num_y_list[x_list.index(
            i)]  # fixed now, don't know for several files yet# fix this, it shows the same and whole num_y_list in front of all aa
        # symbols #!get the same index as [i] in x_list for the y list as well"""


# compare_files(scriptdir)

def math_operation(aa_and_num_dict):
    for key in aa_and_num_dict:
        print("Max amount of ", key, "is: ", max(aa_and_num_dict[key]))
        print("Min amount of ", key, "is: ", min(aa_and_num_dict[key]))
        print("Average of ", key, "is: ", statistics.mean(aa_and_num_dict[key]))
    return max(aa_and_num_dict[key]), min(aa_and_num_dict[key]), statistics.mean(aa_and_num_dict[key])
#could I actually return all three of them? #TODO: check later if this is true




#math_operation(compare_files(scriptdir))
