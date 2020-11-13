from itertools import product as product
seeds = list(product(['A','T','G','C'],repeat = 2))
seeds = list(map(''.join, seeds))
seq = ''
dic = {}
x = ['A','T','G','C']
y = ['T','A','C','G']
'''
for i in x:
      for j in y:
            seeds.append(i+j)
      y.pop(0)
'''
transdic = seq.maketrans("ATGC","TACG")
rcpairs = [[x,x[::-1].translate(transdic)]for x in seeds]

#print(rcpairs)
def read_file(file_name):
    with open(file_name,'r') as file:
      big = []
      flag = int
      seq = ''
      for line in file:
            line = line.rstrip()
            if (line != '' and line[0] != '>'):
                  flag = 1
                  seq = seq+line
            elif(line == ''):
                  big.append([seq])
                  seq = ''
    return (big)

def gen_frequencies(sequence):
    from itertools import product
    import numpy as np
    import re
    seq=sequence[0]
    lst=[]
    mono_freq = {}
    di_freq = {}
    tri_freq = {}


    dilib =([seq[i:i+2] for i in range(0,(len(seq)-1))])
    #generating dinucleotide product
    dinucleo = list(product(['A','T','G','C'],repeat = 2))
    dinucleo = list(map(''.join, dinucleo))

    #generating counts
    mono_count = np.array((list(map(seq.count,['A','T','G','C']))))
    di_count = np.array((list(map(dilib.count,dinucleo))))

    #generating frequencies

    mono_freq['A']  =  ((mono_count[0]/len(seq))+(mono_count[1]/len(seq)))/2
    mono_freq['T']  =  ((mono_count[0]/len(seq))+(mono_count[1]/len(seq)))/2
    mono_freq['G']  =  ((mono_count[2]/len(seq))+(mono_count[3]/len(seq)))/2
    mono_freq['C']  =  ((mono_count[2]/len(seq))+(mono_count[3]/len(seq)))/2


    for x,y in zip(dinucleo,di_count):
          di_freq[x] =  y/(len(seq)-1)

    return [mono_freq,di_freq]


def gen_biasvalues(frequencies):
      data = []
      mono = frequencies[0]
      di = frequencies[1]
      for pair in rcpairs:
            data.append(((di[pair[0]]+di[pair[1]])/2)/(mono[pair[0][0]]*mono[pair[0][1]]))
      return(data)


def save_to_csv(colnames,rownames,data,filename):
      #creating a dataframe to store data
      data_frame = pandas.DataFrame(columns=colnames,index=rownames)
      for col,seq in zip(data_frame,data):
          data_frame[col]=seq
      #print(data_frame)
      data_frame.to_csv(filename)

file1 = read_file('Zika.fasta')
file2 = read_file('Dengue2.fasta')
file3 = read_file('test.fasta')
test = []
file1_file2_bias = []

for genome in file1:
      f = gen_frequencies(genome)
      file1_file2_bias.append(gen_biasvalues(f))

for genome in file2:
      f = gen_frequencies(genome)
      file1_file2_bias.append(gen_biasvalues(f))
for genome in file3:
      f = gen_frequencies(genome)
      test.append(gen_biasvalues(f))

###Support vector code##########################################################
################################################################################
import pandas as pd
import numpy as np
from sklearn import svm
#data visualization libraries
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(font_scale=1.2)


#specifing data for the model
zika_dengue = file1_file2_bias
type_lable = [0 for i in range(20)]
type_lable.extend([1 for i in range(20)])
print(zika_dengue[1:3])
#print(type_lable)
#fit the svm model

model = svm.SVC(kernel = 'linear')
model.fit(zika_dengue,type_lable)

def Zika_or_Dengue(bias):
      if (model.predict(bias))==0:
            print ("prediction : Zika Virus")
      else:
            print ("prediction : Dengue Virus")


Zika_or_Dengue(test)
