

import numpy as np


def cal_frequency(sequence,knumber=1):
      '''cal_frequency(sequence,knumber= default 1){this function returns monopeptide/dipeptide frequency array(numpy)
	for a given sequence}'''
      counts = []
      if knumber == 1:
            Amino = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
            amino_dic = dict.fromkeys(Amino,0)

            for y in range(0,len(sequence)):
                  amino_dic[sequence[y]] = amino_dic.get(sequence[y],0) + (1/len(sequence))
            counts = list(amino_dic.values())

      if knumber == 2:
            Amino = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
            Diamino= [i+j for i in Amino for j in Amino]
            diamino_dic = dict.fromkeys(Diamino,0)


            for y in range(0,len(sequence)):
                  if len(sequence[y:y+2])==2:
                        diamino_dic[sequence[y:y+2]] = diamino_dic.get(sequence[y:y+2],0) + (1/(len(sequence)-1))
            
            counts = list(diamino_dic.values())
            
      counts = np.array(counts,dtype='float32')


      return counts


      
def extract_seq(file):
    ''' this function can be used to extract sequences from a multiple fasta file, sequences are returned in an array'''
    sequences = []
    seq = str()
    with open(file,'r') as file:
        for line in file:
            line = line.rstrip()
            if(line[0]!='>'):
                
                seq = seq+line
            else:
                if seq != '':
                      sequences.append(seq)
    return sequences

def distance(a,b):
      '''This fucntion calculates euclidean distance between 2 numpy arrays, retruns a numpy array'''
      squ_diff = (a-b)**2
      distance = (np.sum(squ_diff,axis=1))**0.5
      return distance
    
def feature_calculation(file,knumber=1):
      '''feature_calculation(file,knumber= default 1) {calculates sequence features from a multifasta sequence file} '''
      features = np.ones([1,20],dtype='float32')
      a =(extract_seq(file))
      
      for seq in a:
            freq=cal_frequency(seq,knumber)
            features=np.append(features,[freq],0)

      return features
    
def optimize_k(train,train_classes,test,test_class,max_k = 9):
	'''optimize_k(train,train_classes,test,test_class,max_k = 9)'''
      predictions = []
      test_class = test_class[1]
      for test_seq in test:
            array2 = []
            
            distance_array = distance(train,test_seq)
            
            sorted_array = np.msort(distance_array)
            #print(sorted_array)
            for k in range(2,max_k,2):
                  sorted_k = sorted_array[k]
                  
                  neighbours = np.ndarray.flatten(np.array(np.where(distance_array <= sorted_k),'int64'))

                  neighbour_classes = np.take(train_classes,neighbours)
                  
                  count = np.count_nonzero(neighbour_classes == test_class)
                  
                  accuracy = (count/k)
                  
                  if (accuracy > 0.5):
                        array2.append(1)
                  else:
                        array2.append(0)

            predictions.append(array2)
      return predictions

def accuracy_k(predictions):
	'''this function calculates the accuracy of the predictions'''
      predict = (predictions)
      k3,k5,k7,k9 = 0,0,0,0
      for a,b,c,d in  predictions:
            k3 += a
            k5 += b
            k7 += c
            k9 += d
      percent_accuracy=(np.array([k3,k5,k7,k9],dtype='float32')/24)*100
      print (percent_accuracy)
      print ("  k3   k5   k7   k9")          

#Data extraction and feature calculation
train_positive = feature_calculation('positive_train_seq.txt')
train_negative = feature_calculation('negative_train_seq.txt')
train_positive = np.delete(train_positive,0,0)
train_negative = np.delete(train_negative,0,0)
train = np.append(train_positive,train_negative,0)

#training classes
ones = np.ones((1,74),dtype='int32')
zeros = np.zeros((1,74))
train_classes = np.append(ones,zeros)

#testing data extraction and feature calculation
test_positive = feature_calculation('positive_test_seq.txt')
test_negative = feature_calculation('negative_test_seq.txt')
test_positive = np.delete(test_positive,0,0)
test_negative = np.delete(test_negative,0,0)

#testing classes
test_negative_classes = np.zeros((len(test_negative)))
test_positive_classes = np.ones((len(test_positive)))


optimize = (optimize_k(train,train_classes,test_positive,test_positive_classes))
accuracy_k(optimize)
print(len(test_positive))
'''
print(np.random.rand(3))
print(np.random.randn(6))
x = (np.random.randint(1,148,148,'int32')[:20])
y = (np.random.randint(1,148,148,'int32'))
print(np.take(y,x))'''

'''this function splits the data setinto test and trian sets'''
'''
#def train_test_split(feature_vector, class_vector):
shuffled = ((np.random.randint(1,148,148,'int32')))
train = shuffled[0:30]
test  = shuffled[30:]
#np.take(array,indices)
print(shuffled)
print(train)'''

