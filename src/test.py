##### This is the file containing much of the code I used for #####
#####   running the keras module in an attempt to model the   #####
#####          colrad function with a neural network          #####

#####    NOTE: Read the comments! Some are quite important    #####

"""Basic imports"""
import numpy as np
from fidanet import *


"""Problematic Imports
   (These can potentially cause issues depending on how everything is installed;
    however, running this file with the following imports commented out before
    runnuing the imports in the console should still work.)"""
import h5py
from keras.models import Sequential, load_model
from keras.layers import Dense, Activation
from keras.optimizers import SGD, Adam
import keras.callbacks


# Just for adding one set of inputs and outputs using inputs taken from
#  uniform distribution. This uses the same original states and dt = 1s for all.
def add_test(ebr = [1e-3,100.], Ter = [1e-3,20.], Tir = [1e-3,20.],
             dener = [1.0e12,2.0e14], Zeffr = [1.,fidanet.impq]):
    eb = np.random.uniform(ebr[0],ebr[1])
    Te = np.random.uniform(Ter[0],Ter[1])
    Ti = np.random.uniform(Tir[0],Tir[1])
    dene = np.random.uniform(dener[0],dener[1])
    Zeff = np.random.uniform(Zeffr[0],Zeffr[1])
    
    denp, denimp = fidanet.getdens(dene,fidanet.impq,Zeff)
    
    fidanet.setplasma(dene,denp,denimp,Te,Ti)
    fidanet.setstates(np.array([1.0e10,0,0,0,0,0]),fidanet.states)
    fidanet.testcol(1,eb,1.0,fidanet.states,fidanet.dens)

    return np.array([eb,Te,Ti,dene,denp,denimp]), fidanet.dens.copy()


# Loops through the single version
def add_tests(num = 1, ebr = [1e-3,100.], Ter = [1e-3,20.], Tir = [1e-3,20.],
             dener = [1.0e12,2.0e14], Zeffr = [1.,fidanet.impq]):
    inputs, outputs  = np.empty((num,6)), np.empty((num,6))
    for i in range(num):
        data = add_test(ebr,Ter,Tir,dener,Zeffr)
        inputs[i] = data[0]
        outputs[i] = data[1]
    return inputs, outputs


# Uses uniform distribution over log space
def add_log_test(ebr = [1e-3,100.], Ter = [1e-3,20.], Tir = [1e-3,20.],
             dener = [1.0e12,2.0e14], Zeffr = [1.,fidanet.impq]):
    eb = np.random.uniform(np.log(ebr[0]),np.log(ebr[1]))
    Te = np.random.uniform(np.log(Ter[0]),np.log(Ter[1]))
    Ti = np.random.uniform(np.log(Tir[0]),np.log(Tir[1]))
    dene = np.random.uniform(np.log(dener[0]),np.log(dener[1]))
    
    eb = np.exp(eb)
    Te = np.exp(Te)
    Ti = np.exp(Ti)
    Zeff = np.random.uniform(Zeffr[0],Zeffr[1])
    
    denp, denimp = fidanet.getdens(dene,fidanet.impq,Zeff)
    
    fidanet.setplasma(dene,denp,denimp,Te,Ti)
    fidanet.setstates(np.array([1.0e10,0,0,0,0,0]),fidanet.states)
    fidanet.testcol(1,eb,1.0,fidanet.states,fidanet.dens)
    
    return np.array([eb,Te,Ti,dene,denp,denimp]), fidanet.dens.copy()


#More looping
def add_log_tests(num = 1, ebr = [1e-3,100.], Ter = [1e-3,20.],
                  Tir = [1e-3,20.], dener = [1.0e12,2.0e14],
                  Zeffr = [1.,fidanet.impq]):
    inputs, outputs  = np.empty((num,6)), np.empty((num,6))
    for i in range(num):
        data = add_log_test(ebr,Ter,Tir,dener,Zeffr)
        inputs[i] = data[0]
        outputs[i] = data[1]
    return inputs, outputs


'''The function below can be commented out due to issues with importing as well.
   This isn't needed if you already have all the data you need in a file.'''

## Generates hdf5 file with input and output data
def save_data(inputs, outputs, file_address):
    f = h5py.File(file_address)
    #groups?
    f.create_dataset('inputs', data = inputs)
    f.create_dataset('outputs', data = outputs)
    #attributes?
    f.close()


## Generates statistics on inputs and outputs
def gen_stats(data, p):
    """
    p is the fraction used for making the model
        p + validation_split = 1
    """
    n =  int(len(data) * p)
    stats = np.array([np.mean(data[:n], axis=0),
                      np.std(data[:n], axis=0)])
    return stats


"""
The following code is commented out as it either can't be run without data,
 or the current values may need changing. What I did was partially guess work.
"""

##### Generating fractional densities
#outputs = np.transpose(np.transpose(outputs) / np.sum(outputs, axis=1))

##### Creates a model
#model = Sequential()
#model.add(Dense(256,input_dim = 6)) #6 for size of inputs
#model.add(Activation('relu'))
#model.add(Dense(128))
#model.add(Activation('tanh'))
#model.add(Dense(64))
#model.add(Activation('tanh'))
#model.add(Dense(6)) # final output must be same size as actual output
#model.compile(optimizer = Adam(), loss = 'mean_absolute_percentage_error',
#              metrics = ['accuracy'])


##Model training
"""
Prior to training, normalize the inputs and outputs.
inputs = (inputs - statin[0]) / statin[1]
outputs = (outputs) / statout[1]
The outputs can be translated first to have a mean of zero, but I avoided that.
"""
#model.fit(inputs,outputs,batch_size=64,epochs=20,
#          callbacks=[keras.callbacks.EarlyStopping(monitor='loss',min_delta=1e-6,patience=5)],
#          validation_split=0.3,verbose=1)



if __name__ == '__main__':
    fidanet.setinputs(fidanet.ai,fidanet.ab,fidanet.impq)
    fidanet.settables()