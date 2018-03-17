from time import time
import sklearn.ensemble as ske
from clustering_cdhit import readJson
import numpy as np
import sys, os
import pickle
        


def train_model(numpy_matrix, region_type, samples_file, samples_sp2): 
        
        samples_dic = readJson(samples_file)
        samplesnumber = len(samples_dic.keys())
        gene_matrix=np.load(numpy_matrix)

        Y = np.zeros((1,samplesnumber))
        Y[:,0:samples_sp2]=1
        print (Y)
        Y=Y.T
        print (Y.shape)

        t0 = time()
        model = ske.RandomForestRegressor(bootstrap=True, n_estimators=1000, n_jobs= 16, max_depth=None, min_samples_split=1.0, random_state=0)
        model.fit(gene_matrix, Y.ravel())
        print("analysis done in in %0.3fs." % (time() - t0)) 

        try:
            os.stat("models")
        except:
            os.mkdir("models")

        # save the model to disk
        filename = 'models/finalized_model_%s.sav' %region_type
        pickle.dump(model, open(filename, 'wb'))



if __name__ == "__main__":

        numpy_matrix = sys.argv[1]
        region_type = sys.argv[2]
        samples_file = sys.argv[3]
        samples_sp2 = sys.argv[4]

        train_model(numpy_matrix, region_type, samples_file, int(samples_sp2))
        #train_model('../matrix/0.7_PERC_gene.npy', 'gene', '../dics/samples.json', 2)