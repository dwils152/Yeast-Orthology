import numpy as np
import pandas as pd
import sys
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    kmer_freqs = np.load(sys.argv[1])
    #print(kmer_freqs.shape)
    #plot a picture of the np array
    #plt.imshow(kmer_freqs, cmap='hot', interpolation='nearest')
    #save the picture
    #plt.savefig('kmer_freqs.png', dpi=600)

    scaler = StandardScaler()
    X = scaler.fit_transform(kmer_freqs)
    tsne = TSNE(n_components=2)
    tsne_components = tsne.fit_transform(X)
    df = pd.DataFrame({'x': tsne_components[:,0], 'y': tsne_components[:,1]})
    sns.scatterplot(x='x', y='y', data=df)
    plt.savefig('tsne.png', dpi=600)

    
    #np.save('pca_components.npy', pca_components)   

    #for i in range(2, 11):
    #    kmeans_model = KMeans(n_clusters=i, random_state=1).fit(kmer_freqs)
    #    labels = kmeans_model.labels_
    #    ss = metrics.silhouette_score(kmer_freqs, labels, metric='euclidean')
    #    print(ss)

if __name__ == "__main__":
    main()
