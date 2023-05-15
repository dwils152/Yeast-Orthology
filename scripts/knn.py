import numpy as np
import pandas as pd
import sys
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from gap_statistic import OptimalK

def main():

    kmer_freqs = np.load(sys.argv[1])

    scaler = StandardScaler()
    X = scaler.fit_transform(kmer_freqs)
    tsne = TSNE(n_components=2)
    tsne_components = tsne.fit_transform(X)
    df = pd.DataFrame({'x': tsne_components[:,0], 'y': tsne_components[:,1]})
    sns.scatterplot(x='x', y='y', data=df, color='#3C5488FF', size=0.05, linewidth=0, legend=False, alpha=0.25)
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.savefig('tsne.png', dpi=600)
    
    optimalK = OptimalK(n_jobs=1)
    n_clusters = optimalK(tsne_components, cluster_array=np.arange(1, 15))

    optimalK.gap_df.to_csv('gap.csv', index=False)

    #for k in range(1, 11):
    #    kmeans_model = KMeans(n_clusters=k, random_state=1).fit(tsne_components)
    #    labels = kmeans_model.labels_
    #    ss = metrics.silhouette_score(kmer_freqs, labels, metric='euclidean')
    #    print('k =', k, 'Silhouette Score:', ss)

if __name__ == "__main__":
    main()
