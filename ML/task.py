import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt


colors = ["red", "green", "blue", "black", "yellow"]
for i in range(6):
    df = pd.read_csv(f"datasets\{i + 1}.csv", index_col=None)
    
    X = df[["x", "y"]].to_numpy()
    y = df["class"].to_numpy()
    n_clusters = np.unique(y).shape[0]
    kmeans = KMeans(n_clusters=n_clusters, random_state=42).fit(X)
    y_pred = kmeans.predict(X)
    
    plt.figure(figsize=(8, 8))
    for j in range(n_clusters):
        color = colors[j]
        res = X[y_pred == j, :]
        plt.scatter(res[:, 0], res[:, 1], c=color)
