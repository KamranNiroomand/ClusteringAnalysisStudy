from sklearn.metrics.pairwise import cosine_distances
from sklearn.cluster import SpectralClustering
import numpy as np
import pandas as pd

csv_files = [
    "pca_thickness_pond.csv", "pca_area_pond.csv", "pca_volume_pond.csv", "pca_subcort_pond.csv",
    "pca_thickness_hbn.csv", "pca_area_hbn.csv", "pca_volume_hbn.csv", "pca_subcort_hbn.csv"
]

data_frames = {}
pdf1, pdf2 = [], []

for file_name in csv_files:
    variable_name = file_name.replace(".csv", "")
    data_frames[variable_name] = pd.read_csv(file_name)
    if "pond" in variable_name.lower():
        pdf1.append(data_frames[variable_name])
    elif "hbn" in variable_name.lower():
        pdf2.append(data_frames[variable_name])

distance_df1 = [cosine_distances(f1) for f1 in pdf1]
distance_df2 = [cosine_distances(f2) for f2 in pdf2]

ALPHAs = [4.2, 4.2, 4.2, 4.2]
ALPHAs2 = [4.4, 4.8, 6.4, 4.4]

sim_matx1 = [np.exp(-(distance_df1[i]) ** 2 / (2 * ALPHAs[i] * np.std(distance_df1[i]) ** 2)) for i in range(4)]
sim_matx2 = [np.exp(-(distance_df2[i]) ** 2 / (2 * ALPHAs2[i] * np.std(distance_df2[i]) ** 2)) for i in range(4)]

np.save("XThick_Pond.npy", sim_matx1[0])
np.save("XThick_HBN.npy", sim_matx2[0])
np.save("XArea_Pond.npy", sim_matx1[1])
np.save("XArea_HBN.npy", sim_matx2[1])
np.save("XCortVol_POND.npy", sim_matx1[2])
np.save("XCortVol_HBN.npy", sim_matx2[2])
np.save("XSubVol_POND.npy", sim_matx1[3])
np.save("XSubVol_HBN.npy", sim_matx2[3])

df1, df2 = pd.read_csv("df1.csv"), pd.read_csv("df2.csv")
df1, df2 = df1.reset_index(drop=True), df2.reset_index(drop=True)

Gs1 = [np.load(f) for f in ["XThick_Pond.npy", "XArea_Pond.npy", "XCortVol_POND.npy", "XSubVol_POND.npy"]]
Gs2 = [np.load(f) for f in ["XThick_HBN.npy", "XArea_HBN.npy", "XCortVol_HBN.npy", "XSubVol_HBN.npy"]]

def spec_cluster(matrix, n):
    spectral = SpectralClustering(n_clusters=n, n_neighbors=30, random_state=42, affinity='nearest_neighbors')
    spectral.fit(matrix)
    return spectral.labels_

for i, label in enumerate(["thickness", "area", "volume", "subcortical"]):
    df1[f'cluster_labels_{label}'] = spec_cluster(Gs1[i], 2)
    df2[f'cluster_labels_{label}'] = spec_cluster(Gs2[i], 2)
    df1[f'cluster_labels_{label}_3'] = spec_cluster(Gs1[i], 3)
    df2[f'cluster_labels_{label}_3'] = spec_cluster(Gs2[i], 3)

df1.to_csv("df1_clustered.csv", index=False)
df2.to_csv("df2_clustered.csv", index=False)