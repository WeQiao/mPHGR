import pandas as pd
import numpy as np
import os

os.chdir('C:\mPHGR\PTB10013')

file_path2=f'control.txt'
datacon=pd.read_csv(file_path2, sep='\t',  header=None)
corr_control=datacon.T.corr()
corr_control.to_numpy()


def calculate_and_save_matrices(file_path,index):
    data = pd.read_csv(file_path, sep='\t',  header=None)
    correlation_matrix = data.T.corr()
    correlation_matrix_ab = np.absolute(correlation_matrix.to_numpy()-corr_control.to_numpy())#差异皮尔逊相关矩阵

    correlation_matrix_H = np.zeros_like(correlation_matrix_ab)
    removed_rows = []

    for i in range(correlation_matrix_ab.shape[0]):

        indices = np.where(correlation_matrix_ab[i, :] >= 0.5)[0]

        if indices.size == 1:

            removed_rows.append(i)
            continue

        if indices.size > 1:

            for j in indices:
                indices3 = indices[indices != j]
                avg_value = np.mean(correlation_matrix_ab[indices3, j])
                correlation_matrix_H[i, j] = avg_value

        else:

            correlation_matrix_H[i, :] = 0

    if removed_rows:
        correlation_matrix_H = np.delete(correlation_matrix_H, removed_rows, axis=0)


    W = np.zeros((correlation_matrix_H.shape[1], correlation_matrix_H.shape[0]))
    for i in range(correlation_matrix_H.shape[1]):
        for j in range(correlation_matrix_H.shape[0]):
            if correlation_matrix_H[j, i] > 0:
                W[i, j] = np.sum(correlation_matrix_H[j, :])

    DV = np.diag(np.sum(W, axis=1))
    DE = np.diag(np.sum(correlation_matrix_H, axis=1))

    epsilon = 1e-10
    DV_reg = DV + epsilon * np.eye(correlation_matrix_H.shape[1])
    DE_reg = DE + epsilon * np.eye(correlation_matrix_H.shape[0])

    P = np.linalg.inv(DV_reg) @ W @ np.linalg.inv(DE_reg) @ correlation_matrix_H

    P_csv_file = f'P_{index}.csv'
    pd.DataFrame(P).to_csv(P_csv_file, index=False, header=False)

k = 17 #The number of single samples (time points) other than the control group
for i in range(1, k + 1):
    file_name = f'test2_data{i}.txt'
    calculate_and_save_matrices(file_name,i)


def read_and_variance(file_name):
    data = np.loadtxt(file_name)
    variances = np.var(data, axis=1)
    variances = variances.reshape(-1, 1)
    return variances


def calculate_d(file_path):
    data = pd.read_csv(file_path, sep='\t',  header=None)
    correlation_matrix = data.T.corr()

    correlation_matrix_ab = np.absolute(correlation_matrix.to_numpy())
    d = np.ones((correlation_matrix_ab.shape[0], 1))

    correlation_matrix_H = np.zeros_like(correlation_matrix_ab)
    removed_rows = []

    for i in range(correlation_matrix_ab.shape[0]):

        indices = np.where(correlation_matrix_ab[i, :] >= 0.3)[0]

        if indices.size <= 29:

            removed_rows.append(i)
            continue

        if indices.size > 1:

            for j in indices:
                indices3 = indices[indices != j]
                avg_value = np.mean(correlation_matrix_ab[indices3, j])
                correlation_matrix_H[i, j] = avg_value

        else:

            correlation_matrix_H[i, :] = 0
    for idx in removed_rows:
        d[idx, 0] = 0
    if removed_rows:
     correlation_matrix_H = np.delete(correlation_matrix_H, removed_rows, axis=0)
    H_csv_file = f'control_H.csv'
    pd.DataFrame(correlation_matrix_H).to_csv(H_csv_file, index=False, header=False)
    pd.DataFrame(d).to_csv('d.csv', index=False, header=False)
file_path = 'control.txt'
calculate_d(file_path)

def calculate_score(file_path_control, file_path_data, file_path_P):
    fcon = read_and_variance(file_path_control)
    f1 = read_and_variance(file_path_data)

    tau = np.abs(f1 - fcon)
    tau = tau / np.sum(tau)

    H1 = pd.read_csv(file_path_P, header=None).values
    H = H1
    n = H.shape[0]

    file_path3 = 'd.csv'
    datad = pd.read_csv(file_path3, sep='\t', header=None)
    d = datad.iloc[:, 0].values
    d2= d.reshape(-1, 1)
    d1 = np.where(np.sum(H, axis=1) == 0, 1, 0)
    d1 = d1.reshape(-1, 1)


    alpha = 0.85
    c = 0.05
    f = 0.1
    e = np.ones(n)
    s = np.ones((n, 1)) / n

    Entropy = np.ones((n, 1)) / n
    indices_to_compute = np.where(d2.flatten() == 1)[0]

    valid_entropies = []
    invalid_nodes = []

    for i in indices_to_compute:
        p_i = H[i, :]
        if np.any(p_i > 0):
            positive_probs = p_i[p_i > 0]
            entropy_i = -np.sum(positive_probs * np.log2(positive_probs))
            Entropy[i] = entropy_i
            valid_entropies.append(entropy_i)
        else:
            invalid_nodes.append(i)


    if invalid_nodes and valid_entropies:

        min_entropy = np.min(valid_entropies)
        for i in invalid_nodes:
            Entropy[i] = min_entropy

    tau2 = d2 * Entropy
    tau2 = tau2 / np.sum(tau2)
    prev_s = s.copy()
    while True:
        s = alpha *( (H + (d1 * e) / n).T).dot(s) + c * tau2 + f * tau

        if np.all(np.abs(s - prev_s) < 1e-8):
            break
        prev_s = s.copy()

    s = s.flatten()

    s_k = np.zeros((n, 2))
    s_k[:, 0] = np.arange(1, n + 1)
    s_k[:, 1] = s

    sorted_indices = np.argsort(-s_k[:, 1])
    s_k_sorted = s_k[sorted_indices, :]

    top_5_percent = int(0.05 * n)
    score = np.sum(s_k_sorted[:top_5_percent, 1])

    return score, s_k_sorted[:top_5_percent, 0]


prefix_data = "test2_data"
prefix_P = "P_"
suffix = ".txt"
control_file = "control.txt"


scores = []
top_species_indices = []
k = 17 #The number of single samples (time points) other than the control group
for i in range(1, k + 1):
    file_name_data = f"{prefix_data}{i}{suffix}"
    file_name_P = f"{prefix_P}{i}.csv"
    score,indices = calculate_score(control_file, file_name_data, file_name_P)
    scores.append(score)

    ## Output the node serial numbers of the top 5%
    # top_species_indices.append(indices)
# top_species_indices_array = np.array(top_species_indices).T
# top_species_indices_array_csv_file = 'top_scores_idE.csv'
# pd.DataFrame(top_species_indices_array).to_csv(top_species_indices_array_csv_file, index=False, header=False)


scores_array = np.array(scores).reshape(-1, 1)
scores_array2 = scores_array / scores_array[0]
np.savetxt("scores_array_10013.txt", scores_array2, fmt='%.6f')
print("Scores_array2 saved to scores_array2.txt")

