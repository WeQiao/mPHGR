import os
import re
import tempfile
import pandas as pd
import numpy as np


def run_mphgr(data_dir, data_type, d, Q, prefix_data, stages=None):
    os.chdir(data_dir)

    file_path2 = 'control.txt'
    datacon = pd.read_csv(file_path2, sep='\t', header=None)
    corr_control = datacon.T.corr()
    corr_control.to_numpy()

    data_type = str(data_type)
    P_folder = tempfile.TemporaryDirectory()
    file_indices = sorted(
        int(m.group(1))
        for file_name in os.listdir('.')
        if (m := re.fullmatch(rf'{prefix_data}(\d+)\.txt', file_name))
    )

    def calculate_and_save_matrices(file_path, index):
        data = pd.read_csv(file_path, sep='\t', header=None)
        correlation_matrix = data.T.corr()
        correlation_matrix_ab = np.absolute(correlation_matrix.to_numpy() - corr_control.to_numpy())
        correlation_matrix_H = np.zeros_like(correlation_matrix_ab)
        removed_rows = []

        for i in range(correlation_matrix_ab.shape[0]):
            indices = np.where(correlation_matrix_ab[i, :] >= Q)[0]

            if data_type == '1' and indices.size == 1:
                removed_rows.append(i)
                continue

            if indices.size > 1:
                for j in indices:
                    indices3 = indices[indices != j]
                    avg_value = np.mean(correlation_matrix_ab[indices3, j])
                    correlation_matrix_H[i, j] = avg_value
            else:
                correlation_matrix_H[i, :] = 0

        if data_type == '1' and removed_rows:
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
        P_csv_file = os.path.join(P_folder.name, f'P_{index}.csv')
        pd.DataFrame(P).to_csv(P_csv_file, index=False, header=False)

    for i in file_indices:
        file_name = f'{prefix_data}{i}.txt'
        calculate_and_save_matrices(file_name, i)

    def read_and_variance(file_name):
        data = np.loadtxt(file_name)
        variances = np.var(data, axis=1)
        variances = variances.reshape(-1, 1)
        return variances

    def calculate_d(file_path):
        data = pd.read_csv(file_path, sep='\t', header=None)
        correlation_matrix = data.T.corr()
        correlation_matrix_ab = np.absolute(correlation_matrix.to_numpy())
        d_vector = np.ones((correlation_matrix_ab.shape[0], 1))
        correlation_matrix_H = np.zeros_like(correlation_matrix_ab)
        removed_rows = []

        for i in range(correlation_matrix_ab.shape[0]):
            indices = np.where(correlation_matrix_ab[i, :] >= 0.3)[0]

            if indices.size < d:
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
            d_vector[idx, 0] = 0

        if data_type == '1':
            if removed_rows:
                correlation_matrix_H = np.delete(correlation_matrix_H, removed_rows, axis=0)

        pd.DataFrame(d_vector).to_csv('d.csv', index=False, header=False)

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
        d_vector = datad.iloc[:, 0].values
        d2 = d_vector.reshape(-1, 1)
        d1 = np.where(np.sum(H, axis=1) == 0, 1, 0)
        d1 = d1.reshape(-1, 1)

        alpha = 0.85
        c = 0.05
        f = 0.1
        e = np.ones(n)
        p = np.ones((n, 1)) / n
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

        prev_p = p.copy()
        while True:
            p = alpha * ((H + (d1 * e) / n).T).dot(p) + c * tau2 + f * tau
            if np.all(np.abs(p - prev_p) < 1e-8):
                break
            prev_p = p.copy()

        p = p.flatten()
        p_k = np.zeros((n, 2))
        p_k[:, 0] = np.arange(1, n + 1)
        p_k[:, 1] = p
        sorted_indices = np.argsort(-p_k[:, 1])
        p_k_sorted = p_k[sorted_indices, :]

        top_5_percent = int(0.05 * n)
        score = np.mean(p_k_sorted[:top_5_percent, 1])
        return score, p_k_sorted[:top_5_percent, 0]

    prefix_P = 'P_'
    suffix = '.txt'
    control_file = 'control.txt'
    scores = []
    top_species_indices = []

    for i in file_indices:
        file_name_data = f'{prefix_data}{i}{suffix}'
        file_name_P = os.path.join(P_folder.name, f'{prefix_P}{i}.csv')
        score, indices = calculate_score(control_file, file_name_data, file_name_P)
        scores.append(score)
        top_species_indices.append(indices)

    scores_array = np.array(scores).reshape(-1, 1)

    if data_type == '2':
        stage_scores = {}
        for stage, indices in stages.items():
            stage_scores[stage] = np.mean(
                [scores[file_indices.index(i)] for i in indices]
            )

        reference_score = next(iter(stage_scores.values()))
        for stage in stage_scores:
            stage_scores[stage] = stage_scores[stage] / reference_score

        column_vector = np.array(list(stage_scores.values())).reshape(-1, 1)
        np.savetxt('Score.csv', column_vector, fmt='%.6f', delimiter=',')

        import matplotlib.pyplot as plt
        labelx = list(stage_scores.keys())
        labely = list(stage_scores.values())
        fig, ax = plt.subplots(figsize=(10, 6))
        pic = ax.plot(labelx, labely, marker='o')
        # ax.set_title('The Scores for Gastric Cancer')
        ax.set_xlabel('Stage')
        ax.set_ylabel('Score')
        plt.show()

        P_folder.cleanup()
        return {
            'stage_scores': stage_scores,
        }

    scores_array2 = scores_array / scores_array[0]
    np.savetxt('Score.csv', scores_array2, fmt='%.6f', delimiter=',')
    P_folder.cleanup()
    return {
        'normalized_scores': scores_array2,
    }
