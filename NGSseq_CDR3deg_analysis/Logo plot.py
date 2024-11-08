import logomaker
import pandas as pd
import matplotlib.pyplot as plt


def create_motif_dict(seq, seq_type="prot", weight=1):
    if seq_type == 'dna':
        variants = ["A", "T", "G", "C"]
    elif seq_type == 'rna':
        variants = ["A", "U", "G", "C"]
    elif seq_type == 'prot':
        variants = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    else:
        print(f"Unknown seq type '{seq_type}'")
        return
    motif_dict = {var: [weight if seq[i] == var else 0 for i in range(len(seq))] for var in variants}
    return motif_dict


def sum_motif_dicts(list_of_dicts):
    variants = list(list_of_dicts[0].keys())
    seq_len = len(list_of_dicts[0][variants[0]])
    sum_motif_dict = {var: [0 for _ in range(seq_len)] for var in variants}
    for motif_dict in list_of_dicts:
        for var in variants:
            sum_motif_dict[var] = [a + b for a, b in zip(sum_motif_dict[var], motif_dict[var])]
    return sum_motif_dict


def normalize_motif_dict(motif_dict):
    variants = list(motif_dict.keys())
    seq_len = len(motif_dict[variants[0]])
    sum_list = [sum([motif_dict[var][i] for var in variants]) for i in range(seq_len)]
    norm_motif_dict = dict()
    for var in variants:
        norm_motif_dict[var] = [a / b if b != 0 else 0 for a, b in zip(motif_dict[var], sum_list)]
    return norm_motif_dict


def get_logo_for_sequences(data, seq_type="prot"):
    motif_dicts = []
    for i, r in data.iterrows():
        seq = r['Sequence']
        if pd.isna(seq):
            continue
        weight = r['Frequency']
        if weight == 0:
            continue
        motif_dict = create_motif_dict(seq, seq_type=seq_type, weight=weight)
        motif_dicts.append(motif_dict)

    if not motif_dicts:
        print("No valid sequences found.")
        return

    motif_dict_sum = sum_motif_dicts(motif_dicts)
    motif_dict_sum = normalize_motif_dict(motif_dict_sum)

    info_matrix = pd.DataFrame(motif_dict_sum)
    logo = logomaker.Logo(info_matrix, stack_order='big_on_top')
    logo.ax.set_ylabel('Weight', labelpad=-1, fontsize=10)
    logo.ax.set_title("Specific", fontsize=12)
    plt.show()


# Загрузка данных
data_filename = r'C:\Users\Роман\Documents\Science\Moscow\Seq\NGS\YFV_model (Summer_24)\Union (Alpha_pos).csv'
data = pd.read_csv(data_filename)

# Суммирование частот для одинаковых последовательностей по всем репликам
grouped_data = data.groupby('Sequence').agg({'Frequency': 'sum'}).reset_index()

# Построение взвешенного logo plot
get_logo_for_sequences(grouped_data, seq_type="prot")
