from Bio import motifs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import scipy as sp
from itertools import combinations
import seaborn as sns
import statistics as stat

# функция для создания списка всех последовательностей, входящих в датасет с учетом числа их повторов:

def sequence_list_maker(df):
    aaSeqCDR3_series = df['aaSeqCDR3']
    readcount_series = df['readCount']
    aaSeqCDR3_len_list = [len(seq) for seq in aaSeqCDR3_series]
    len_mode = stat.mode(aaSeqCDR3_len_list)
    for index in range(len(aaSeqCDR3_series)):
        if len(aaSeqCDR3_series[index]) != len_mode:
            aaSeqCDR3_series = aaSeqCDR3_series.drop(index=index)
            readcount_series = readcount_series.drop(index=index)
    new_aaSeqCDR3_list = []
    for num, seq in zip(readcount_series, aaSeqCDR3_series):
        i = int(num)
        while i != 0:
            new_aaSeqCDR3_list.append(seq)
            i -= 1
    return new_aaSeqCDR3_list

# функция для объединения списков всех последовательностей из разных датасетов (для объединения последовательностей из реплик):

def merge_lists(*sequence_lists):
    summary_seq_list = []
    for seq_list in sequence_lists:
        aaSeqCDR3_len_list = [len(seq) for seq in seq_list]
        summary_seq_list += seq_list
    for i in summary_seq_list:
        if len(i) != stat.mode(aaSeqCDR3_len_list):
            summary_seq_list.remove(i)
    return summary_seq_list

# функция для создания и визуализации PWM:

def pmw_maker(sequence_list):
    sequence_list = motifs.create(sequence_list, alphabet='GLYSEQUDNFAKRHCVPWIMT_*')
    pwm = sequence_list.counts.normalize(pseudocounts=0.5)
    color_code = {'G': '#ff0000',
                  'L': '#ff4d00',
                  'Y': '#ff7514',
                  'S': '#ff8800',
                  'E': '#ffb300',
                  'Q': '#fff44f',
                  'U': '#ccff00',
                  'D': '#7fff00',
                  'N': '#0bda51',
                  'F': '#008000',
                  'A': '#2effaf',
                  'K': '#1dacd6',
                  'R': '#122faa',
                  'H': '#735ec4',
                  'C': '#9966cc',
                  'V': '#7442c8',
                  'P': '#9400d3',
                  'W': '#440066',
                  'I': '#240935',
                  'M': '#641349',
                  'T': '#ca2c92',
                  '_': '#000000',
                  '*': '#525252'}
    pwm_df = pd.DataFrame(pwm)
    pwm_logo = logomaker.Logo(pwm_df,
                              stack_order='small_on_top',
                              color_scheme=color_code,
                              shade_below=.5,
                              fade_below=.5,
                              font_name='Consolas')

    pwm_logo.style_spines(visible=False)
    pwm_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    pwm_logo.ax.set_ylabel("Aminoacid probability", labelpad=-1, fontsize=10)
    pwm_logo.ax.xaxis.set_ticks_position('none')
    pwm_logo.ax.xaxis.set_tick_params(pad=-1)
    return plt.show()

# функция для записи аминокислот и их частот в вырожденных позициях в словарь:

def degenereate_pos_find(sequence_list):
    sequence_list = motifs.create(sequence_list, alphabet='GLYSEQUDNFAKRHCVPWIMT_*')
    pwm = sequence_list.counts.normalize(pseudocounts=0.5)
    pwm_df = pd.DataFrame(pwm)
    pwm_df_transposed = pwm_df.T
    list_aa = list(pwm_df_transposed.index)
    columns_dict = {}
    for column in pwm_df_transposed:
        column_series = pwm_df_transposed[column]
        aa_runner = 0
        new_keys_aa_prob_dict = []
        new_values_aa_prob_dict = []
        for value in column_series:
            value = round(value, 3)
            if 0.009 < value < 0.99:
                new_keys_aa_prob_dict.append(list_aa[aa_runner])
                new_values_aa_prob_dict.append(value)
            aa_runner += 1
        degenerate_aa_prob_dict = {i: j for (i, j) in zip(new_keys_aa_prob_dict, new_values_aa_prob_dict)}
        columns_dict.update({column: degenerate_aa_prob_dict})
    for i in list(columns_dict):
        if not columns_dict[i]:
            del columns_dict[i]
    return columns_dict

# функция для проверки корреляции между частотами допустимых аминок-т в вырожденных позициях по критерию независимости хи-квадрат (критерий Пирсона):

def chi_square_calc(dict1, dict2):
    chi2_analysis_dict = {}
    for i, j in zip(dict1, dict2):
        row_i = dict1[i]
        row_j = dict2[j]
        a = []
        b = []
        if len(row_i) != len(row_j):
            for key1 in row_i.keys():
                if key1 not in row_j.keys():
                    b.append((key1, float(0.000)))
            for key2 in row_j.keys():
                if key2 not in row_i.keys():
                    a.append((key2, float(0.000)))
        row_i.update(a)
        row_i.update(b)
        row_i = dict(sorted(row_i.items()))
        row_j = dict(sorted(row_j.items()))
        d1_row_i_key_list = []
        d1_row_i_val_list = []
        for key in row_i.keys():
            d1_row_i_key_list.append(key)
        for value in row_i.values():
            d1_row_i_val_list.append(value)
        '''
         создали списки с ключами и со значениями для каждого словаря с (аминокислотами: частотами),
         содержащегося в первом словаре с (позициями: частотным распределением аминокислот) 
        '''
        d2_row_j_key_list = []
        d2_row_j_val_list = []
        for key in row_j.keys():
            d2_row_j_key_list.append(key)
        for value in row_j.values():
            d2_row_j_val_list.append(value)
        '''
         то же самое для второго словаря
        '''
        pos_i_aa_freq_lists = [d1_row_i_val_list, d2_row_j_val_list]
        chi2_position_i = sp.stats.chi2_contingency(pos_i_aa_freq_lists)
        chi2_analysis_dict.update({i: chi2_position_i})
    return chi2_analysis_dict

# функция для построения тепловых карт условных вероятностей:

def internal_correlation_hm(seq_list):
    deg_pos_aa_dict = degenereate_pos_find(seq_list)
    positions_df = pd.DataFrame(deg_pos_aa_dict)
    positions_df.fillna(0, inplace=True)
    a = positions_df.columns
    for pair in combinations(a, 2):
        current_position_1, current_position_2 = pair[0], pair[1]
        mid_dict = {}
        for current_aa_1 in positions_df.index:
            init_dict = {}
            for current_aa_2 in positions_df.index:
                counter = 0
                for i in seq_list:
                    if ((i[current_position_1] == current_aa_1) and (i[current_position_2] == current_aa_2)):
                        counter += 1
                init_dict.update({current_aa_2: counter})
            mid_dict.update({current_aa_1: init_dict})
        aa_df = pd.DataFrame(mid_dict)
        pos_df_series = positions_df[current_position_1]
        alt_series = positions_df[current_position_2]
        for col in aa_df.columns:
            aa_df[col] = aa_df[col]/len(seq_list)
        aa_corr_df = aa_df.div(pos_df_series, axis=1)
        aa_corr_df.fillna(0, inplace=True)
        for col in aa_corr_df.columns:
            aa_corr_df[col].replace([np.inf, -np.inf], alt_series[col], inplace=True)

        figsize = (8, 8)
        fontsize_annotation = 8
        fig, ax = plt.subplots(figsize=figsize)
        plt.title(f'Positions {current_position_1}-{current_position_2} Correlation Matrix',
                  fontsize=16)
        aa_corr_heatmap = sns.heatmap(aa_corr_df,
                                      annot=True,
                                      square=True,
                                      fmt='.1f',
                                      annot_kws={'size': str(fontsize_annotation)},
                                      cmap=sns.color_palette("viridis", as_cmap=True))
        aa_corr_heatmap.set(xlabel=f'Aminoacids in {current_position_1} position',
                            ylabel=f'Aminoacids in {current_position_2} position')
        aa_corr_heatmap.set_xticklabels(aa_corr_heatmap.get_xmajorticklabels(), fontsize=12)
        aa_corr_heatmap.set_yticklabels(aa_corr_heatmap.get_xmajorticklabels(), fontsize=12)
        aa_corr_heatmap.xaxis.tick_top()
        aa_corr_heatmap.xaxis.set_label_position('top')
        aa_corr_heatmap.yaxis.tick_left()
        plt.show()
    return 0

YFV_deg_a_pos_1 = pd.read_csv('/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/MiXCR/YF_deg_alpha_pos_01.clones_TRAD.tsv', sep='\t')

YFV_deg_a_neg_1 = pd.read_csv('/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/MiXCR/YF_deg_alpha_neg_01.clones_TRAD.tsv', sep='\t')
YFV_deg_a_neg_2 = pd.read_csv('/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/MiXCR/YF_deg_alpha_neg_02.clones_TRAD.tsv', sep='\t')

YFV_deg_b_pos_1 = pd.read_csv('/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/MiXCR/YF_deg_beta_pos_01.clones_TRB.tsv', sep='\t')
YFV_deg_b_pos_2 = pd.read_csv('/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/MiXCR/YF_deg_beta_pos_02.clones_TRB.tsv', sep='\t')

YFV_deg_b_neg_1 = pd.read_csv('/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/MiXCR/YF_deg_beta_neg_01.clones_TRB.tsv', sep='\t')
YFV_deg_b_neg_2 = pd.read_csv('/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/MiXCR/YF_deg_beta_neg_02.clones_TRB.tsv', sep='\t')

YFV_deg_a_pos_seq_list = merge_lists(sequence_list_maker(YFV_deg_a_pos_1))

YFV_deg_a_neg_seq_list = merge_lists(sequence_list_maker(YFV_deg_a_neg_1),
                                     sequence_list_maker(YFV_deg_a_neg_2))

YFV_deg_b_pos_seq_list = merge_lists(sequence_list_maker(YFV_deg_b_pos_1),
                                     sequence_list_maker(YFV_deg_b_pos_2))

YFV_deg_b_neg_seq_list = merge_lists(sequence_list_maker(YFV_deg_b_neg_1),
                                     sequence_list_maker(YFV_deg_b_neg_2))

YFV_deg_a_pos_deg_pos = degenereate_pos_find(YFV_deg_a_pos_seq_list)
YFV_deg_a_pos_deg_pos_df = pd.DataFrame(YFV_deg_a_pos_deg_pos)
YFV_deg_a_pos_deg_pos_df.fillna(0, inplace=True)
print(YFV_deg_a_pos_deg_pos_df)

YFV_deg_a_neg_deg_pos = degenereate_pos_find(YFV_deg_a_neg_seq_list)
YFV_deg_a_neg_deg_pos_df = pd.DataFrame(YFV_deg_a_neg_deg_pos)
YFV_deg_a_neg_deg_pos_df.fillna(0, inplace=True)
print(YFV_deg_a_neg_deg_pos_df)

YFV_deg_b_pos_deg_pos = degenereate_pos_find(YFV_deg_b_pos_seq_list)
YFV_deg_b_pos_deg_pos_df = pd.DataFrame(YFV_deg_b_pos_deg_pos)
YFV_deg_b_pos_deg_pos_df.fillna(0, inplace=True)
print(YFV_deg_b_pos_deg_pos_df)

YFV_deg_b_neg_deg_pos = degenereate_pos_find(YFV_deg_b_neg_seq_list)
YFV_deg_b_neg_deg_pos_df = pd.DataFrame(YFV_deg_b_neg_deg_pos)
YFV_deg_b_neg_deg_pos_df.fillna(0, inplace=True)
print(YFV_deg_b_neg_deg_pos_df)

a = internal_correlation_hm(YFV_deg_a_pos_seq_list)
b = internal_correlation_hm(YFV_deg_a_neg_seq_list)
c = internal_correlation_hm(YFV_deg_b_pos_seq_list)
d = internal_correlation_hm(YFV_deg_b_neg_seq_list)
