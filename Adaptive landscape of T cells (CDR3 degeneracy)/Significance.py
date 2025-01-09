import pandas as pd
import matplotlib.pyplot as plt

# Загрузка данных из файла CSV
file_path = r'C:\Users\Роман\Documents\Science\Moscow\Seq\NGS\CMV_model (Summer_23)_Exp2\Annealed_data (mut-1).csv'  # Укажите путь к вашему файлу CSV
data = pd.read_csv(file_path)

# Усреднение значений частот клонотипов по репликам для каждой фракции
pos_data = data[data['fraction'] == 'pos'].groupby('sequence')['frequency'].mean()
neg_data = data[data['fraction'] == 'neg'].groupby('sequence')['frequency'].mean()

# Удаление значений NaN
pos_data.dropna(inplace=True)
neg_data.dropna(inplace=True)

# Исключение консенсусного клонотипа
consensus_sequence = 'CASSLAPGATNEKLFF'
if consensus_sequence in pos_data.index:
    pos_data.drop(consensus_sequence, inplace=True)
if consensus_sequence in neg_data.index:
    neg_data.drop(consensus_sequence, inplace=True)

# Соотнесение сиквенсов между pos и neg фракциями
common_sequences = pos_data.index.intersection(neg_data.index)
pos_data = pos_data.loc[common_sequences]
neg_data = neg_data.loc[common_sequences]

# Загрузка результатов значимых клонотипов
significant_clonotypes = pd.read_csv('significant_clonotypes.csv')['Sequence']

# Построение scatter plot
plt.figure(figsize=(10, 6))

# Отрисовка всех точек
plt.scatter(pos_data, neg_data, alpha=0.7, s=3, color='blue', label='Non-significant')

# Выделение значимых клонотипов
for sequence in significant_clonotypes:
    if sequence in common_sequences:
        idx = common_sequences.get_loc(sequence)
        plt.scatter(pos_data.iloc[idx], neg_data.iloc[idx], color='red', s=10)

plt.title('Fraction Distribution')
plt.xlabel('Specific')
plt.ylabel('Non-Specific')
plt.legend(['Non-significant', 'Significant'], loc='upper right')
plt.grid(True)
plt.show()

# Вывод таблицы с значениями частот для каждого сиквенса
df = pd.DataFrame({'Seq': common_sequences, 'Specific': pos_data.values, 'Non-Specific': neg_data.values})
csv_file_path = 'Fraction_scores (mut-1).csv'  # Укажите путь к файлу CSV для сохранения таблицы
df.to_csv(csv_file_path, index=False)
