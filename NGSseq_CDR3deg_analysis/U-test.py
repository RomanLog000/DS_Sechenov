import pandas as pd
from scipy.stats import mannwhitneyu

# Загрузка данных из CSV файла
file_path = r'C:\Users\Роман\Documents\Science\Moscow\Seq\NGS\YFV_model (April_24)\YF_LLW_Alpha\Pattern_8.csv'
data = pd.read_csv(file_path)

# Удаление консенсусного клонотипа
consensus_sequence = 'CASSLAPGATNEKLFF'
data = data[data['sequence'] != consensus_sequence]

# Получение уникальных значений последовательностей клонотипов
unique_sequences = data['sequence'].unique()

# Создание пустого списка для хранения значимых клонотипов
significant_clonotypes = []

# Создание списка для хранения p-value для каждого клонотипа
p_values = []

# Цикл для проведения U-теста Манна-Уитни для каждого клонотипа
for sequence in unique_sequences:
    # Получение данных для текущего клонотипа
    sequence_data_pos = data[(data['sequence'] == sequence) & (data['fraction'] == 'pos')]['frequency']
    sequence_data_neg = data[(data['sequence'] == sequence) & (data['fraction'] == 'neg')]['frequency']

    # Проверка на количество значений и исключение клонотипов с недостаточным размером выборки
    if len(sequence_data_pos) > 1 and len(sequence_data_neg) > 1:
        # Выполнение U-теста Манна-Уитни
        statistic, p_value = mannwhitneyu(sequence_data_pos, sequence_data_neg)

        # Оценка статистической значимости
        alpha = 0.05
        if p_value < alpha:
            # Добавление значимого клонотипа в список
            significant_clonotypes.append(sequence)
            # Добавление p-value в список
            p_values.append(p_value)

# Создание DataFrame для хранения результатов
results_df = pd.DataFrame(columns=['Sequence', 'Specific (pos)', 'Non-specific (neg)', 'p-value'])

# Заполнение DataFrame данными о значимых клонотипах
for sequence, p_value in zip(significant_clonotypes, p_values):
    # Получение данных для текущего клонотипа
    sequence_data_pos = data[(data['sequence'] == sequence) & (data['fraction'] == 'pos')]['frequency']
    sequence_data_neg = data[(data['sequence'] == sequence) & (data['fraction'] == 'neg')]['frequency']

    # Округление частот клонотипов до одного знака после запятой
    mean_frequency_pos = round(sequence_data_pos.mean(), 1)
    mean_frequency_neg = round(sequence_data_neg.mean(), 1)

    # Добавление строки в DataFrame с информацией о клонотипе, его частотах и p-value
    results_df = results_df.append({'Sequence': sequence,
                                    'Specific (pos)': mean_frequency_pos,
                                    'Non-specific (neg)': mean_frequency_neg,
                                    'p-value': round(p_value, 2)},  # Округление p-value до двух знаков
                                    ignore_index=True)

# Сохранение результатов в CSV файл
results_df.to_csv('significant_clonotypes.csv', index=False)
