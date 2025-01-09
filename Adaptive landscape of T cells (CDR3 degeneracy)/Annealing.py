import pandas as pd

# Загрузка данных из загруженного файла
file_path = r'C:\Users\Роман\Documents\Science\Moscow\Seq\NGS\YFV_model (Summer_24)\Processed_a.csv'  # Укажите путь к загруженному файлу
data = pd.read_csv(file_path)

# Параметры таблицы
replicas = 2
patterns = 8

# Подготовка списков столбцов
pos_freq_columns = [f'pos_freq_{i}_{pattern}' for i in range(1, replicas + 1) for pattern in range(1, patterns + 1)]
pos_seq_columns = [f'pos_seq_{i}_{pattern}' for i in range(1, replicas + 1) for pattern in range(1, patterns + 1)]
neg_freq_columns = [f'neg_freq_{i}_{pattern}' for i in range(1, replicas + 1) for pattern in range(1, patterns + 1)]
neg_seq_columns = [f'neg_seq_{i}_{pattern}' for i in range(1, replicas + 1) for pattern in range(1, patterns + 1)]

# Создание DataFrame для хранения соответствия "сиквенс-частота" для каждой реплики и каждой фракции
sequences_freq = pd.DataFrame(columns=['sequence', 'frequency', 'replica', 'fraction', 'pattern'])

# Итерация по всем комбинациям реплик и паттернов для pos и neg фракций
for i in range(replicas):
    for j in range(patterns):
        # Добавление данных pos фракции
        pos_pattern_data = pd.DataFrame({
            'sequence': data[pos_seq_columns[i * patterns + j]],
            'frequency': data[pos_freq_columns[i * patterns + j]],
            'replica': i + 1,
            'fraction': 'pos',
            'pattern': j + 1
        })
        sequences_freq = pd.concat([sequences_freq, pos_pattern_data], ignore_index=True)

        # Добавление данных neg фракции
        neg_pattern_data = pd.DataFrame({
            'sequence': data[neg_seq_columns[i * patterns + j]],
            'frequency': data[neg_freq_columns[i * patterns + j]],
            'replica': i + 1,
            'fraction': 'neg',
            'pattern': j + 1
        })
        sequences_freq = pd.concat([sequences_freq, neg_pattern_data], ignore_index=True)

# Сохранение результата в новый CSV файл
output_file_path = 'Annealed_a.csv'
sequences_freq.to_csv(output_file_path, index=False)

# Вывод первых нескольких строк для проверки
print(sequences_freq.head())
