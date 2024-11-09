import pandas as pd

# Загрузка данных из исходного файла
file_path = r'C:\Users\Роман\Documents\Science\Moscow\Seq\NGS\YFV_model (Summer_24)\Assembled_b.csv'  # Укажите путь к вашему исходному файлу
data = pd.read_csv(file_path)

# Создание списка столбцов с аминокислотными последовательностями
sequence_columns = [col for col in data.columns if 'seq_' in col]

# Копируем данные, чтобы не изменять исходные
filtered_data = data.copy()

# Фильтрация строк с аминокислотными последовательностями длиной 13 и без инделов для каждого столбца с последовательностями
for column in sequence_columns:
    # Оставляем только те строки в столбце, которые соответствуют условию
    filtered_data[column] = filtered_data[column].apply(
        lambda x: x if pd.notna(x) and isinstance(x, str) and len(x) == 13 and not any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in x) else pd.NA
    )

# Сохранение отфильтрованных данных в новый CSV файл
filtered_file_path = 'Processed_data (beta).csv'  # Укажите путь для сохранения нового файла
filtered_data.to_csv(filtered_file_path, index=False)

# Вывод информации о количестве строк и столбцов после фильтрации
print("Количество строк после фильтрации:", len(filtered_data))
print("Количество столбцов после фильтрации:", len(filtered_data.columns))
