import pandas as pd

def collapse_sequences(input_file, output_file):
    # Загрузка данных из файла CSV
    data = pd.read_csv(input_file)

    # Преобразование типа данных столбца 'Frequency' в числовой формат (float)
    data['Frequency'] = data['Frequency'].astype(float)

    # Схлопывание одинаковых последовательностей в пределах каждой реплики и суммирование частот
    collapsed_data = data.groupby(['Sequence', 'Replica']).agg({'Frequency': 'sum'}).reset_index()

    # Округление значений частоты до одного знака после запятой
    collapsed_data['Frequency'] = collapsed_data['Frequency'].round(1)

    # Сохранение объединенных данных в файл CSV
    collapsed_data.to_csv(output_file, index=False)

# Пример использования функции
input_file = r'C:\Users\Роман\Documents\Science\Moscow\Seq\NGS\YFV_model (Summer_24)\Negative (beta).csv'
output_file = 'Collapsed_b.csv'
collapse_sequences(input_file, output_file)
