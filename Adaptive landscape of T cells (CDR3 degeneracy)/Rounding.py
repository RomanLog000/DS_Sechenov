import pandas as pd

def round_frequency(input_file, output_file):
    # Загрузка данных из файла CSV
    data = pd.read_csv(input_file)

    # Округление значений столбцов 'frequency_pos' и 'frequency_neg' до одного знака после запятой
    data['frequency_pos'] = data['frequency_pos'].round(1)
    data['frequency_neg'] = data['frequency_neg'].round(1)

    # Сохранение округленных данных в файл CSV
    data.to_csv(output_file, index=False)

# Пример использования функции
input_file = r'C:\Users\Роман\Documents\Science\Moscow\Seq\NGS\YFV_model (April_24)\YF_LLW_Beta\Collapsed_data (pattern-1).csv'
output_file = 'rounded_data.csv'
round_frequency(input_file, output_file)
