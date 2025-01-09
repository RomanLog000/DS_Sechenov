import pandas as pd

# Загрузка данных из CSV файла
df = pd.read_csv('/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/communities (beta).csv')

# Подсчет количества узлов в каждом сообществе
community_sizes = df.groupby('Community').size().reset_index(name='Size')

# Сохранение результата в новый CSV файл
output_path = '/Users/Роман/Documents/Science/Moscow/Seq/NGS/YFV_model (Summer_24)/community_size.csv'
community_sizes.to_csv(output_path, index=False)
