import pandas as pd

# Загрузка данных
pos_file = r"C:\Users\Роман\Documents\Science\Moscow\Ph.D\Olya\YFV_cdr3_sequences_pos_A_filtered.csv"
neg_file = r"C:\Users\Роман\Documents\Science\Moscow\Ph.D\Olya\YFV_cdr3_sequences_neg_A_filtered.csv"

# Чтение файлов
pos_data = pd.read_csv(pos_file, header=None, names=["Sequence"])
neg_data = pd.read_csv(neg_file, header=None, names=["Sequence"])

# Удаление общих строк
unique_pos_data = pos_data[~pos_data["Sequence"].isin(neg_data["Sequence"])]

# Сохранение результата
unique_pos_data.to_csv("YFV_cdr3_sequences_pos_A_filtered_unique.csv", index=False, header=False)
