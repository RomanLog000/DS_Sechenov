import pandas as pd
import matplotlib.pyplot as plt

# Загрузка данных из файла CSV
file_path = '/Users/Роман/Documents/Science/Moscow/Seq/NGS/CMV_model (Summer_23)/Exp2/Collapsed_data (mut-2).csv'
data = pd.read_csv(file_path)

plt.figure(figsize=(10, 6))
plt.scatter(data['frequency_pos'], data['frequency_neg'], alpha=0.7, s=10)
plt.xlabel('Specific')
plt.ylabel('Non-Specific')
plt.axline((0, 0), slope=1, color='red', label='by slope')
plt.grid(True)
plt.show()
