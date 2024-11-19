import pandas as pd


student_data = pd.read_csv(r'students_performance\students_performance.csv')
print(student_data.loc[155, 'writing score'])

# print(melb_data.info())

# df.info()
# df.describe(include=['type'])
# select_dtypes(include=['type'])
# df.shape
# df.head()
# df.tail()
# df.astype(type)
# value_counts(normalize=True)

# & - конъюнкция
# | - дизъюнкция
