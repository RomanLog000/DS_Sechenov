import pandas as pd


student_data = pd.read_csv(r'data\students_performance.csv')
mask1 = student_data['race/ethnicity'] == 'group A'
mask2 = student_data['race/ethnicity'] == 'group C'
print(abs(student_data[mask1]['writing score'].median() - student_data[mask2]['writing score'].mean()))




# df.info()
# df.describe(include=['type'])
# df.select_dtypes(include=['type'])
# df.shape
# df.head()
# df.tail()
# df.astype(type)
# df.value_counts(normalize=True)
# df.insert[position, name, data]
# df.merge(df, type of junction, on='column')

# & - конъюнкция
# | - дизъюнкция


# Data aggregation

# shop	country	    name	    qty
# 347	Украина	    Киев	    NaN
# 427	РФ	        Самара	    3
# 707	Беларусь	Минск	    4
# 957	РФ  	    Иркутск	    2
# 437	РФ	        Москва	    1
# 345	Украина	    Киев	    NaN

# res.pivot_table(['qty'],['country'], aggfunc='sum', fill_value = 0)

# country	qty
# Беларусь	4
# РФ	    6
# Украина	0
