import pandas as pd


countries_df = pd.DataFrame(
    data = [
        ['Англия', 56.29, 133396],
        ['Канада', 38.05, 9984670],
        ['США', 322.28, 9826630],
        ['Россия', 146.24, 17125191],
        ['Украина', 45.5, 603628],
        ['Беларусь', 9.5, 207600],
        ['Казахстан', 17.04, 2724902]
    ],
    columns= ['country', 'population', 'square'],
    index = ['UK', 'CA', 'US', 'RU', 'UA', 'BY', 'KZ']
)


countries_df.to_csv('data/countries_df.csv', index=False, sep=';')
melb_data = pd.read_csv('data/melb_data.csv', sep=',')

melb_data['Car'] = melb_data['Car'].astype('int64')
melb_data['Bedroom'] = melb_data['Bedroom'].astype('int64')
melb_data['Bathroom'] = melb_data['Bathroom'].astype('int64')
melb_data['Propertycount'] = melb_data['Propertycount'].astype('int64')
melb_data['YearBuilt'] = melb_data['YearBuilt'].astype('int64')

filtered_data = melb_data[melb_data['Price'] < 10**6]
filtered_data = filtered_data[(filtered_data['Rooms'] > 5) | (filtered_data['YearBuilt'] > 2015)]
average_price = filtered_data['Price'].mean()
print(average_price)

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
