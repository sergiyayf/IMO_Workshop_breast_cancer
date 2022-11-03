import pandas as pd

df = pd.read_csv('Initialization_For_PhysiCell.csv')

print(df)

df['x'] *= 15
df['y'] *= 15

df2 = df.replace(4,7)
df3 = df2.replace(3,4)
df4 = df3.replace(0,3)
df5 = df4.replace(2,5)
df6 = df5.replace(1,2)

df7 = df6.replace(7,0)

df7.to_csv('cells.csv',header = False, index = False)
print(df7)