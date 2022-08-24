#!/datb1/wuyuxuan/anaconda3/envs/anaconda/bin/python

# conda activate anaconda
import pandas as pd

df = pd.read_csv("Zygote_merged_gch.bedgraph",sep=' ', header=None)
df_sum = df.groupby(0)[[3,4]].sum()
df_sum['propotion'] = df_sum[4] / (df_sum[3]+df_sum[4])
df_sum.to_csv('Zygote_methylation_fraction.csv')