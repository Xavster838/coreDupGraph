'''
cluster mapping table and output annotations as tab delimited txt file.
'''
import numpy as np 
import pandas as pd
import networkx as nx 
from community import community_louvain

df = pd.read_csv(snakemake.input.mapping_stats , sep = "\t")
#make network
dup_net = nx.from_pandas_edgelist(df, source = "ref_loc_name", target="q_dup_name", edge_attr="score")
#get louvain clusterings.
comms = community_louvain.best_partition(dup_net)
#output to dataframe
out_df = pd.DataFrame(data = [comms.keys(), comms.values()]).transpose()
out_df.columns = ['dup', 'cluster']
out_df.to_csv(snakemake.output.cluster_annotation , sep = "\t", header = True, index = False)