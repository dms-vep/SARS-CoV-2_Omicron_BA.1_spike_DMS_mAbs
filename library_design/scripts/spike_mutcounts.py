import pandas as pd
import re

# get variables from `snakemake`
mutcounts_tsv = snakemake.input.mutcounts_tsv
spike_mutcounts_csv = snakemake.output.spike_mutcounts

annotation_df = pd.read_csv(mutcounts_tsv, sep="\t", low_memory=False)

# parse the annotations in **Spike** into a dictionary object
annotation_dict = {
    k: [sub.split(":")[1] for sub in v.split(";") if sub.split(":")[0] == "S"]
    for k, v in annotation_df.set_index("node_id")["aa_mutations"].to_dict().items()
}

# Convert dictionary to df with were each node mutation is in a row
spike_mutcounts = pd.DataFrame.from_dict(annotation_dict, orient='index')
spike_mutcounts.reset_index(level=0, inplace=True)
spike_mutcounts = spike_mutcounts.melt(id_vars=['index'], value_name='spike_mut')
spike_mutcounts = spike_mutcounts[~spike_mutcounts['spike_mut'].isnull()].drop(['variable', 'index'], 1)

#count spike mutations reoccuring on nodes
spike_mutcounts = spike_mutcounts.spike_mut.value_counts().rename_axis('spike_mut').reset_index(name='counts')
spike_mutcounts

#split spike_mut string
spike_mutcounts['wt'] = spike_mutcounts['spike_mut'].str[0]
spike_mutcounts['amino_acid'] = spike_mutcounts['spike_mut'].str[-1]
spike_mutcounts['site'] = spike_mutcounts['spike_mut'].str[1:-1]
spike_mutcounts = spike_mutcounts.loc[spike_mutcounts['amino_acid'] != '*']

spike_mutcounts = spike_mutcounts[['site','amino_acid', 'counts']].copy()
spike_mutcounts = spike_mutcounts.rename(columns={'counts':"n_mutations_to"})



print(f"Writing to {spike_mutcounts_csv}")
spike_mutcounts.to_csv(spike_mutcounts_csv, index=False)
