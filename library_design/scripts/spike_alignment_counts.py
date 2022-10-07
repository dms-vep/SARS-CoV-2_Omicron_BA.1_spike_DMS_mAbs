import pandas as pd
import requests

table_url = snakemake.params.table_url
alignment_counts = snakemake.output.alignment_counts

print(f"mutation count dataframe from {table_url}")
html_text = requests.get(table_url, verify=False).content.decode("utf-8")
df = pd.read_html(html_text, encoding="utf-8")[0]

print(f"parsing mutation count dataframe ")
cols_of_interest = {
    "WildtypeAA": "wildtype",
    "Position": "site",
    "MutatedAA": "mutant",
    "#Occurrence": "count",
}

mut_counts = (
    df
    .query("MutatedAA.notnull()", engine="python")
    .query("WildtypeAA.notnull()", engine="python")
    .query("Position.notnull()", engine="python")
    .query("Protein == 'Spike'")
    .rename(columns=cols_of_interest)
    [cols_of_interest.values()]
    .assign(site=lambda x: x["site"].astype(int))
    .sort_values("count", ascending=False)
    .reset_index(drop=True)
)

mut_counts = mut_counts[(mut_counts.mutant != 'J') & (mut_counts.mutant != 'Z') & (mut_counts.mutant != 'B')]

#save mutation counts
print(f"saving {alignment_counts}")
mut_counts.to_csv(alignment_counts, index=False)