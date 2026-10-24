import sys
import pandas as pd

# python3 this.py tpm.txt_list TPM output.txt

# reading list
with open(sys.argv[1], 'r') as f:
    file_list = [line.strip() for line in f]

# set dict
gene_data = {}

# 
for file_name in file_list:
    print(f"Processing {file_name}...")  # 
    for chunk in pd.read_csv(file_name, sep='\t', usecols=["Gene ID", sys.argv[2]], chunksize=1000):
        for _, row in chunk.iterrows():
            gene_id = row["Gene ID"]
            value = row[sys.argv[2]]
            if gene_id not in gene_data:
                gene_data[gene_id] = {}
            gene_data[gene_id][file_name] = value

#  DataFrame
merged_df = pd.DataFrame.from_dict(gene_data, orient="index").reset_index()
merged_df.rename(columns={"index": "Gene ID"}, inplace=True)

# 
output_file = sys.argv[3]
merged_df.to_csv(output_file, sep='\t', index=False)
print(f"Merge completed. Output saved to {output_file}.")
