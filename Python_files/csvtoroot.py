import pandas as pd
import uproot
import numpy as np

csv_file = "../CSV_files/Bplus_to_JpsiKplus_signal.csv"
root_file = "../Root_files/Bplus_to_JpsiKplus_signal.root"

df = pd.read_csv(csv_file)

# convert all columns to numeric
df = df.apply(pd.to_numeric, errors="coerce")
df = df.dropna()

# structured array (TTree-compatible)
dtype = [(col, "f8") for col in df.columns]
structured = np.zeros(len(df), dtype=dtype)

for col in df.columns:
    structured[col] = df[col].to_numpy()

with uproot.recreate(root_file) as f:
    f.mktree(
        "Events",
        {col: "float64" for col in df.columns}
    )
    f["Events"].extend(structured)

print("Written:", len(df), "events as TTree")