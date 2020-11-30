import pandas as pd


def parse_alignment_descriptions(handle: str) -> pd.DataFrame:
    """ Parse the horrible file provided by the NCBI """
    # skip header line
    rows = pd.read_csv(handle, skiprows=[0], header=None)
    # add columns
    rows.columns = (
        pd.read_csv(handle, nrows=0).columns.to_series().apply(lambda x: x.strip())
    )
    rows.loc[:, "Accession"] = rows["Accession"].apply(
        lambda x: x.split(",")[1].replace(")", "").replace('"', "")
    )

    return rows
