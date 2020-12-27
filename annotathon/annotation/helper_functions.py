import pathlib

from typing import List, Dict, NoReturn, Any, Callable, Tuple, Optional, Union

from Bio import SeqIO as seqio, SearchIO as searchio, Entrez as entrez

import pandas as pd
from ..utils.customobjs import Path as path


def add_unkown_description_tag(description: pd.Series) -> pd.Series:
    """ If no description is provided, add [unkown species] tag to description """
    no_species = (
        description.apply(lambda x: x if len(x.split("[")) == 1 else np.nan)
        .dropna()
        .index
    )

    return description[no_species].apply(
        lambda x: f"{x} [unknown species]" if len(x.split("[")) == 1 else x
    )


def add_function(description: pd.Series) -> pd.Series:
    """ Add a function column based on a description column """
    return description.apply(lambda x: x.split("[")[0])


def add_species(description: pd.Series) -> pd.Series:
    """ Add a species column based on a description column """
    return description.apply(
        lambda x: x.split("[")[1].replace("]", "")
        if len(x.split("[")) == 2
        else x.split("[")[0]
    )


def add_taxonomy(
    df: pd.DataFrame,
    file: Union[str, pathlib.Path],
    fformat: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """

    Params:
          df : a pandas.DataFrame, created by parsing a
               NCBI BLAST "Alignment-Descriptions" file.
               This means that the dataframe is expected
               to contain a `Accession` column.

        file : a file containing the taxonomy info

     fformat : The format of the file, defaults to "gb" (genbank)
               For more info on this subject please read BioPython's documentation:
               https://biopython.readthedocs.io/en/latest/api/Bio.SeqIO.html?highlight=Seqio%20formats#file-formats

    Returns:
        A dataframe containing a new column "taxonomy" with taxonomy info.
    """

    assert isinstance(df, pd.DataFrame), f"param `df` has invalid type {type(df)}"
    assert isinstance(file, str) or isinstance(
        file, pathlib.Path
    ), f"param `cache_file` has invalid type {type(file)}"
    assert (
        "Accession" in df.columns
    ), f"DataFrame `df` has no column `Accession` which is needed for this function"

    if isinstance(file, str):
        file = path(file)

    if not file.exists():
        raise FileNotFoundError(f"{file.absolute().as_posix()}")

    # Begin processing
    _df = df.copy()

    # Read taxonomy info
    with open(file, "r") as f:
        seq_dict = {sq.id: sq for sq in seqio.parse(f, format=fformat)}

    # Create accession lookup, taxonomy-yeilding function
    _tax_by_id = (
        lambda x: ":".join(seq_dict[x].annotations["taxonomy"])
        if x in seq_dict.keys()
        else ""
    )

    # Add the taxonomy column
    _df.loc[:, "taxonomy"] = ""
    _df.loc[:, "taxonomy"] = _df.Accession.apply(_tax_by_id)

    # Count the number of entries on each one
    n_newcols = _df.taxonomy.apply(lambda x: np.nan if not x else x).dropna().shape[0]
    n_oldcols = df.shape[0]

    # Let the user now not all records got taxonomy info
    if (n_newcols != n_oldcols) and verbose:
        print("Some entries have no taxonomy info!")
        print(f"Original number of entries = {n_oldcols}")
        print(f"Entries with non-null taxonomy info = {n_newcols}")

    return _df


def download_and_cache_genbank(
    accessions: List[str],
    cache_file: Union[str, pathlib.Path],
    efetch_kw: Optional[Dict[str, str]] = None,
    overwrite: Optional[bool] = None,
    verbose: Optional[bool] = True,
) -> bool:
    """
    Download and cache data via Biopython - NCBI's API

    This was thought as a way of downloading data directly from the NCBI, in order to create
    a unified interface to perform multicriteria filtering and selection of hit table results,
    also caled "Alignment descriptions.csv", when performing a BLAST search.

    When setting `rettype` in `efetch_kw` param, make sure that you use a valid format.
    For more info on this subject please read BioPython's documentation:
    https://biopython.readthedocs.io/en/latest/api/Bio.SeqIO.html?highlight=Seqio%20formats#file-formats

    It is recommended to use "gb" (shorthand for genbank) as `rettype` as this format preserves the most
    information. You could opt for "fasta", but this would leave you with only the sequence and
    a not so descriptive name.

    Parameters:
        accessions: a list of strings, accession identifiers.

        cache_file: either a string or a pathlib.Path instance.
                    This will be the file used to store the results.

         efetch_kw: A dictionnary containing valid keyword arguments
                    to call Bio.Entrez.efetch(). See `help(Bio.Entrez.efetch)`
                    for more details.
                    Defaults to:
                    _efetch_kw = {
                        "db": "protein",
                        "rettype": "gb",
                        "retmode": "text"
                    }

         overwrite: Boolean indicating if it is ok to overwrite the cache file.
                    Defaults to True.

           verbose: Print status messages whilst executing the function ?
                    Defaults to True.

    Returns:
        True  : if `cache_file` was successfully created and has size > 0 (not blank)

        False : otherwise
    """
    _efetch_kw = {"db": "protein", "rettype": "gb", "retmode": "text"}
    assert isinstance(
        accessions, list
    ), f"param `accessions` has invalid type {type(accessions)}"
    assert isinstance(cache_file, str) or isinstance(
        cache_file, pathlib.Path
    ), f"param `cache_file` has invalid type {type(cache_file)}"

    efetch_kw = efetch_kw or _efetch_kw
    overwrite = overwrite or False

    if isinstance(cache_file, str):
        cache_file = path(cache_file)

    if cache_file.exists() and not overwrite:
        raise FileExistsError(
            f"You are trying to overwrite {cache_file.absolute().as_posix()}"
        )

    if (not cache_file.exists()) or overwrite:

        if verbose:
            print(f"Number of accesions to fetch : {len(accessions)}")
            print("Querying entrez and fetching results...")
            print("Be patient, this might take a while...")
        with entrez.efetch(id=accessions, **efetch_kw) as in_handle:
            sequences = seqio.parse(in_handle, format=efetch_kw["rettype"])
            if verbose:
                print(f"Finished quering for {len(accessions)} accession numbers")
            with open(cache_file, "w") as out_handle:
                if verbose:
                    print(
                        f"Writting entries to cache file {cache_file.absolute().as_posix()}"
                    )
                seqio.write(sequences, out_handle, format=efetch_kw["rettype"])

        # os.stat_result(st_mode, st_ino, st_dev, st_nlink,
        # st_uid, st_gid, st_size, st_atime, st_mtime, st_ctime)
        if cache_file.exists():
            if cache_file.stat()[6] > 0:
                return True
            else:
                return False
        else:
            return False
