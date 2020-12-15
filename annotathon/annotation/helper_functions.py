import pandas as pd


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
