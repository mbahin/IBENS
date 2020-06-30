import pandas as pd
import numpy as np
from dask import dataframe as dd
from dask.multiprocessing import get

df = pd.DataFrame(np.array([np.arange(1, 99999), np.arange(1, 99999)]).T, columns=["Val1", "Val2"])

df2 = dd.from_pandas(df, npartitions=20).map_partitions(lambda df: df.apply(lambda row: row.Val1 * row.Val2, axis=1)).compute(scheduler=get)
print(df2.head())
