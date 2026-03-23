# Build CDF Correction Models

CDF-based Probability matching, also known as histogram matching or quantile mapping, is a statistical technique used to adjust the distribution of a dataset (e.g., a forecast distribution) to match that of another dataset (e.g., an observed distribution). The primary objective of this method is not to directly correct individual data points but to ensure that the overall statistical properties, such as the frequency of occurrence of specific values, match between the two datasets. In radiosonde observation, CDF-based probability matching is commonly used as a QC tool to ensure data quality consistency for field campaigns ([Nuret et al. 2008](https://www.jstage.jst.go.jp/article/jmsj/103/5/103_2025-029/_html/-char/en#ref22); [Ciesielski et al. 2009](https://www.jstage.jst.go.jp/article/jmsj/103/5/103_2025-029/_html/-char/en#ref1)).

More details about the model construction procedure can be found in [our published work](https://www.jstage.jst.go.jp/article/jmsj/103/5/103_2025-029/_html/-char/en). Here is the script, `build_cdf_corrections.py`, to build the model from paired co-launch data.

## Prerequisites

- Have a compatable python interpreter (developed and tested with 3.9 ~ 3.12)
- Install required packages: `numpy`, `pandas`, and `matplotlib`.

## Data Preprocessing

Make sure you already complete the [data preprocessing] procedure and put the 

Aggregate all paired colaunch data from `.csv` to one single dataframe using the following code:

``` python
from datetime import datetime
import os
import pandas as pd

DATAPATH = 'YOUR_PAIRED_DATA_DIRECTORY/' # The directory of paired data in CSV format
files = []
# List all files in the DATAPATH
for root, dirs, files in os.walk(DATAPATH):
    for file_name in files:
        # Construct the full file path
        full_path = os.path.join(root, file_name)
        files.append(full_path)

# Aggregate files
data = []
for f in paired_files:
    tmp = pd.read_csv(f)
    data.append(tmp)

df = pd.concat(data, axis=0)

# Write the aggregate data to a single picke file (binary)
df.to_pickle("data.pkl") 
```

Note that the 10 paired files in the `sample_data/paired/` can be aggregated but they are not sufficient to build the CDF models. If you want to test the script, you may use [this file](https://drive.google.com/file/d/1dQAINO47POVuiS_iLUJEYVwROwPMe6X4/view?usp=sharing) (~112MB) which is similar to what we used in our study.

With the aggregated data, then use the script:

``` bash
> python build_cdf_corrections.py -i [AGGREGATED_DATA] -o [OUTPUT_DIRECTORY]
``` 

The output directory contains 4 files:

- `cdf_T_correction.csv`: the model for temperature (T) correction.
- `cdf_RH_correction.csv`: the model for relative humidity (RH) correction.
- `summary_cdf_T_correction.csv`: statistics of temperature (T) correction.
- `summary_cdf_RH_correction.csv`: statistic of relative humidity (RH) correction.

