# CDF-based Correction of Storm Tracker
Workflow and code
2024.09.19

## Workflow
- Collect raw data
- From level 0 to level 1: format unification.
- From level 1 to level 2: correct known errors such as missing values and impossible data.
- Perform CDF-based correction (and generate level 3 data) 
<img src="https://github.com/tingsyo/storm-tracker-calibration-published/blob/main/docs/images/images/ver01_fig01.jpg" />
Image source: [Developing High-Quality Field Program Sounding Datasets (Ciesielski et al., 2011)](https://journals.ametsoc.org/view/journals/bams/93/3/bams-d-11-00091_1.xml)

## Example: 2021-0803-12Z

- Files:
    - [ST_L0_2021080312_3014.csv](https://github.com/tingsyo/storm-tracker-calibration-published/blob/main/release/ver_01_20240919/ST_L0_2021080312_3014.csv)
    - [ST_L1_2021080312_3014.csv](https://github.com/tingsyo/storm-tracker-calibration-published/blob/main/release/ver_01_20240919/ST_L1_2021080312_3014.csv)
    - [ST_L2_2021080312_3014.csv](https://github.com/tingsyo/storm-tracker-calibration-published/blob/main/release/ver_01_20240919/ST_L2_2021080312_3014.csv)
    - [ST_L3_2021080312_3014.csv](https://github.com/tingsyo/storm-tracker-calibration-published/blob/main/release/ver_01_20240919/ST_L3_2021080312_3014.csv)
    - [ST_L2_2021080312_3014_cdf_corrected.csv](https://github.com/tingsyo/storm-tracker-calibration-published/blob/main/release/ver_01_20240919/ST_L2_2021080312_3014_cdf_corrected.csv)

<img src="https://github.com/tingsyo/storm-tracker-calibration-published/blob/main/docs/images/images/ver01_fig02.jpg" />


## Prerequisites of CDF-based Correction
- Python (free)
    - You may download Python through the official website or Anaconda (which provides pre-installed packages for various purposes).

- Python packages used (free):
    - NumPy: the fundamental package for scientific computing with Python
    - Pandas: pandas is a fast, powerful, flexible, and easy-to-use open-source data analysis and manipulation tool.

- **Note:** the two packages are pre-installed if you use Anaconda. Otherwise, you may install it with the commands:
``` bash
pip install numpy pandas
```


## Before Using CDF-based Correction
- Preprocess your ST observation data to level 2.
- Define the timestamp of the observation in `YYYYMMDDHH` (in UTC). For example, 2024091902 is the timestamp for 2024-09-19 10:00 AM Taipei Time.
- Define `dp0`, the ground check pressure reference. This value is the difference between the pressure measure by ST at launch `(P_st_0)` and the reference pressure sensor `(P_0)`. That is to say, `dp0 = P_st_0 – P_0`.


## Using CDF-based Correction
- In the directory of the code and model files:
    - perfrom_cdf_correction.py
    - cdf_T_correction.csv
    - cdf_RH_correction.csv

- In the command line, type:
``` bash
python perform_cdf_correction.py -i ST_L2_2021080312_3014.csv -t 2021080312 --dp0 1.5941
```
- The command results in `ST_L2_2021080312_3014_cdf_corrected.csv` and some statistics.
- You may type `python perform_cdf_correction.py -h` for more options.

