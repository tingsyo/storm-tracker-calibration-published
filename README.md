# Storm Tracker Calibration

This repository hosted the code and documentations for the calibration of the upper-air radiosonde instrument, “Storm Tracker” (ST). The work has been reviewed and published on [Journal of the Meteorological Society of Japan](https://www.jstage.jst.go.jp/article/jmsj/103/5/103_2025-029/_article) in 2025 (DOI: [https://doi.org/10.2151/jmsj.2025-029](https://doi.org/10.2151/jmsj.2025-029))

The details of the calibration algorithms can be found in the oublished paper. Here we host the following tools:

- `sample_data`: 10 samples of the co-launch dataset. Data in different preprocessing stages are tored as `raw`, `L0`, `L1`, `L2`, and `paired`.
- `docs/CDF_based_calibration/correction_tool/`: All 663 colaunches profiles used in the study.
- `CDF_based_calibration/correction_tool/`: the CDF-based correction tools.
- `CDF_based_calibration/model_builder/`: `python` scripts to build CDF models from paired data.


## Citation
```
Hung-Chi KUO, Ting-Shuo YO, Hungjui YU, Shih-Hao SU, Ching-Hwang LIU, Po-Hsiung LIN, Data Quality Control and Calibration for Mini-Radiosonde System “Storm Tracker” in Taiwan, Journal of the Meteorological Society of Japan. Ser. II, 2025, Volume 103, Issue 5, Pages 573-593
```

## Links
- Online ISSN: 2186-9057, 
- Print ISSN: 0026-1165, 
- Article (open access): https://doi.org/10.2151/jmsj.2025-029

