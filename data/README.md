# Clinical example data

`cav.csv` and `psor.csv` are unmodified CSV exports of the public example objects in `msm` 1.8.2, by Christopher Jackson. They contain 2,846 cardiac allograft vasculopathy observations and 806 psoriatic arthritis observations, respectively. The numeric patient identifiers are the identifiers supplied in those public example objects.

Source: [the msm package](https://chjackson.github.io/msm/), including its [CAV](https://chjackson.github.io/msm/reference/cav.html) and [PsA](https://chjackson.github.io/msm/reference/psor.html) documentation. The original clinical references are Sharples et al. (2003) and Gladman and Farewell (1999), as cited in the manuscript.

The source package is distributed under GPL version 2 or later. The included `GPL-2` gives the licence text for these data exports. `Export_Data.R` recreates the CSVs from the package without applying patient exclusions or changing values. The analysis programs retain the package objects as their inputs.
