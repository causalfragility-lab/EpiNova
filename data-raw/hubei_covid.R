## Run this script once to regenerate data/hubei_covid.rda
## Place in data-raw/ and call: source("data-raw/hubei_covid.R")

hubei_covid <- list(
  NI = c(41,41,41,45,62,131,200,270,375,444,549,729,
         1052,1423,2714,3554,4903,5806,7153,9074,11177,
         13522,16678,19665,22112,24953,27100,29631,31728,33366),
  RI = c(1,1,7,10,14,20,25,31,34,45,55,71,94,121,152,213,
         252,345,417,561,650,811,1017,1261,1485,1917,2260,
         2725,3284,3754),
  N  = 58.5e6,
  begin_date = "2020-01-13",
  end_date   = "2020-02-11",
  description = paste(
    "Daily cumulative confirmed (NI) and removed (RI) COVID-19 cases",
    "in Hubei Province, China, from Jan 13 to Feb 11 2020.",
    "Source: DXY.cn. Population N = 58.5 million."
  )
)

usethis::use_data(hubei_covid, overwrite = TRUE, compress = "xz")
