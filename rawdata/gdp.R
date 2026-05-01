# preso https://fred.stlouisfed.org/series/GDPC1 il 2026-05-01
gdp_full <- read.csv(file = "gdp.csv", as.is = TRUE, header = TRUE)
gdp_full$observation_date <- as.Date(gdp_full$observation_date)
gdpraw <- gdp_full$GDPC1
gdp <- ts(data = gdp_full$GDPC1, frequency=4, start=c(1947, 01))
save("gdp", file = "../data/gdp.rda", compress = "bzip2")
save("gdpraw", file = "../data/gdpraw.rda", compress = "bzip2")
