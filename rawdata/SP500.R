sp500 <- read.table(file = "SP500.txt", as.is = TRUE, header = TRUE)[, 1]
save("sp500", file = "../data/sp500.rda", compress = "bzip2")
