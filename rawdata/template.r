dataset <- read.table(file = "dataset.txt", as.is = TRUE, header = TRUE)
names(dataset) <- tolower(names(dataset))
save("dataset", file = "../data/dataset.rda", compress = "bzip2")
