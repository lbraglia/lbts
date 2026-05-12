load("fMRI.RData")
FMRI <- Y
save("FMRI", file = "../data/fmri.rda", compress = "bzip2")
