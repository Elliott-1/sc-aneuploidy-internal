library(Matrix)

full_mtx <- readRDS("/net/noble/vol1/home/es1/aneuploidy/data/full_mtx.RDS")
# Write the matrix to a Matrix Market file
writeMM(full_mtx, file = "/net/noble/vol1/home/es1/aneuploidy/data/full_mtx.mtx")