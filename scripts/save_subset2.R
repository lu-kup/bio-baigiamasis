t1 <- readRDS("../inputs/exampleForLukas.RDS")
head(t1)

cat("Number of columns: ", ncol(t1))
print("")
cat("Number of rows: ", nrow(t1))

filtered_t1 <- t1[t1$seqnames == 18, ]
sorted_t1 <- filtered_t1[order(filtered_t1$start), ]

cat("Number of columns: ", ncol(sorted_t1))
print("")
cat("Number of rows: ", nrow(sorted_t1))

saveRDS(sorted_t1, file = "../inputs/subset2.rds")
