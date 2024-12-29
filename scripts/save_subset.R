t1 <- readRDS("../inputs/exampleForLukas.RDS")
head(t1)

cat("Number of columns: ", ncol(t1))
print("\n")
cat("Number of rows: ", nrow(t1))

filtered_t1 <- t1[t1$seqnames == 18, ]
sorted_t1 <- filtered_t1[order(filtered_t1$start), ]

subset1 <- sorted_t1[1:700000, ]

cat("Number of columns: ", ncol(subset1))
print("\n")
cat("Number of rows: ", nrow(subset1))

saveRDS(subset1, file = "../inputs/subset1.rds")
