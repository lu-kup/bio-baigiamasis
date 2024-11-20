t1 <- readRDS("exampleForLukas.RDS")
head(t1)

cat("Number of columns: ", ncol(t1))
print("\n")
cat("Number of rows: ", nrow(t1))

subset1=t1[1:700000,]

cat("Number of columns: ", ncol(subset1))
print("\n")
cat("Number of rows: ", nrow(subset1))

saveRDS(subset1, file = "subset1.rds")
