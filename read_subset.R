t1 <- readRDS("subset1.rds")
head(t1)
tail(t1)

cat("Number of columns: ", ncol(t1))
print("\n")
cat("Number of rows: ", nrow(t1))

min_start <- min(t1$start)
max_start <- max(t1$start)

cat("\nMinimum start: ", min_start)
cat("\nMaximum start: ", max_start)
