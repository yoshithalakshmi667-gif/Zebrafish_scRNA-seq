library(Matrix)

expr_file <- "exprMatrix.tsv.gz.tmp"

# =========================================================
# 1. Read cell names
# =========================================================

cat("Reading header...\n")

con <- gzfile(expr_file, "rt")

header <- readLines(con, n = 1)

close(con)

cell_names <- strsplit(
    header,
    "\t",
    fixed = TRUE
)[[1]][-1]

cat("Cells:", length(cell_names), "\n")


# =========================================================
# 2. Read expression matrix
# =========================================================

cat("\nReading expression matrix...\n")

con <- gzfile(expr_file, "rt")

# Skip header
readLines(con, n = 1)

gene_names <- character(0)

row_list <- list()
col_list <- list()
value_list <- list()

gene_index <- 0

while (TRUE) {

    lines <- readLines(con, n = 500)

    if (length(lines) == 0) {
        break
    }

    for (line in lines) {

        # Skip blank lines
        if (nchar(line) == 0) {
            next
        }

        gene_index <- gene_index + 1

        fields <- strsplit(
            line,
            "\t",
            fixed = TRUE
        )[[1]]

        # Gene name
        gene_names[gene_index] <- fields[1]

        # Expression values
        values <- as.numeric(fields[-1])

        # Non-zero values
        nz <- which(values != 0)

        if (length(nz) > 0) {

            row_list[[gene_index]] <- rep.int(
                gene_index,
                length(nz)
            )

            col_list[[gene_index]] <- nz

            value_list[[gene_index]] <- values[nz]
        }

        # Progress
        if (gene_index %% 500 == 0) {

            cat(
                "\rGenes processed:",
                gene_index
            )

            flush.console()
        }
    }
}

close(con)

cat("\n\nGenes processed:", gene_index, "\n")
cat("Cells:", length(cell_names), "\n")


# =========================================================
# 3. Combine non-zero entries
# =========================================================

cat("\nCombining non-zero entries...\n")

i <- unlist(
    row_list,
    use.names = FALSE
)

j <- unlist(
    col_list,
    use.names = FALSE
)

x <- unlist(
    value_list,
    use.names = FALSE
)

cat(
    "Non-zero entries:",
    length(x),
    "\n"
)


# =========================================================
# 4. Create sparse matrix
# =========================================================

cat("\nCreating sparse matrix...\n")

expr_matrix <- sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(
        gene_index,
        length(cell_names)
    ),
    dimnames = list(
        gene_names,
        cell_names
    )
)


# =========================================================
# 5. Verify
# =========================================================

cat("\n====================================\n")
cat("Matrix successfully created!\n")
cat("====================================\n")

cat("Genes:", nrow(expr_matrix), "\n")
cat("Cells:", ncol(expr_matrix), "\n")
cat("Non-zero entries:", length(expr_matrix@x), "\n")
cat("Matrix class:", class(expr_matrix), "\n")

cat("====================================\n")

# Check dimensions
print(dim(expr_matrix))

# Check first 5 genes
cat("\nFirst 5 genes:\n")
print(head(rownames(expr_matrix), 5))

# Check first 5 cells
cat("\nFirst 5 cells:\n")
print(head(colnames(expr_matrix), 5))

# Check small section
cat("\nFirst 5 genes × first 5 cells:\n")
print(expr_matrix[1:5, 1:5])

saveRDS(
    expr_matrix,
    file = "zebrafish_expression_sparse.rds"
)

cat("Expression matrix saved successfully!\n")
