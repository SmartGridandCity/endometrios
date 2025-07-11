# =============================================================================
# 1. INSTALL & LOAD LIBRARIES
# =============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required <- c("affy","limma","ggplot2","randomForest","pheatmap","ggrepel",
              "caret","iml","data.table","dplyr","pROC",
              "VennDiagram","clusterProfiler","org.Hs.eg.db")
bioc_pkgs <- c("affy","limma","clusterProfiler","org.Hs.eg.db")
cran_pkgs <- setdiff(required, bioc_pkgs)

for (pkg in bioc_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
for (pkg in cran_pkgs)  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
invisible(lapply(required, library, character.only = TRUE))

# =============================================================================
# 2. SET PATHS AND READ CEL FILES
# =============================================================================
BASE <- normalizePath("C:/Users/guera/Documents/Python Scripts/R", winslash = "/")
DATA <- file.path(BASE, "data")
RESULTS <- file.path(BASE, "results_final")
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)

cel_files <- list.files(DATA, pattern = "\\.CEL$", ignore.case = TRUE)
stopifnot(length(cel_files) > 0)

affy_data <- ReadAffy(celfile.path = DATA)
eset <- rma(affy_data)
expr <- exprs(eset)
labels <- factor(c(rep("Eutopic", 9), rep("Ectopic", 9)))

# =============================================================================
# 3. QC PLOTS
# =============================================================================
png(file.path(RESULTS, "boxplot_expression.png"))
boxplot(expr, main = "Expression Distribution", las = 2)
dev.off()

top_genes <- names(sort(apply(expr, 1, var), decreasing = TRUE))[1:10]
annot <- data.frame(Group = labels); rownames(annot) <- colnames(expr)

png(file.path(RESULTS, "heatmap_top_genes.png"))
pheatmap(expr[top_genes, ], scale = "row", annotation_col = annot)
dev.off()

# =============================================================================
# 4. DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================
design <- model.matrix(~labels)
fit <- eBayes(lmFit(expr, design))
de_table <- topTable(fit, coef = 2, adjust = "fdr", number = Inf)
de_table$logFC <- fit$coefficients[, 2]
de_table$pval  <- fit$p.value[, 2]

sig <- subset(de_table, pval < 0.05 & abs(logFC) > 1)
features_basic <- rownames(sig)
write.csv(sig, file.path(RESULTS, "significant_genes.csv"))

# Volcano plots
png(file.path(RESULTS, "volcano_plot.png"))
ggplot(de_table, aes(x = logFC, y = -log10(pval), color = pval < 0.05 & abs(logFC) > 1)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "Volcano Plot", x = "LogFC", y = "-log10(p-value)")
dev.off()

png(file.path(RESULTS, "volcano_annotated.png"))
ggplot(de_table, aes(x = logFC, y = -log10(pval),
                     label = ifelse(pval < 0.05 & abs(logFC) > 1, rownames(de_table), ""))) +
  geom_point(aes(color = pval < 0.05 & abs(logFC) > 1), alpha = 0.6) +
  geom_text_repel(size = 3, max.overlaps = 10) + theme_minimal() +
  labs(title = "Volcano Plot w/ Labels")
dev.off()

# =============================================================================
# 5. FEATURE SELECTION METHODS
# =============================================================================

# ---- RFECV ----
X_basic <- as.data.frame(t(expr[features_basic, ]))
rfe_ctrl <- rfeControl(functions = rfFuncs, method = "repeatedcv", number = 10, repeats = 5)
rfecv_result <- rfe(x = X_basic, y = labels, sizes = c(5, 10, 20, 30, 50), rfeControl = rfe_ctrl)
features_rfecv <- rfecv_result$optVariables
writeLines(features_rfecv, file.path(RESULTS, "RFECV_features.txt"))

png(file.path(RESULTS, "RFECV_performance.png"))
plot(rfecv_result, type = c("g", "o"), main = "RFECV Accuracy")
dev.off()

# ---- Gini Importance ----
rf_model_gini <- randomForest(X_basic, labels, importance = TRUE)
rf_gini <- importance(rf_model_gini, type = 2)[, 1]
features_gini <- names(sort(rf_gini, decreasing = TRUE)[1:20])
writeLines(features_gini, file.path(RESULTS, "Gini_features.txt"))

# ---- Lasso Selection via glmnet ----
library(glmnet)
X_mat <- as.matrix(X_basic)
cv_fit <- cv.glmnet(X_mat, labels, alpha = 1, family = "binomial")
lasso_coef <- coef(cv_fit, s = "lambda.min")
features_lasso <- rownames(lasso_coef)[lasso_coef[, 1] != 0][-1]  # Remove intercept
writeLines(features_lasso, file.path(RESULTS, "Lasso_features.txt"))

# ---- Boruta Selection ----
if (!requireNamespace("Boruta", quietly = TRUE)) install.packages("Boruta")
library(Boruta)

boruta_result <- Boruta(x = X_basic, y = labels, doTrace = 0)
features_boruta <- getSelectedAttributes(boruta_result, withTentative = FALSE)
writeLines(features_boruta, file.path(RESULTS, "Boruta_features.txt"))

# ---- Filter-Based Selection: Top Differential t-statistics ----
top_stats <- rownames(de_table)[order(abs(de_table$logFC), decreasing = TRUE)][1:20]
features_filter <- intersect(top_stats, rownames(expr))
writeLines(features_filter, file.path(RESULTS, "Filter_features.txt"))


# ---- SHAP (Commented Out: unstable with small sample size) ----
# rf_model_shap <- randomForest(X_basic, labels)
# predictor <- Predictor$new(rf_model_shap, data = X_basic, y = labels)
# 
# shap_values <- list()
# for (i in 1:nrow(X_basic)) {
#   shap <- tryCatch(Shapley$new(predictor, x.interest = X_basic[i,, drop = FALSE]), error = function(e) NULL)
#   if (!is.null(shap)) shap_values[[i]] <- shap$results$feature
# }
# shap_freq <- sort(table(unlist(shap_values)), decreasing = TRUE)
# features_shap <- intersect(names(shap_freq)[1:20], colnames(X_basic))
# writeLines(as.character(features_shap), file.path(RESULTS, "SHAP_features.txt"))
# 
# shap_df <- data.frame(Gene = features_shap, Count = as.numeric(shap_freq[features_shap]))
# pdf(file.path(RESULTS, "SHAP_frequency_barplot.pdf"), width = 9, height = 6)
# p <- ggplot(shap_df, aes(x = reorder(Gene, Count), y = Count)) +
#   geom_bar(stat = "identity", fill = "tomato") +
#   coord_flip() +
#   labs(title = "SHAP Frequency Across Samples", x = "Gene", y = "Frequency") +
#   theme_minimal()
# print(p)
# dev.off()

# ---- Permutation Importance via vip (Commented Out: crashes due to small sample size or dispatch failure) ----
# if (!requireNamespace("vip", quietly = TRUE)) install.packages("vip")
# library(vip)
# library(caret)
# 
# rf_model_train <- train(
#   x = X_basic,
#   y = labels,
#   method = "rf",
#   trControl = trainControl(method = "none")
# )
# 
# rf_pred_wrapper <- function(object, newdata) {
#   predict(object, newdata = newdata, type = "prob")[, "Ectopic"]
# }
# 
# perm <- tryCatch({
#   vi_permute(
#     object = rf_model_train,
#     train = X_basic,
#     target = labels,
#     metric = "accuracy",
#     pred_wrapper = rf_pred_wrapper,
#     nsim = 10,
#     parallel = FALSE,
#     progress = FALSE
#   )
# }, error = function(e) {
#   message("Permutation importance failed: ", e$message)
#   NULL
# })
# 
# if (!is.null(perm)) {
#   features_perm <- perm$Variable[1:20]
#   writeLines(as.character(features_perm), file.path(RESULTS, "Permutation_features.txt"))
#   
#   pdf(file.path(RESULTS, "Permutation_importance_barplot.pdf"), width = 9, height = 6)
#   print(vip(perm, num_features = 20, bar_type = "ggplot"))
#   dev.off()
# }



# =============================================================================
# 6. COMBINE FEATURE SETS
# =============================================================================

features_list <- list(
  Basic   = features_basic,
  RFECV   = features_rfecv,
  Gini    = features_gini,
  Lasso   = features_lasso,
  Boruta  = features_boruta,
  Filter  = features_filter
)


if (exists("features_perm") && length(features_perm) > 0) {
  features_list$Permutation <- features_perm
}

features_list <- Filter(function(f) length(f) > 0, features_list)

feature_summary <- do.call(rbind, lapply(names(features_list), function(method) {
  data.frame(Method = method, Feature = features_list[[method]])
}))
write.csv(feature_summary, file.path(RESULTS, "FeatureSelection_Comparison.csv"), row.names = FALSE)

# =============================================================================
# 7. GENE INCLUSION DIAGRAM
# =============================================================================
# Build binary matrix of gene inclusion across methods
all_genes <- unique(unlist(features_list))
binary_matrix <- sapply(features_list, function(set) all_genes %in% set)
rownames(binary_matrix) <- all_genes

# Filter to genes selected by ???2 methods
selected_gene_mask <- rowSums(binary_matrix) >= 2
binary_matrix_filtered <- binary_matrix[selected_gene_mask, ]


if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)
binary_matrix_filtered <- binary_matrix_filtered[, sort(colnames(binary_matrix_filtered))]

# Convert logical matrix to numeric (0 = FALSE, 1 = TRUE)
binary_numeric <- apply(binary_matrix_filtered, 2, as.numeric)
rownames(binary_numeric) <- rownames(binary_matrix_filtered)

# Optional: sort methods alphabetically
binary_numeric <- binary_numeric[, sort(colnames(binary_numeric))]

# Plot heatmap
pheatmap(binary_numeric,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = c("grey95", "#1f78b4"),
         border_color = NA,
         main = "Feature Inclusion Across Selection Methods",
         fontsize_row = 6)

pdf(file.path(RESULTS, "FeatureInclusion_Heatmap.pdf"), width = 10, height = 10)
pheatmap(binary_numeric,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = c("grey95", "#1f78b4"),
         border_color = NA,
         main = "Feature Inclusion Across Selection Methods",
         fontsize_row = 6)
dev.off()


# =============================================================================
# 8. MODEL TRAINING & EVALUATION
# =============================================================================
X_all <- as.data.frame(t(expr))
rf_models <- list()
for (name in names(features_list)) {
  f <- intersect(features_list[[name]], colnames(X_all))
  if (length(f) > 0) {
    mod <- tryCatch(randomForest(X_all[, f], labels, ntree = 500), error = function(e) NULL)
    if (!is.null(mod)) rf_models[[name]] <- mod
  }
}
rf_models <- Filter(Negate(is.null), rf_models)

accuracy <- sapply(rf_models, function(m) mean(m$predicted == labels))
oob_error <- sapply(rf_models, function(m) m$err.rate[nrow(m$err.rate), 1])
comparison <- data.frame(Method = names(rf_models),
                         Accuracy = round(accuracy, 3),
                         OOB_Error = round(oob_error, 3))
write.csv(comparison, file.path(RESULTS, "ModelComparison.csv"), row.names = FALSE)

png(file.path(RESULTS, "OOB_Error_Comparison.png"), width = 900, height = 600)
print(
  ggplot(comparison, aes(x = Method, y = OOB_Error, fill = Method)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = OOB_Error), vjust = -0.5) +
    labs(title = "OOB Error Comparison", y = "OOB Error") +
    theme_minimal()
)
dev.off()

# =============================================================================
# 9. ROC & CONFUSION MATRICES
# =============================================================================
roc_results <- list()
for (name in names(rf_models)) {
  prob <- predict(rf_models[[name]], type = "prob")[, "Ectopic"]
  roc_results[[name]] <- roc(labels, prob, levels = rev(levels(labels)))
}
png(file.path(RESULTS, "ROC_Curve_Comparison.png"), width = 900, height = 700)
plot(roc_results[[1]], col = "red", lwd = 2, main = "ROC Curves")
for (i in 2:length(roc_results)) {
  plot(roc_results[[i]], add = TRUE, col = i + 1, lwd = 2)
}
legend("bottomright", legend = names(roc_results), col = 2:(length(roc_results) + 1), lwd = 2)
dev.off()

for (name in names(rf_models)) {
  cm <- confusionMatrix(predict(rf_models[[name]]), labels)
  sink(file.path(RESULTS, paste0("ConfusionMatrix_", name, ".txt")))
  cat("Confusion Matrix -", name, "\n\n")
  print(cm)
  sink()
}
