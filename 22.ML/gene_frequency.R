library(ggplot2)

d <- read.table("model.genes.txt", header=T, sep="\t", stringsAsFactors=F, check.names=F)
total_models <- length(unique(d$algorithm))
cat("有效模型总数:", total_models, "\n")

freq_table <- as.data.frame(table(d$features))
colnames(freq_table) <- c("gene", "frequency")

freq_table <- freq_table[order(-freq_table$frequency), ]
freq_table$percentage <- freq_table$frequency / total_models * 100
freq_table$pct_label <- sprintf("%.0f%%", freq_table$percentage)

write.table(freq_table, "model.gene_frequency.txt", sep="\t", row.names=F, quote=F, col.names=T)

high_freq_genes <- freq_table[freq_table$frequency >= 10, ]
write.table(high_freq_genes, "model.high_freq_genes.txt", sep="\t", row.names=F, quote=F, col.names=T)

cat("=== 基因频率统计结果 ===\n")
cat("总基因数: ", nrow(freq_table), "\n")
cat("被 >= 50% 模型选中: ", sum(freq_table$percentage >= 50), "\n")
cat("被 >= 30% 模型选中: ", sum(freq_table$percentage >= 30), "\n")
cat("被 >= 10% 模型选中: ", sum(freq_table$percentage >= 10), "\n")
cat("最高频率: ", max(freq_table$frequency), " (", max(freq_table$percentage), "%)\n")
cat("输出文件: model.gene_frequency.txt, model.high_freq_genes.txt, model.gene_frequency.pdf\n")

top_n <- min(30, nrow(freq_table))
top_genes <- freq_table[1:top_n, ]
top_genes$gene <- factor(top_genes$gene, levels = rev(top_genes$gene))

p <- ggplot(top_genes, aes(x = gene, y = frequency)) +
  geom_bar(stat = "identity", fill = "#3B8FC2", width = 0.7) +
  geom_text(aes(label = pct_label), hjust = -0.1, size = 3.2) +
  coord_flip() +
  labs(title = "Gene selection frequency across 113 models (Top 30)",
       subtitle = NULL,
       x = "Gene", y = "Frequency (# of models)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_family = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("model.gene_frequency.pdf", p, width = 8, height = max(6, top_n * 0.25), device = "pdf")
cat("绘图完成: model.gene_frequency.pdf\n")
