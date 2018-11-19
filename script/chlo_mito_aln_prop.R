reads_len <- read.table('chlo_mito/Gossypium_trilobum.trimmedReads.chlo.mito.reads.length')
colnames(reads_len) <- c('id', 'tot.length')

mito <- read.table('chlo_mito/Gossypium_trilobum.trimmedReads.mito.blast.out.aln.length')
colnames(mito) <- c('id', 'mito.aln.length')

chlo <- read.table('chlo_mito/Gossypium_trilobum.trimmedReads.chlo.blast.out.aln.length')
colnames(chlo) <- c('id', 'chlo.aln.length')

result <- merge(reads_len, mito, by.x = 1, by.y = 1, all.x = T)
result <- merge(result, chlo, by.x = 1, by.y = 1, all.x = T)
result[is.na(result)] <- 0

result$mito.aln.percent <- result$mito.aln.length / result$tot.length * 100
result$chlo.aln.percent <- result$chlo.aln.length / result$tot.length * 100

pdf('chlo_mito/chlo_mito_aln_prop.pdf', hei=4)
barplot(t(as.matrix(result[order(-result$mito.aln.percent, result$chlo.aln.percent), c('mito.aln.percent', 'chlo.aln.percent')])), beside = T, border = NA, col = c('red', 'blue'), xaxt = 'n', ylab='Percentage of Read Length (%)')
legend('topright', c('mitochondria', 'chloroplast'), fill = c('red', 'blue'), border = NA, bty='n')
#
par(mfrow = c(1,2))
hist(result$mito.aln.percent, main = 'mitochondria', xlab = 'Percentage of Read Length (%)')
hist(result$chlo.aln.percent, main = 'chloroplast', xlab = 'Percentage of Read Length (%)')
dev.off()

