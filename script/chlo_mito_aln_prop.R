reads_len <- read.table('chlo_mito/Gossypium_trilobum.trimmedReads.chlo.mito.reads.length')

mito <- read.table('chlo_mito/Gossypium_trilobum.trimmedReads.mito.blast.out.aln.length')

chlo <- read.table('chlo_mito/Gossypium_trilobum.trimmedReads.chlo.blast.out.aln.length')

mito_aln <- merge(reads_len, mito, by.x=1, by.y=1)
mito_aln_prop <- mito_aln$V2.y/ mito_aln$V2.x

chlo_aln <- merge(reads_len, chlo, by.x=1, by.y=1)
chlo_aln_prop <- chlo_aln$V2.y/ chlo_aln$V2.x

pdf('chlo_mito/chlo_mito_aln_prop.pdf', hei=4)
par(mfrow = c(1,2))
hist(mito_aln_prop)
hist(chlo_aln_prop)
dev.off()

