library(RColorBrewer)
library(scales)
col_set <- brewer.pal(n = 3, name = 'Set1')

argv <- commandArgs(T)
input <- argv[1]
input2 <- argv[2]
output <- paste0(input, '.pdf')

bt <- read.table(input)
bt_lst <- split(bt, bt$V1)
dfm <- read.table(input2, row.names = 1)

pdf(output, hei=5)
for (n in 1:length(bt_lst)) {
  m <- bt_lst[[n]]
  plot(0, xlim = c(0, max(m[7:10])), ylim = c(0, 10), type = 'n', yaxt='n', ylab='', xlab='', main = 'alignment length > 10000 bp')
  for(i in 1:nrow(m)) {
    if (m[i, 4] > 10000) {
      polygon(x = c(m[i, c(7,9,10,8)]), y = c(8,2,2,8), col = alpha(col_set[1], 0.8))
    }
  }
  text(0, 8.5, paste0(m[1, 1], ' (', dfm[m[1, 1], ], ' bp)'), pos = 4)
  text(0, 1.5, m[1, 2], pos = 4)
}
dev.off()

