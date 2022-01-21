library(mice)
library(HardyWeinberg)

options(expressions = 100000)
args <- commandArgs(trailingOnly = T)
DATA <- read.table(file = args[1], sep = "\t", header = T, check.names = F)

th<-0.05
for (i in 1: nrow(DATA)) {
  rr<-DATA$`0/0`[i]
  ra<-DATA$`0/1`[i]
  aa<-DATA$`1/1`[i]
  sum<-rr + ra + aa
  p<-ifelse(sum == 0, 0, (rr * 2 + ra) / (2 * sum))
  q<-ifelse(sum == 0, 0, (aa * 2 + ra) / (2 * sum))
  if (p < th || q < th) {
    T<-HWExact(X = c(AA = rr, AB = ra, BB = aa))
  } else {
    T<-HWChisq(X = c(AA = rr, AB = ra, BB = aa))
  }
  DATA$freq0[i] = sprintf("%.7f", p)
  DATA$HWpvalue[i] = T$pval
}

write.table(x = DATA, file = args[1], quote = F, sep = "\t", row.names = F)