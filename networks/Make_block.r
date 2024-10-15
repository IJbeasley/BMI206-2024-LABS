rm(list=ls())

VEGAS1 <- read.table("HT.pvals.out", header=T, as.is=T)
VEGAS2 <- read.table("MS.pvals.out", header=T, as.is=T)
##### Sort VEGAS output by chromosome and position
VEGAS1 <- VEGAS1[order(VEGAS1$Chr,VEGAS1$Start,VEGAS1$Stop),]
VEGAS2 <- VEGAS2[order(VEGAS2$Chr,VEGAS2$Start,VEGAS2$Stop),]



##Make block column
p.BlockDefine.threshold <- 0.001
VEGAS1$block <- NA
current.block <- 0
for(current.row in 1:nrow(VEGAS1)){
    ## move on if pval is not good
    if(VEGAS1[current.row,"GenePvalue"] >= p.BlockDefine.threshold) next
    ## if first row, or new block, then bump to next block
    if(current.row == 1 ||
         VEGAS1[current.row - 1,"GenePvalue"] >= p.BlockDefine.threshold)
      current.block <- current.block + 1
    ## put the block in the matrix
    VEGAS1[current.row,"block"] <- current.block
} ## 
rm(current.row,current.block)
  

VEGAS2$block <- NA
current.block <- 0
for(current.row in 1:nrow(VEGAS2)){
    ## move on if pval is no good
    if(VEGAS2[current.row,"GenePvalue"] >= p.BlockDefine.threshold) next
    ## if first row, or new block, then bump to next block
    if(current.row == 1 ||
         VEGAS2[current.row - 1,"GenePvalue"] >= p.BlockDefine.threshold)
      current.block <- current.block + 1
    ## put the block in the matrix
    VEGAS2[current.row,"block"] <- current.block
} ##
rm(current.row,current.block)

# Create a data frame of blocks from the two VEGAS at matching gene names 
VEGAS1 = VEGAS1[!duplicated(VEGAS1$Gene),]
VEGAS2 = VEGAS2[!duplicated(VEGAS2$Gene),]
rownames(VEGAS1) = VEGAS1$Gene
rownames(VEGAS2) = VEGAS2$Gene
blocks = data.frame("VEGAS1_block" = VEGAS1$block, "VEGAS2_block"= VEGAS2[rownames(VEGAS1), "block"])

# Plot segments where there is a block
png("VEGAS_blocks.png")
k = 0.1
plot(NA, NA, xlim=c(0,4), ylim=c(0,14000), ylab = "Numer of Genes", xlab = "", bty="n", xaxt="n", yaxt="n")
axis(side=1, at = c(0,1,2))
axis(side=2, las=2)
mtext(text = "HT", side=1, at=0.5, col="blue")
mtext(text = "MS", side=1, at=1.5, col="blue")
for(i in 1:nrow(blocks))
{
    if(is.na(blocks$VEGAS1_block[i]))  segments(0,i, 1,i, col=NA, lwd = k)
    else segments(0,i, 1,i, col="red", lwd=k)
}

for(i in 1:nrow(blocks))
{
    if(is.na(blocks$VEGAS2_block[i])) segments(1,i, 2,i, col=NA, lwd =k)
    else segments(1,i, 2,i, col="red", lwd = k)
}
dev.off()

