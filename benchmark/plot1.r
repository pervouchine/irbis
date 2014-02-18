library(Hmisc)
library(ggplot2)

cmd_args = commandArgs();

data = read.delim(cmd_args[4], header=F)

colnames(data) = c('id','Threshold','Length','value')

pdf(cmd_args[5],width=9, height=8,onefile=T)

data$Threshold = as.factor(data$Threshold)
data$Length = as.factor(data$Length)

data1 = data[data$id==0,]
data2 = data[data$id>0,]

se <- function(x){1.64*sd(x)/sqrt(length(x))}

aggregate(data2$value,list(Threshold=data2$Threshold,Length=data2$Length), FUN=mean) -> data3a
aggregate(data2$value,list(Threshold=data2$Threshold,Length=data2$Length), FUN=se  ) -> data3b
merge(data3a, data3b, by=c('Threshold','Length')) -> data3
colnames(data3)=c('Threshold','Length','mean','se')

merge(data1, data3, by=c('Threshold','Length')) -> data4
data4$fdr = with(data4, mean/value)
data4$fdr.SE = with(data4, se/value)

q = theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),legend.text = element_text(size=16),axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))


df = data4[,c('Threshold','Length','fdr','fdr.SE')]
p <- ggplot(df, aes(x=Threshold, y=fdr)) + geom_line(aes(group=Length)) + geom_point(aes(shape=Length),size=4)
p + geom_errorbar(width=.1, aes(ymin=fdr-fdr.SE, ymax=fdr+fdr.SE, shape=Length)) + scale_shape_manual(values=c(0,5,6,15,17)) + xlab(expression(paste("Conservation threshold (",t[2],")",sep=""))) + ylab("FDR") + q


p <- ggplot(df, aes(x=Length, y=fdr)) + geom_line(aes(group=Threshold)) + geom_point(aes(shape=Threshold),size=4)
p + geom_errorbar(width=.1, aes(ymin=fdr-fdr.SE, ymax=fdr+fdr.SE, shape=Threshold)) + scale_shape_manual(values=c(0,5,6,8,15,17))    + xlab("Minimum structure length (L)") + ylab("FDR") + q


dev.off()

