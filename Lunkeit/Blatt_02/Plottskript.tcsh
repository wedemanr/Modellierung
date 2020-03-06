#!/bin/tcsh -evx
#
set FN = Evolutionsgleichung.pdf
#
cat > rplot.r << ENDE
in1 = read.table("Evolutionsgleichung.txt", sep=";")
in2 = read.table("Output.txt", sep=";")
pdf("${FN}")
t1 = in1\$V1
T1 = in1\$V2
t2 = in2\$V1
T2 = in2\$V2
plot(t2,T2, type='l', col="blue",
main = "Evolutionsgleichung",
xlab = "Time",
ylab = "Temperature in °C", xlim = c(0,200))
lines(t1,T1, type="l2",col="red")
legend("bottomright",legend=c("Analytisch","Numerisch"),col=c("red","blue"), lty=1:1, cex=0.8)
dev.off()
q()
ENDE
#
R --vanilla < rplot.r
exit

