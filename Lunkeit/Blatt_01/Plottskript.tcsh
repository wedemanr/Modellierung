#!/bin/tcsh -evx
#
set FN = Evolutionsgleichung.pdf
#
cat > rplot.r << ENDE
input = read.table("Evolutionsgleichung.txt", sep=";")
pdf("${FN}")
time = input\$V1
temp = input\$V2
plot(time,temp, type='l', col="red",
main = "Evolutionsgleichung",
xlab = "Timestep",
ylab = "Temperature in °C")

dev.off()
q()
ENDE
#
R --vanilla < rplot.r
exit

