input = read.table("Evolutionsgleichung.txt", sep=";")
pdf("Evolutionsgleichung.pdf")
time = input$V1
temp = input$V2
plot(time,temp, type='l', col="red",
main = "Evolutionsgleichung",
xlab = "Timestep",
ylab = "Temperature in °C")

dev.off()
q()
