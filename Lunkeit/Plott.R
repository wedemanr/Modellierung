in.data <- read.table(file = "Evolutionsgleichung.txt",sep = ";")

plot(in.data$V1, in.data$V2, type = "l", xlab = "Timesteps", ylab = "Temperature in °C")
