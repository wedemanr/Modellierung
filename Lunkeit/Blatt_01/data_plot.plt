set xlabel "t"
set ylabel "T"
m="./Evolutionsgleichung.txt"
set terminal x11 0
set nokey
set grid
set title 'Die Evolutionsgleichung'
plot m
