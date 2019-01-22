set logscale y
set ytics
set logscale y2
plot 'out.log' using 4:5 with lines axes x1y1, 'out.log' using 4:1 with lines axes x1y2
set grid
while (1) {
	replot
	pause 1
}
