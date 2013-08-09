#!/bin/sh
cat > temp1.gpl <<-EOI
        set term post eps enh colour solid "Helvetica" 12
	set pointsize 0.2
	set xlabel "Position (a.u.)"
	set ylabel "Energy"
	set autoscale
        set xrange [0.0:1000.0]
#        set yrange [-0.03:0.03]
        set yrange [0.0:0.004]

        set data style lines
        set style line 1 linetype 1 lw 3
        set style line 2 linetype 2 lw 3
        set style line 3 linetype 3 lw 3
        set style line 4 linetype 4 lw 3

	set output "pwell.eps"
        plot "fort.30" u 1:2 ti "Potential" ls 1, \
         "fort.40" u 1:2 ti "E1" ls 2, \
         "fort.40" u 1:3 ti "E2" ls 3, \
         "fort.40" u 1:4 ti "E3" ls 4
#
#         "fort.30" u 1:3 ti "|1>" ls 2, \
#         "fort.30" u 1:4 ti "|2>" ls 3, \
#         "fort.30" u 1:5 ti "|3>" ls 4 , \
#
EOI
gnuplot temp1.gpl
rm temp1.gpl
#convert gpop.png gpop.ps
echo Finished.
