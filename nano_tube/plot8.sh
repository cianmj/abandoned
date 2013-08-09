#!/bin/sh
cat > temp1.gpl <<-EOI
        set term post eps enh colour solid "Helvetica" 12
	set pointsize 0.2
	set xlabel "Position (a.u.)"
	set ylabel "Energy"
	set autoscale
        set xrange [0.0:1000.0]
        set yrange [-0.03:0.03]

        set data style lines
        set style line 1 linetype 1 lw 3
        set style line 2 linetype 2 lw 3
        set style line 3 linetype 3 lw 3
        set style line 4 linetype 4 lw 3
        set style line 5 linetype 5 lw 3
        set style line 6 linetype 6 lw 3
        set style line 7 linetype 7 lw 3
        set style line 8 linetype 8 lw 3
        set style line 9 linetype 9 lw 3

	set output "wfn.eps"
        plot "fort.30" u 1:2 ti "Potential" ls 1, \
         "fort.30" u 1:3 ti "|1>" ls 2, \
         "fort.30" u 1:4 ti "|2>" ls 3, \
         "fort.30" u 1:5 ti "|3>" ls 4 , \
         "fort.30" u 1:6 ti "|4>" ls 5 , \
         "fort.30" u 1:7 ti "|5>" ls 6 , \
         "fort.30" u 1:8 ti "|6>" ls 7 , \
         "fort.30" u 1:9 ti "|7>" ls 8 , \
         "fort.30" u 1:10 ti "|8>" ls 9
EOI
cat > temp2.gpl << EOI
        set term post eps enh colour solid "Helvetica" 12
	set pointsize 0.2
	set xlabel "Position (a.u.)"
	set ylabel "Energy"
	set autoscale
        set xrange [0.0:1000.0]
        set yrange [0.0:0.006]

        set data style lines
        set style line 1 linetype 1 lw 3
        set style line 2 linetype 2 lw 3
        set style line 3 linetype 3 lw 3
        set style line 4 linetype 4 lw 3
        set style line 5 linetype 5 lw 3
        set style line 6 linetype 6 lw 3
        set style line 7 linetype 7 lw 3
        set style line 8 linetype 8 lw 3
        set style line 9 linetype 9 lw 3

	set output "nrg.eps"
        plot "fort.30" u 1:2 ti "Potential" ls 1, \
         "fort.40" u 1:2 ti "E1" ls 2, \
         "fort.40" u 1:3 ti "E2" ls 3, \
         "fort.40" u 1:4 ti "E3" ls 4, \
         "fort.40" u 1:5 ti "E4" ls 5, \
         "fort.40" u 1:6 ti "E5" ls 6, \
         "fort.40" u 1:7 ti "E6" ls 7, \
         "fort.40" u 1:8 ti "E7" ls 8, \
         "fort.40" u 1:9 ti "E8" ls 9
EOI
gnuplot temp1.gpl
gnuplot temp2.gpl
rm temp1.gpl temp2.gpl
echo Finished.
