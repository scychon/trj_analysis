set term png enhanced font ',15'
set xlabel 'Simulation Time (ns)'
set ylabel 'Probability'
set key font ',12'
set output 'dpdoo_from_pr_pump.png'
set title 'P_{R}'' from P_{R}'
set xrange[0:25]
set yrange[0:1]
p 'dpdoo_from_pr_hbond_avg.dat' u ($1/1000):7 w l lw 5 title 'PUMP'
set output 'dpdoo_from_pr_chem.png'
p 'dpdoo_from_pr_hbond_avg.dat' u ($1/1000):13 w l lw 5 title 'CHEM'
#    EOF
