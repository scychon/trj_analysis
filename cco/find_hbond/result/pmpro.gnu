set term png enhanced font ',15'
set xlabel 'Simulation Time (ns)'
set ylabel 'Probability'
set key font ',12'
set output 'pmpro_ppdoo_solvate_pump.png'
set title 'P_{M,R} from P_{R}"'
set xrange[0:50]
set yrange[0:1]
p 'pmp_ppdoo_solvate_hbond_avg.dat' u ($1/1000):7 w l lw 5 title 'PUMP'
set output 'pmpro_ppdoo_solvate_chem.png'
p 'pmp_ppdoo_solvate_hbond_avg.dat' u ($1/1000):13 w l lw 5 title 'CHEM'
#    EOF
