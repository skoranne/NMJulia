set term postsc eps color
set output "M8_D9.eps"
set ylabel "Energy"
set xlabel "X"
set title "1D Maxwell's Equation M=8 Delta=1e-9 DT=DX"
plot './M.txt' using 1:2 with linespo title "Exact E(z)", './M.txt' using 1:3 with linespo title "Numerical E(z)", './M.txt' using 1:4 with linespo title "H(z)"
quit
