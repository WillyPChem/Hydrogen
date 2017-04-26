#!/usr/bin/gnuplot


set terminal png

## The name of the file you want your plot saved as... must end in .png
set output 'Dipole_Plot_120.png'

## The command to plot your datafile
## Syntax is plot followed by filename with your data in quotation marks
## To specify which columns of data you want, immediately follow
## filename with u c1:c2 where c1 is the first column you want to plot
## and c2 is the second column you want to plot
## In these data files, first column is time, second column is dipole moment
## third column is electric field... to plot dipole moment on y axis and time on## x -axis, type u 1:2
## Style of plot can be specified... 'w l' means 'with lines'
## 'lw 2' means 'linewidth 2' which is wider than 'lw 1'
## You can specify a label for the data by typing title followed by 
## the label you want in quotation marks
plot 'dipoleMoment_2_sigma_120.txt' u 1:2 w l lw 2 title 'Dipole Moment', \
'dipoleMoment_2_sigma_120.txt' u 1:3 w l lw 2 title 'Electric Field'

set output 'Absorption_sigma_90.png'
plot 'AbsorptionSpectrum_sigma_90.txt' w l lw 2 title 'Absorption of H'
