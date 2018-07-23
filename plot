set key off

set term png
set output 'vector.png'
set size square
set autoscale fix
set contour base
set cntrparam level incremental -100, 0.5, 100
#unset surface

    set size square 
    unset title 
    set pm3d map 
    set xrange [0:16] 
    set yrange [0:16]
     set cbrange[0:1] 
   set palette rgbformulae 22,13,10 
    splot 'Resultant_velocty.txt' 
