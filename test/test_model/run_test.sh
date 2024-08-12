#!/bin/bash


set -e

#---------------------------------------------------------------------
function run_prog ()
{
    prog=$1

    echo "Compiling..."
    make $prog

    echo "Running program..."
    ./${prog} > ${prog}_eqcat.out

    echo "Plotting..."
    gnuplot <<EOF
set terminal pngcairo enhanced font 'Verdana,8'
set output '${prog}_timemag.png'
set xrange [0:10]
set xlabel "Time [yr]"
set ylabel "M_W"
plot '${prog}_eqcat.out' u 1:5 w points
EOF

    echo "Displaying"
    if [ -x /usr/bin/eog ]; then
        eog ${prog}_timemag.png
    elif [ -x /usr/bin/ristretto ]; then
        ristretto ${prog}_timemag.png
    else
        echo Output check plot written to $prog.png
    fi

}
#---------------------------------------------------------------------

run_prog test_model