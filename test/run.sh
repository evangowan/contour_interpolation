#! /bin/bash


first_file=test_file2.gmt
second_file=test_file.gmt

xmin=-250
xmax=1250
ymin=2500
ymax=4500

R_option="-R${xmin}/${xmax}/${ymin}/${ymax}"

x_axis_length=15
y_axis_length=15

J_option="-JX${x_axis_length}/${y_axis_length}"

label_interval=250
sub_tick_interval=500

x_label="x"
y_label="y"

B_option="-Bf${label_interval}a${sub_tick_interval} -Bx+l${x_label} -By+l${y_label} -BWeSn"

plot=outlines_only.ps


psxy ${second_file}  ${J_option} ${R_option} ${B_option} -Wthick,blue -K -P  > ${plot}
psxy ${second_file} -R -J  -K -O -P -Sc0.1 -Gblack >> ${plot}
psxy ${first_file} -R -J -Wthick,red -K -O -P >> ${plot}
psxy ${first_file} -R -J -Sc0.1 -Gblack -K -O -P >> ${plot}

first_polygons=$( grep '>' ${first_file}  | wc -l )
second_polygons=$(  grep '>' ${second_file}  | wc -l )

cat << END_CAT > params.txt
${first_polygons} ${second_polygons}
${first_file}
${second_file}
END_CAT

cat params.txt


./../find_direction







triangulate direction_file.bin  -bi3d -boI -V > direction_file.tin

grid_spacing=20

tin_file_size=$( du -b direction_file.tin | awk '{print $1}' )

./../interpolate_direction ${grid_spacing} ${tin_file_size}

psxy direction_grid.txt -R -J -Sv0.15+e -Ggrey -K -O -P -Wthinnest >> ${plot}

psxy direction_file.txt -R -J -Sv0.18+e -Gblack -O -P -Wthinnest >> ${plot}


./../boundary_mask ${grid_spacing}

plot="mask_test1.ps"

psxy boundary_mask_1.txt  ${J_option} ${R_option} ${B_option} -Ss0.22 -Wthick -K -P  > ${plot}

psxy ${second_file}  -R -J -Wthick,blue -O -K -P  >> ${plot}
psxy ${first_file} -R -J -Wthick,red  -O -P >> ${plot}

plot="mask_test2.ps"

psxy boundary_mask_2.txt  ${J_option} ${R_option} ${B_option} -Ss0.22 -Wthick -K  -P  > ${plot}



psxy ${second_file}  -R -J -Wthick,blue -O -K -P  >> ${plot}
psxy ${first_file} -R -J -Wthick,red  -O -P >> ${plot}

plot="flowline_test.ps"

./../flowlines ${grid_spacing}

psxy ${second_file}  ${J_option} ${R_option} ${B_option} -Wthick,blue -K -P  > ${plot}

psxy ${first_file} -R -J -Wthick,red -K -O -P >> ${plot}

psxy direction_grid.txt -R -J -Sv0.15+e -Ggrey -K -O -P -Wthinnest >> ${plot}

psxy flowlines.txt  -R -J -Wthickest,green -K -O -P >> ${plot}

psxy flowlines.txt  -R -J -Wthin,black -K -O -P  >> ${plot}

#psxy fort.545 -R -J -Sc0.22 -Wthick,pink -O -K -P >> ${plot}


./../contour_creation 0.25

psxy final_contour.txt -R -J -Sd0.15 -Gpurple -O -K -P >> ${plot}

./../contour_creation 0.5

psxy final_contour.txt -R -J -Sd0.15 -Gbrown -O -K -P >> ${plot}


./../contour_creation 0.75

psxy final_contour.txt -R -J -Sd0.15 -Gmagenta -O -K -P >> ${plot}

psxy fort.546 -R -J -Sc0.22 -Wthick,yellow -O -K -P >> ${plot}

psxy fort.547 -R -J -Sc0.22 -Wthick,pink -O  -P >> ${plot}



./../contour_creation 0.1

surface final_contour.txt -Gsurface.nc -I${grid_spacing} ${R_option} -Ll-0.1 -Lu1.1

makecpt -Cjet -T-0.1/1.1/0.1  > shades.cpt

plot=contour_plot.ps
grdimage surface.nc  ${J_option} ${R_option} ${B_option} -K  -Cshades.cpt -P > $plot

psxy ${second_file} -J -R -Wthick,blue -K -O -P  >> ${plot}

psxy ${first_file} -R -J -Wthick,red  -K -O -P >> ${plot}

grdcontour surface.nc -J -R  -O -Cshades.cpt -W0.75p,black -A0.2+f9p,black+gwhite >> ${plot}





