FC = gfortran

FCFLAGS = -g -fbacktrace -fbounds-check  -ffpe-trap='zero'
#FCFLAGS = -O2 

# if compiling with the Intel Fortran compiler, you need to add an extra flag
#FCFLAGS = -O2 -assume byterecl






#####################

# modules
point_in_polygon.o: point_in_polygon.f90
	$(FC) -o point_in_polygon.o -c  $(FCFLAGS) point_in_polygon.f90

direction_mod.o: direction_mod.f90
	$(FC) -o direction_mod.o -c  $(FCFLAGS) direction_mod.f90

read_polygons.o: read_polygons.f90
	$(FC) -o read_polygons.o -c  $(FCFLAGS) read_polygons.f90

read_minmax.o: read_minmax.f90
	$(FC) -o read_minmax.o -c  $(FCFLAGS) read_minmax.f90

direction_grid_mod.o: direction_grid_mod.f90
	$(FC) -o direction_grid_mod.o -c  $(FCFLAGS) direction_grid_mod.f90

boundary_mask_mod.o: boundary_mask_mod.f90
	$(FC) -o boundary_mask_mod.o -c  $(FCFLAGS) boundary_mask_mod.f90


spline.o: spline.f90
	$(FC) -o spline.o -c  $(FCFLAGS) spline.f90

crossover.o: crossover.f90
	$(FC) -o crossover.o -c $(FCFLAGS) crossover.f90

tin.o: tin.f90
	$(FC) -o tin.o -c  $(FCFLAGS) tin.f90 

# programs

grid_direction_test: grid_direction_test.f90
	$(FC) -o grid_direction_test  $(FCFLAGS) grid_direction_test.f90

direction_flow: direction_flow.f90
	$(FC) -o direction_flow  $(FCFLAGS) direction_flow.f90

find_direction: find_direction.f90 point_in_polygon.o direction_mod.o read_polygons.o
	$(FC) -o find_direction  $(FCFLAGS) find_direction.f90 point_in_polygon.o direction_mod.o read_polygons.o


shrinking: shrinking.f90 read_polygons.o direction_mod.o point_in_polygon.o
	$(FC) -o shrinking  $(FCFLAGS) shrinking.f90 read_polygons.o direction_mod.o point_in_polygon.o


interpolate_direction: interpolate_direction.f90 point_in_polygon.o direction_mod.o read_minmax.o direction_grid_mod.o 
	$(FC) -o interpolate_direction  $(FCFLAGS) interpolate_direction.f90 point_in_polygon.o direction_mod.o read_minmax.o direction_grid_mod.o

flowline_direction: flowline_direction.f90 point_in_polygon.o direction_mod.o read_minmax.o crossover.o tin.o boundary_mask_mod.o read_polygons.o
	$(FC) -o flowline_direction  $(FCFLAGS) flowline_direction.f90 point_in_polygon.o direction_mod.o read_minmax.o crossover.o tin.o boundary_mask_mod.o read_polygons.o

boundary_mask: boundary_mask.f90 read_polygons.o read_minmax.o boundary_mask_mod.o point_in_polygon.o
	$(FC) -o boundary_mask  $(FCFLAGS) boundary_mask.f90 read_polygons.o read_minmax.o boundary_mask_mod.o point_in_polygon.o



flowlines: flowlines.f90 read_polygons.o read_minmax.o boundary_mask_mod.o direction_grid_mod.o crossover.o
	$(FC) -o flowlines  $(FCFLAGS) flowlines.f90 read_polygons.o read_minmax.o boundary_mask_mod.o direction_grid_mod.o crossover.o

flowlines3: flowlines3.f90 read_polygons.o read_minmax.o boundary_mask_mod.o direction_grid_mod.o crossover.o
	$(FC) -o flowlines3  $(FCFLAGS) flowlines3.f90 read_polygons.o read_minmax.o boundary_mask_mod.o direction_grid_mod.o crossover.o

contour_creation: contour_creation.f90 spline.o read_polygons.o
	$(FC) -o contour_creation  $(FCFLAGS) contour_creation.f90 spline.o read_polygons.o

eliminate_outside: eliminate_outside.f90 point_in_polygon.o  read_polygons.o
	$(FC) -o eliminate_outside  $(FCFLAGS) eliminate_outside.f90 point_in_polygon.o  read_polygons.o


merge_polygons: merge_polygons.f90
	$(FC) -o merge_polygons  $(FCFLAGS) merge_polygons.f90
