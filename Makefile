FC = gfortran
FCFLAGS = -fopenmp

input.o: input.f90
	$(FC) $(FCFLAGS) -c input.f90

assign.o: assign.f90
	$(FC) $(FCFLAGS) -c assign.f90

density.o: density.f90
	$(FC) $(FCFLAGS) -c density.f90

surface_wrap.o: surface_wrap.f90
	$(FC) $(FCFLAGS) -c surface_wrap.f90

water_angle.o: water_angle.f90
	$(FC) $(FCFLAGS) -c water_angle.f90

assign: input.o assign.o
	$(FC) $(FCFLAGS) input.o assign.o -o assign

density: input.o density.o
	$(FC) $(FCFLAGS) input.o density.o -o density

surface_wrap: input.o surface_wrap.o
	$(FC) $(FCFLAGS) input.o surface_wrap.o -o surface_wrap

water_angle: input.o water_angle.o
	$(FC) $(FCFLAGS) input.o water_angle.o -o water_angle

all: input.o density.o assign.o surface_wrap.o water_angle.o
	$(FC) $(FCFLAGS) input.o density.o -o density
	$(FC) $(FCFLAGS) input.o assign.o -o assign
	$(FC) $(FCFLAGS) input.o surface_wrap.o -o surface_wrap
	$(FC) $(FCFLAGS) input.o water_angle.o -o water_angle

clean:
	rm -f *.o

realclean:
	rm -f *.o assign density surface_wrap water_angle