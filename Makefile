FC = gfortran
FCFLAGS = -fopenmp

input.o: input.f90
	$(FC) $(FCFLAGS) -c input.f90

assign.o: assign.f90
	$(FC) $(FCFLAGS) -c assign.f90

density.o: density.f90
	$(FC) $(FCFLAGS) -c density.f90

order_layer.o: order_layer.f90
	$(FC) $(FCFLAGS) -c order_layer.f90

fluctuation.o: fluctuation.f90
	$(FC) $(FCFLAGS) -c fluctuation.f90

surface_wrap.o: surface_wrap.f90
	$(FC) $(FCFLAGS) -c surface_wrap.f90

water_angle.o: water_angle.f90
	$(FC) $(FCFLAGS) -c water_angle.f90

Hbonds.o: Hbonds.f90
	$(FC) $(FCFLAGS) -c Hbonds.f90

vvcf.o: vvcf.f90
	$(FC) $(FCFLAGS) -c vvcf.f90

vdos.o: vdos.f90
	$(FC) $(FCFLAGS) -c vdos.f90

assign: input.o assign.o
	$(FC) $(FCFLAGS) input.o assign.o -o assign

density: input.o density.o
	$(FC) $(FCFLAGS) input.o density.o -o density

order_layer: input.o order_layer.o
	$(FC) $(FCFLAGS) input.o order_layer.o -o order_layer

fluctuation: input.o fluctuation.o
	$(FC) $(FCFLAGS) input.o fluctuation.o -o fluctuation

surface_wrap: input.o surface_wrap.o
	$(FC) $(FCFLAGS) input.o surface_wrap.o -o surface_wrap

water_angle: input.o water_angle.o
	$(FC) $(FCFLAGS) input.o water_angle.o -o water_angle

Hbonds: input.o Hbonds.o
	$(FC) $(FCFLAGS) input.o Hbonds.o -o Hbonds

vvcf: input.o vvcf.o
	$(FC) $(FCFLAGS) input.o vvcf.o -o vvcf

vdos: input.o vdos.o
	$(FC) $(FCFLAGS) input.o vdos.o -o vdos

all: input.o density.o order_layer.o fluctuation.o assign.o surface_wrap.o water_angle.o Hbonds.o vvcf.o vdos.o
	$(FC) $(FCFLAGS) input.o density.o -o density
	$(FC) $(FCFLAGS) input.o order_layer.o -o order_layer
	$(FC) $(FCFLAGS) input.o fluctuation.o -o fluctuation
	$(FC) $(FCFLAGS) input.o assign.o -o assign
	$(FC) $(FCFLAGS) input.o surface_wrap.o -o surface_wrap
	$(FC) $(FCFLAGS) input.o water_angle.o -o water_angle
	$(FC) $(FCFLAGS) input.o Hbonds.o -o Hbonds
	$(FC) $(FCFLAGS) input.o vvcf.o -o vvcf
	$(FC) $(FCFLAGS) input.o vdos.o -o vdos

clean:
	rm -f *.o 

realclean:
	rm -f *.o *.mod assign density order_layer fluctuation surface_wrap water_angle Hbonds vvcf vdos
