FC = gfortran
FCFLAGS = -fopenmp

input.o: input.f90
	$(FC) $(FCFLAGS) -c input.f90

sb_go.o: sb_go.f90
	$(FC) $(FCFLAGS) -c sb_go.f90

assign.o: assign.f90
	$(FC) $(FCFLAGS) -c assign.f90

surface_wrap.o: surface_wrap.f90
	$(FC) $(FCFLAGS) -c surface_wrap.f90

density.o: density.f90
	$(FC) $(FCFLAGS) -c density.f90

order_layer.o: order_layer.f90
	$(FC) $(FCFLAGS) -c order_layer.f90

fluctuation.o: fluctuation.f90
	$(FC) $(FCFLAGS) -c fluctuation.f90

water_angle.o: water_angle.f90
	$(FC) $(FCFLAGS) -c water_angle.f90

hbonds.o: hbonds.f90
	$(FC) $(FCFLAGS) -c hbonds.f90

vvcf.o: vvcf.f90
	$(FC) $(FCFLAGS) -c vvcf.f90

vdos.o: vdos.f90
	$(FC) $(FCFLAGS) -c vdos.f90

assign: input.o assign.o
	$(FC) $(FCFLAGS) input.o assign.o -o assign

surface_wrap: input.o surface_wrap.o
	$(FC) $(FCFLAGS) input.o surface_wrap.o -o surface_wrap

density: input.o sb_go.o density.o
	$(FC) $(FCFLAGS) input.o sb_go.o density.o -o density

order_layer: input.o order_layer.o
	$(FC) $(FCFLAGS) input.o order_layer.o -o order_layer

fluctuation: input.o sb_go.o fluctuation.o
	$(FC) $(FCFLAGS) input.o sb_go.o fluctuation.o -o fluctuation

water_angle: input.o sb_go.o water_angle.o
	$(FC) $(FCFLAGS) input.o sb_go.o water_angle.o -o water_angle

hbonds: input.o sb_go.o hbonds.o
	$(FC) $(FCFLAGS) input.o sb_go.o hbonds.o -o hbonds

vvcf: input.o sb_go.o vvcf.o
	$(FC) $(FCFLAGS) input.o sb_go.o vvcf.o -o vvcf

vdos: input.o sb_go.o vdos.o
	$(FC) $(FCFLAGS) input.o sb_go.o vdos.o -o vdos

all: input.o sb_go.o density.o order_layer.o fluctuation.o assign.o surface_wrap.o water_angle.o hbonds.o vvcf.o vdos.o
	$(FC) $(FCFLAGS) input.o assign.o -o assign
	$(FC) $(FCFLAGS) input.o surface_wrap.o -o surface_wrap
	$(FC) $(FCFLAGS) input.o sb_go.o density.o -o density
	$(FC) $(FCFLAGS) input.o sb_go.o order_layer.o -o order_layer
	$(FC) $(FCFLAGS) input.o sb_go.o fluctuation.o -o fluctuation
	$(FC) $(FCFLAGS) input.o sb_go.o water_angle.o -o water_angle
	$(FC) $(FCFLAGS) input.o sb_go.o hbonds.o -o hbonds
	$(FC) $(FCFLAGS) input.o sb_go.o vvcf.o -o vvcf
	$(FC) $(FCFLAGS) input.o sb_go.o vdos.o -o vdos

clean:
	rm -f *.o 

realclean:
	rm -f *.o *.mod assign surface_wrap density order_layer fluctuation water_angle hbonds vvcf 
