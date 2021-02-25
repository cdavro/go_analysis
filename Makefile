FC = gfortran
FCFLAGS = -fopenmp
#FC = ifort
#FCFLAGS = -qopenmp

#Objects
input.o: input.f90
	$(FC) $(FCFLAGS) -c input.f90

sb_go.o: sb_go.f90
	$(FC) $(FCFLAGS) -c sb_go.f90

assign.o: assign.f90
	$(FC) $(FCFLAGS) -c assign.f90

assign_ff.o: assign_ff.f90
	$(FC) $(FCFLAGS) -c assign_ff.f90

extract.o: extract.f90
	$(FC) $(FCFLAGS) -c extract.f90

surface_wrap.o: surface_wrap.f90
	$(FC) $(FCFLAGS) -c surface_wrap.f90

dist.o: dist.f90
	$(FC) $(FCFLAGS) -c dist.f90

react_event.o: react_event.f90
	$(FC) $(FCFLAGS) -c react_event.f90

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

proton_hop.o: proton_hop.f90
	$(FC) $(FCFLAGS) -c proton_hop.f90

hydroxyde_hop.o: hydroxyde_hop.f90
	$(FC) $(FCFLAGS) -c hydroxyde_hop.f90

#Individual programs
assign: input.o sb_go.o assign.o
	$(FC) $(FCFLAGS) input.o sb_go.o assign.o -o assign

assign_ff: input.o sb_go.o assign_ff.o
	$(FC) $(FCFLAGS) input.o sb_go.o assign_ff.o -o assign_ff

surface_wrap: input.o sb_go surface_wrap.o
	$(FC) $(FCFLAGS) input.o sb_go.o surface_wrap.o -o surface_wrap

extract: input.o sb_go.o extract.o
	$(FC) $(FCFLAGS) input.o sb_go.o extract.o -o extract

density: input.o sb_go.o density.o
	$(FC) $(FCFLAGS) input.o sb_go.o density.o -o density

dist: input.o sb_go.o dist.o
	$(FC) $(FCFLAGS) input.o sb_go.o dist.o -o dist

react_event: input.o sb_go.o react_event.o
	$(FC) $(FCFLAGS) input.o sb_go.o react_event.o -o react_event

order_layer: input.o sb_go.o order_layer.o
	$(FC) $(FCFLAGS) input.o order_layer.o -o order_layer

fluctuation: input.o sb_go.o fluctuation.o
	$(FC) $(FCFLAGS) input.o sb_go.o fluctuation.o -o fluctuation

water_angle: input.o sb_go.o water_angle.o
	$(FC) $(FCFLAGS) input.o sb_go.o water_angle.o -o water_angle

hbonds: input.o sb_go.o hbonds.o
	$(FC) $(FCFLAGS) input.o sb_go.o hbonds.o -o hbonds

vvcf: input.o sb_go.o vvcf.o
	$(FC) $(FCFLAGS) input.o sb_go.o vvcf.o -o vvcf

proton_hop: input.o sb_go.o proton_hop.o
	$(FC) $(FCFLAGS) input.o sb_go.o proton_hop.o -o proton_hop

hydroxyde_hop: input.o sb_go.o hydroxyde_hop.o
	$(FC) $(FCFLAGS) input.o sb_go.o hydroxyde_hop.o -o hydroxyde_hop

all: input.o sb_go.o density.o order_layer.o fluctuation.o assign.o assign_ff.o extract.o surface_wrap.o dist.o react_event.o water_angle.o hbonds.o vvcf.o proton_hop.o hydroxyde_hop.o
	$(FC) $(FCFLAGS) input.o sb_go.o assign.o -o assign
	$(FC) $(FCFLAGS) input.o sb_go.o assign_ff.o -o assign_ff
	$(FC) $(FCFLAGS) input.o sb_go.o surface_wrap.o -o surface_wrap
	$(FC) $(FCFLAGS) input.o sb_go.o extract.o -o extract
	$(FC) $(FCFLAGS) input.o sb_go.o density.o -o density
	$(FC) $(FCFLAGS) input.o sb_go.o dist.o -o dist
	$(FC) $(FCFLAGS) input.o sb_go.o react_event.o -o react_event
	$(FC) $(FCFLAGS) input.o sb_go.o order_layer.o -o order_layer
	$(FC) $(FCFLAGS) input.o sb_go.o fluctuation.o -o fluctuation
	$(FC) $(FCFLAGS) input.o sb_go.o water_angle.o -o water_angle
	$(FC) $(FCFLAGS) input.o sb_go.o hbonds.o -o hbonds
	$(FC) $(FCFLAGS) input.o sb_go.o vvcf.o -o vvcf
	$(FC) $(FCFLAGS) input.o sb_go.o proton_hop.o -o proton_hop
	$(FC) $(FCFLAGS) input.o sb_go.o hydroxyde_hop.o -o hydroxyde_hop

clean:
	rm -f *.o *.mod

realclean:
	rm -f *.o *.mod assign assign_ff extract surface_wrap dist react_event density order_layer fluctuation water_angle hbonds vvcf proton_hop hydroxyde_hop
