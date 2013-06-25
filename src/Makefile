OCTFILES=mshm_refine.oct mshm_dolfin_read.oct mshm_dolfin_write.oct

all: $(OCTFILES)

%.oct:  %.cc
	mkoctfile $< -ldolfin

clean:
	-rm -f *.o core octave-core *.oct *~
# DO NOT DELETE