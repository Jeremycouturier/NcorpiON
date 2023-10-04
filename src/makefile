FF = gcc
OPTF = -O2

all: ncorpion

structure.o: structure.c structure.h parameters.h
	$(FF) $(OPTF) -c structure.c
	
physics.o: physics.c physics.h structure.o parameters.h collision.h rk4.h
	$(FF) $(OPTF) -c physics.c
	
collision.o: collision.c collision.h physics.o structure.o parameters.h ffm.h
	$(FF) $(OPTF) -c collision.c
	
ffm.o: ffm.c ffm.h structure.o parameters.h
	$(FF) $(OPTF) -c ffm.c

rk4.o: rk4.c rk4.h physics.o structure.o collision.o parameters.h ffm.h
	$(FF) $(OPTF) -c rk4.c
	
ncorpion: main.c structure.o physics.o rk4.o collision.o ffm.o parameters.h
	$(FF) main.c $(OPTF) structure.o physics.o collision.o rk4.o ffm.o -lm -o ncorpion

clean: 
	\rm -f *.o *.log
