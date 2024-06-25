export openGL_bool = 0
#Make sure to match openGL_bool to the same value in parameters.h

ifeq ($(openGL_bool), 1)
	FF = mpicc
else
	FF = gcc
endif

OPTF = -O2 -Wall -Wextra -pedantic -g

all: ncorpion

structure.o: structure.c structure.h parameters.h
	$(FF) $(OPTF) -c structure.c
	
physics.o: physics.c physics.h structure.o parameters.h
	$(FF) $(OPTF) -c physics.c
	
collision.o: collision.c collision.h physics.o structure.o parameters.h
	$(FF) $(OPTF) -c collision.c
	
ffm.o: ffm.c ffm.h collision.o physics.o structure.o parameters.h
	$(FF) $(OPTF) -c ffm.c

rk4.o: rk4.c rk4.h ffm.o collision.o physics.o structure.o parameters.h
	$(FF) $(OPTF) -c rk4.c
	
display.o: display.c display.h rk4.o ffm.o collision.o physics.o structure.o parameters.h
	$(FF) $(OPTF) -c display.c
	
spring.o: spring.c spring.h display.o rk4.o ffm.o collision.o physics.o structure.o parameters.h
	$(FF) $(OPTF) -c spring.c
	
ncorpion: main.c spring.o display.o rk4.o ffm.o collision.o physics.o structure.o parameters.h
	$(FF) main.c $(OPTF) structure.o physics.o collision.o ffm.o rk4.o display.o spring.o -lm -o ncorpion

clean: 
	\rm -f *.o *.log
