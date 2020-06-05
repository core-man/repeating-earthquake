
CC = gcc
#CC = icc

SACLIB = /opt/SAC/lib
LIB = -lm

ALL : ccor clean


ccor : cross_correlation_sod.o
	$(CC) -o ccor_sac cross_correlation_sod.c -Llibtau -lfftw3 -ltau $(SACLIB)/libsac.a $(SACLIB)/sacio.a $(LIB)

clean :
	rm -rf *.o

