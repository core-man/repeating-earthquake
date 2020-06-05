
CC = gcc
#CC = icc

SACLIB = /opt/SAC/lib
LIB = -lm


ALL : cc clean

cc : cc_sac.o
	$(CC) -o ../bin/cc_sac cc_sac.c -lfftw3 $(SACLIB)/libsac.a $(SACLIB)/sacio.a $(LIB)

clean :
	rm -rf *.o

