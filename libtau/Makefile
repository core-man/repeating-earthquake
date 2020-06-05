MOD = iasp91
SRC = libtau.c tauint.c tauspl.c em$(MOD).c
OBJ = $(subst .c,.o,$(SRC))

all: libtau remodl setbrn 

remodl: remodl.o 
	$(CC) -o remodl remodl.o -L. -ltau -lm

setbrn: setbrn.o 
	$(CC) -o setbrn setbrn.o -L. -ltau -lm

libtau: $(OBJ)
	ar rcs libtau.a $(OBJ)

clean:
	rm -f $(OBJ) remodl.o setbrn.o remodl setbrn libtau.a

install:
	./remodl 
	./setbrn
	rm -f *.lis remodl.hed

