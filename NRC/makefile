all: labelling

exo1.o: exo1.c 
	gcc -c exo1.c nrio.c nralloc.c nrarith.c -lm -w 

exo1: exo1.o 
	gcc -o exo1 exo1.o nrio.c nralloc.c nrarith.c -lm -w 

binarization.o: binarization.c 
	gcc -c binarization.c nrio.c nralloc.c nrarith.c -lm -w 

binarization: binarization.o 
	gcc -o binarization binarization.o nrio.c nralloc.c nrarith.c -lm -w 	

labelling.o: labelling.c 
	gcc -c labelling.c nrio.c nralloc.c nrarith.c -lm -w 

labelling: labelling.o 
	gcc -o labelling labelling.o nrio.c nralloc.c nrarith.c -lm -w 	

clean: 
	rm -f -rf *.o