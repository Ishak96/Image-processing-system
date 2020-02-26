all: extraction movementDetection morphoMath_EDOF hough labelling exo1 binarization

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

hough.o: hough.c 
	gcc -c hough.c nrio.c nralloc.c nrarith.c -lm -w 

hough: hough.o 
	gcc -o hough hough.o nrio.c nralloc.c nrarith.c -lm -w 	

morphoMath_EDOF.o: morphoMath_EDOF.c 
	gcc -c morphoMath_EDOF.c nrio.c nralloc.c nrarith.c -lm -w 

morphoMath_EDOF: morphoMath_EDOF.o 
	gcc -o morphoMath_EDOF morphoMath_EDOF.o nrio.c nralloc.c nrarith.c -lm -w 	

movementDetection.o: movementDetection.c 
	gcc -c movementDetection.c nrio.c nralloc.c nrarith.c -lm -w 

movementDetection: movementDetection.o 
	gcc -o movementDetection movementDetection.o nrio.c nralloc.c nrarith.c -lm -w 	

extraction.o: extraction.c 
	gcc -c extraction.c nrio.c nralloc.c nrarith.c -lm -w 

extraction: extraction.o 
	gcc -o extraction extraction.o nrio.c nralloc.c nrarith.c -lm -w 	

clean: 
	rm -f -rf *.o