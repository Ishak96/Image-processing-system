#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <gtk/gtk.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "corner_detection.h"
#include "utils.h"

void destroy(void) {
  gtk_main_quit();
}

int main(int argc, char** argv) {
  	GtkWidget* window;
  	GtkWidget* image;

	long nrh, nrl,
         nch, ncl;
	byte **I;

	gtk_init (&argc, &argv);

	I = LoadPGM_bmatrix("../Images/src/cubesx3.pgm", &nrl, &nrh, &ncl, &nch);

	int** harris_i = harris(I, nrl, nrh, ncl, nch);

	rgb8** red_corner = display_corner(I, harris_i, nrl, nrh, ncl, nch);

	SavePPM_rgb8matrix(red_corner, nrl, nrh, ncl, nch, "../Images/results/corner_detection/gradient_corner_detection_cubesx3.ppm");

	/* Free the byte sum matrix */
	free_bmatrix(I, nrl, nrh, ncl, nch);
	free_rgb8matrix(harris_i, nrl, nrh, ncl, nch);

	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	image  = gtk_image_new_from_file("../Images/results/corner_detection/gradient_corner_detection_cubesx3.ppm");

  	gtk_signal_connect(GTK_OBJECT (window), "destroy",
             GTK_SIGNAL_FUNC (destroy), NULL);

  	gtk_container_add(GTK_CONTAINER (window), image);

  	gtk_widget_show_all(window);

  	gtk_main();	
	
	return 0;
}