#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <ctype.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "utils.h"
#include "erosion_dilatation.h"
#include "extraction.h"
#include "labelling.h"
#include "etiquettage_caracterisation.h"
#include "track.h"

extern int CEIL = 30;

int** mask = NULL;

const char* extention = "ppm";
const char* dir = "../Images/Sequences/Lbox/ppm/";
const char* label = "lbox";

const char* ref_image_file = "../Images/results/morphoMath/extraxtion/lbox_medianFilter.pgm";

const char* extract_image_file = "../Images/results/object_tracking/e.ppm";
const char* result_dir_extract = "../Images/results/object_tracking/";

const char* result_dir = "../Images/results/object_tracking/Lbox/";

int maskH = 3;
int maskW = 3;

void fill_maskR() {
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 1;
	mask[1][0] = 1; mask[1][1] = 1; mask[1][2] = 1;
	mask[2][0] = 1; mask[2][1] = 1; mask[2][2] = 1;
}

void fill_maskD() {
	mask[0][0] = 1; mask[0][1] = 1; mask[0][2] = 1;
	mask[1][0] = 1; mask[1][1] = 1; mask[1][2] = 1;
	mask[2][0] = 1; mask[2][1] = 1; mask[2][2] = 1;
}

void cleanImage(byte** f, long nrl, long nrh, long ncl, long nch)
{	
	int n = 1;
	
	byte** out = bmatrix(nrl, nrh, ncl, nch);
	
	//debruitage
	fill_maskD();
	n_erosion(f, out, nrl, nrh, ncl, nch, mask, maskW, maskH, n);

	n_dilatation(out, f, nrl, nrh, ncl, nch, mask, maskW, maskH, n);

	//remplissage
	fill_maskR();

	n_dilatation(f, out, nrl, nrh, ncl, nch, mask, maskW, maskH, n);

	n_erosion(out, f, nrl, nrh, ncl, nch, mask, maskW, maskH, n);

	free_bmatrix(out, nrl, nrh, ncl, nch);
}

void best_parameter(Parameters* param, int size, Parameters* prev_param)
{
	int h_err = 10000;
	int dist_min = 10000;
	int winer = -1;
	for(int i = 0; i < size; i++) {
		
		int curr_h_err = 0;
		for(int j = 0; j < 255; j++) {
			if(prev_param->H[j] != param[i].H[j]) {
				curr_h_err++;
			}
		}

		int curr_dist = abs(prev_param->OX - param[i].OX) + abs(prev_param->OY - param[i].OY);
		if(curr_h_err < h_err && curr_dist < dist_min && curr_dist < 10) {
			winer = i;
			h_err = curr_h_err;
			dist_min = curr_dist;
		}
	}

	if(winer == -1)
		prev_param->area_label = 0;
	else
		*prev_param = param[winer];
}

void track_object(byte** I, int** E, byte** byteI1, 
				  int track_label, int initial_pic, 
				  long nrl, long nrh, long ncl, long nch)
{
	int n = strlen(dir) + strlen(label) + 3 + strlen(extention);
	int m = strlen(result_dir) + strlen(label) + 3 + strlen(extention);

	int n_seq = fileCount(dir, label);

	int step = 1;
	int displayNext = step;
	int percent = 0;
	int inc = 1;
	int sum = 1;
	
	int vignette_size = 30;

	int nb_label;

	const char file_name[n];
	const char file_name_result[m];

	Vignette vignette;
	Vignette init_vignette;
	Parameters prev_param;
	Parameters init_param;

	int** harris_i = imatrix(nrl, nrh, ncl, nch);
	rgb8** I1;

	int init = 1;
	while(sum <= n_seq) {
		for(int i = initial_pic; i <= n_seq && i >= 1; i+=inc) {

			sprintf(file_name, "%s%s%03d.%s", dir, label, i, extention);
			I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);

			//traitement
			rgb8_minus(I1, I, byteI1, nrl, nrh, ncl, nch);
			
			cleanImage(byteI1, nrl, nrh, ncl, nch);

			LOOKUP_TABLE_LABELLING(byteI1, E, &nb_label, nrl, nrh, ncl, nch);

			Parameters* param = compute_parameters(I1, E, nb_label,
								  				   nrl, nrh, ncl, nch);

			if(init){
				prev_param = param[track_label - 1];
				init_param = prev_param;

				harris(harris_i, byteI1, nrl, nrh, ncl, nch);
				init_system(&vignette, I1, harris_i, E, 
							track_label, vignette_size, nrl, nrh, ncl, nch);
				init_vignette = vignette;
				init = 0;
			}

			track_label = track(E, I1, &vignette, nrl, nrh, ncl, nch);
			best_parameter(param, nb_label, &prev_param);

			if(track_label != 0)
				extract_zone(I1, byteI1, vignette, 1, nrl, nrh, ncl, nch);
			if(prev_param.area_label != 0)
				extract_zone_p(I1, byteI1, prev_param, 1, nrl, nrh, ncl, nch);

			sprintf(file_name_result, "%s%s%03d.%s", result_dir, label, i, extention);
			SavePPM_rgb8matrix(I1, nrl, nrh, ncl, nch, file_name_result);

			free_rgb8matrix(I1, nrl, nrh, ncl, nch);
			for(int i = 0; i < nb_label; i++)
				free_parameters(param+i);

			percent = 100 * (sum / (float)n_seq);
			print_progress(percent, 100);

			sum++;
		}
		inc = -1;
		initial_pic = initial_pic - 1;
		vignette = init_vignette;
		prev_param = init_param;
	}

	//Free matrix
	free_bmatrix(byteI1, nrl, nrh, ncl, nch);
	free_imatrix(harris_i, nrl, nrh, ncl, nch);
	destroy_mask();
}

byte** compute_ref_image(int* nrl, int* nrh, int* ncl, int* nch)
{	
	byte** I;

	int n_seq = fileCount(dir, label);

	if(access(ref_image_file, F_OK) != -1) {
		I = LoadPGM_bmatrix(ref_image_file, nrl, nrh, ncl, nch);
	} else {
		I = medianFilter(dir, label, extention, nrl, nrh, ncl, nch, n_seq);
		SavePGM_bmatrix(I, nrl, nrh, ncl, nch, ref_image_file);
	}

	return I;
}

void destroy(void) {
  gtk_main_quit();
}

int main(int argc, char** argv) {
  	GtkWidget* window;
  	GtkWidget* image;

	int track_label = 1;
	int initial_pic = 127;
	int start_track = 0;
	int nb_label;

	char request[3];
	long nrl, nrh, ncl, nch;

	if(argc < 2) {
		fprintf(stderr, "Please start the program like this -> $%s [LABEL]\n", argv[0]);
		printf("LABEL = Lbox | Fomd\n");
		return -1;
	}

	if(!strcmp(argv[1], "Fomd")) {
		dir = "../Images/Sequences/Fomd/ppm/";
		label = "fomd";

		ref_image_file = "../Images/results/morphoMath/extraxtion/fomd_medianFilter.pgm";

		result_dir = "../Images/results/object_tracking/Fomd/";	
	}


	int n = strlen(dir) + strlen(label) + 3 + strlen(extention);
	int m = strlen(result_dir_extract) + 1 + strlen(extention);
	
	const char file_name[n];
	const char file_name_result[m];
	
	byte** I = compute_ref_image(&nrl, &nrh, &ncl, &nch);
	rgb8** I1;
	
	int** E = imatrix(nrl, nrh, ncl, nch);
	byte** byteI1 = bmatrix(nrl, nrh, ncl, nch);
	mask = imatrix(0, maskH, 0, maskW);
	
	do {
		gtk_init (&argc, &argv);
		window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
		
		printf("Please enter the picture number in the directory :\n");
		scanf("%d", &initial_pic);

		printf("Please enter the label of the object :\n");
		scanf("%d", &track_label);

		if(initial_pic > 0 && track_label > 0) {
			sprintf(file_name, "%s%s%03d.%s", dir, label, initial_pic, extention);
			I1 = LoadPPM_rgb8matrix(file_name, &nrl, &nrh, &ncl, &nch);

			rgb8_minus(I1, I, byteI1, nrl, nrh, ncl, nch);
			
			cleanImage(byteI1, nrl, nrh, ncl, nch);

			LOOKUP_TABLE_LABELLING(byteI1, E, &nb_label, nrl, nrh, ncl, nch);
			EXTRACTZONE_RGB8(I1, E, track_label, nrl, nrh, ncl, nch);
			
			sprintf(file_name_result, "%se.%s", result_dir_extract, extention);
			SavePPM_rgb8matrix(I1, nrl, nrh, ncl, nch, extract_image_file);
		}
		
		image  = gtk_image_new_from_file(extract_image_file);
  		
  		gtk_signal_connect(GTK_OBJECT (window), "destroy", GTK_SIGNAL_FUNC (destroy), NULL);
		
		gtk_container_add(GTK_CONTAINER (window), image);

		gtk_widget_show_all(window);

  		gtk_main();

		printf("Do you want to start tracking : \n");
		
		if(scanf("%3s", request)) {
			for(int i = 0; request[i]; i++) request[i] = tolower(request[i]);
			
			if(!strcmp(request, "oui") || !strcmp(request, "yes"))
				start_track = 1;
		}
	} while(!start_track);

	track_object(I, E, byteI1, track_label, initial_pic, nrl, nrh, ncl, nch);

	free_bmatrix(I, nrl, nrh, ncl, nch);
	free_rgb8matrix(I1, nrl, nrh, ncl, nch);
	free_imatrix(E, nrl, nrh, ncl, nch);
	free_imatrix(mask, 0, maskH, 0, maskW);

	return 0;
}