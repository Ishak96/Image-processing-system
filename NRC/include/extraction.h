#ifndef __EXTRACTION_H__
#define __EXTRACTION_H__

byte** temporalAverage(const char* dir_name, const char* label, const char* extention,
					   long* nrl, long* nrh, long* ncl, long* nch, int n_seq);

byte** medianFilter(const char* dir_name, const char* label, const char* extention,
					   long* nrl, long* nrh, long* ncl, long* nch, int n_seq);

#endif