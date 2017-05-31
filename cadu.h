#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DEBUG 0

#define LIBAIFF_NOCOMPAT 1 // do not use LibAiff 2 API compatibility
#include <libaiff/libaiff.h>
#include <fftw3.h>

// Complex number with integer coefficients: z = a + j*b
typedef struct {
	int32_t a;
	int32_t b;
} complex_int32;

// Functions to compress the AIFF file
void read_aiff_file ( char *aiff_input_path, float ***samples, uint64_t *nSamples, 
					  uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample );

void calc_channels_fft ( float **samples, complex_int32 ***ffts, 
						 uint64_t nSamples, uint16_t channels );

void crop_ffts ( complex_int32 **ffts, uint64_t **num_terms, uint16_t percentage,
				 uint64_t nSamples, uint16_t uchannels );

void write_cadu_file ( char *cadu_output_path, complex_int32 **ffts, uint64_t *num_terms, uint64_t nSamples, 
					   uint16_t channels, double samplingRate, uint16_t bitsPerSample );


// Function to uncompress the CADU file
void read_cadu_file ( char *cadu_input_path, complex_int32 ***ffts, uint64_t *nSamples, 
					  uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample );

void inverse_fft ( complex_int32 **ffts, float ***samples, uint64_t nSamples, uint16_t channels );

void convert_to_write_format ( float **samples, int32_t **samples_final, uint64_t nSamples, uint16_t channels );

void write_aiff_file ( char *aiff_output_path, int32_t *samples_final, uint64_t nSamples,
					   uint16_t channels, double samplingRate, uint16_t bitsPerSample );


// Auxiliar functions
void free_arrays ( float **samples, complex_int32 **ffts, int32_t *samples_final,
				   uint64_t *num_terms, uint64_t nSamples, uint16_t channels );

void calc_compression_rate ( char *aiff_input_path, char *cadu_output_path );
