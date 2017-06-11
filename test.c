#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "cadu.h"

#define COMPRESS 16


int main ( int argc, char **argv )
{
	// Input parameters
	char *aiff_input_path, *cadu_output_path;
	uint16_t percentage, i;
	uint16_t percentages[4] = {70, 80, 90, 95};

	// AIFF file parameters
	uint64_t nSamples;
	uint16_t channels;
	double   samplingRate;
	uint16_t bitsPerSample;

	// Array of samples
	float **samples;

	// Array of FFTs
	complex_int32 **ffts;
	complex_int16  **ffts16 = NULL;
	complex_int8  **ffts8 = NULL;

	// Array with numbers of elements to keep in each FFT
	uint64_t *num_terms;
	int32_t  *max_terms;

	if ( argc != 3 )
	{
		printf( "Usage: ./test <AIFF input file> <CADU output file>\n" );
		printf( "Example: ./test song.aif song.cadu\n" );
		exit(EXIT_FAILURE);
	}

	// Parse input parameters
	aiff_input_path = argv[1];
	cadu_output_path = argv[2];

	// Crop FFTs to keep <Percentage of information to keep>
	for (i = 0; i < 4; ++i)
	{
		// Read <AIFF input file>
		read_aiff_file( aiff_input_path, &samples, &nSamples, 
						&channels, &samplingRate, &bitsPerSample );

		// Calculate FFT of channels
		calc_channels_fft( samples, &ffts, nSamples, channels );

		percentage = percentages[i];
		crop_ffts( ffts, &num_terms, percentage, nSamples, channels );

		#if (COMPRESS == 16)
			// Discretize FFTs to 16 bits
			discretize_ffts_int16( ffts, &ffts16, &max_terms, num_terms, nSamples, channels );
			write_cadu_file_int16( cadu_output_path, ffts16, num_terms, max_terms,
								   nSamples, channels, samplingRate, bitsPerSample );
		#elif (COMPRESS == 8)
			// Discretize FFTs to 8 bits
			discretize_ffts_int8( ffts, &ffts8, &max_terms, num_terms, nSamples, channels );
			// Write <CADU output file>
			write_cadu_file_int8( cadu_output_path, ffts8, num_terms, max_terms,
								  nSamples, channels, samplingRate, bitsPerSample );
		#else
			write_cadu_file_int32( cadu_output_path, ffts, num_terms, max_terms,
								   nSamples, channels, samplingRate, bitsPerSample );
		#endif	

		free_arrays( samples, ffts, ffts16, ffts8, NULL, num_terms, nSamples, channels );

		calc_compression_rate( aiff_input_path, cadu_output_path );
	}
	
	

	printf("#ALL #DONE xD \n\n");

	return 0;
}
