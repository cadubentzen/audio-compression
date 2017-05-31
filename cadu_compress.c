#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "cadu.h"


int main ( int argc, char **argv )
{
	// Input parameters
	char *aiff_input_path, *cadu_output_path;
	uint16_t percentage;

	// AIFF file parameters
	uint64_t nSamples;
	uint16_t channels;
	double   samplingRate;
	uint16_t bitsPerSample;

	// Array of samples
	float **samples;

	// Array of FFTs
	complex_int32 **ffts;

	// Array with numbers of elements to keep in each FFT
	uint64_t *num_terms;

	if ( argc != 4 )
	{
		printf( "Usage: ./audio_compression <AIFF input file> <CADU output file> <Percentage of information to keep>\n" );
		printf( "Example: ./audio_compression song.aif song.cadu 90\n" );
		exit(EXIT_FAILURE);
	}

	// Parse input parameters
	aiff_input_path = argv[1];
	cadu_output_path = argv[2];
	percentage = (uint16_t) atoi(argv[3]);

	// Read <AIFF input file>
	read_aiff_file( aiff_input_path, &samples, &nSamples, 
					&channels, &samplingRate, &bitsPerSample );

	// Calculate FFT of channels
	calc_channels_fft( samples, &ffts, nSamples, channels );

	// Crop FFTs to keep <Percentage of information to keep>
	crop_ffts( ffts, &num_terms, percentage, nSamples, channels );


	// Write <CADU output file>
	write_cadu_file( cadu_output_path, ffts, num_terms, nSamples, channels,
	 				 samplingRate, bitsPerSample );

	free_arrays( samples, ffts, NULL, num_terms, nSamples, channels );

	calc_compression_rate( aiff_input_path, cadu_output_path );

	printf("#ALL #DONE xD \n\n");

	return 0;
}
