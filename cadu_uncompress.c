#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "cadu.h"

#define COMPRESS 16


int main(int argc, char **argv)
{
	// Input parameters
	char *cadu_input_path, *aiff_output_path;

	// AIFF file parameters
	uint64_t nSamples;
	uint16_t channels;
	double   samplingRate;
	uint16_t bitsPerSample;

	// Array of FFTs
	complex_int32 **ffts;

	// Array of samples and unified samples
	float **samples;
	int32_t *samples_final;

	if ( argc != 3 )
	{
		printf( "Usage: ./audio_compression <CADU input file> <AIFF output file>\n" );
		printf( "Example: ./audio_compression song.cadu song_uncompressed.aiff\n" );
		exit(EXIT_FAILURE);
	}

	cadu_input_path = argv[1];
	aiff_output_path = argv[2];

#if (COMPRESS == 8)
	read_cadu_file_int8( cadu_input_path, &ffts, &nSamples, &channels,
				    	 &samplingRate, &bitsPerSample );
#elif (COMPRESS == 16)
	read_cadu_file_int16( cadu_input_path, &ffts, &nSamples, &channels,
				    	  &samplingRate, &bitsPerSample );
#else
	read_cadu_file_int32( cadu_input_path, &ffts, &nSamples, &channels,
				    	  &samplingRate, &bitsPerSample );
#endif

	inverse_fft( ffts, &samples, nSamples, channels );

	convert_to_write_format( samples, &samples_final, nSamples, channels );

	write_aiff_file( aiff_output_path, samples_final, nSamples,
					 channels, samplingRate, bitsPerSample );

	free_arrays( samples, ffts, NULL, NULL, samples_final, NULL, nSamples, channels );

	printf("#ALL #DONE xD\n");

	return 0;
}
