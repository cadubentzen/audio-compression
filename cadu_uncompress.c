#include <stdio.h>
#include <stdint.h>
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

void read_cadu_file ( char *cadu_input_path, complex_int32 ***ffts, uint64_t *nSamples, 
					  uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample );

void inverse_fft ( complex_int32 **ffts, float ***samples, uint64_t nSamples, uint16_t channels );

void convert_to_write_format ( float **samples, int32_t **samples_final, uint64_t nSamples, uint16_t channels );

void write_aiff_file ( char *aiff_output_path, int32_t *samples_final, uint64_t nSamples,
					   uint16_t channels, double samplingRate, uint16_t bitsPerSample );

void free_arrays( complex_int32 **ffts, float **samples, int32_t *samples_final, uint64_t nSamples,
				  uint16_t channels );

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

	read_cadu_file( cadu_input_path, &ffts, &nSamples, &channels,
				    &samplingRate, &bitsPerSample );

	inverse_fft( ffts, &samples, nSamples, channels );

	convert_to_write_format( samples, &samples_final, nSamples, channels );

	write_aiff_file( aiff_output_path, samples_final, nSamples,
					 channels, samplingRate, bitsPerSample );

	free_arrays( ffts, samples, samples_final, nSamples, channels );

	printf("All done xD\n");

	return 0;
}

void read_cadu_file ( char *cadu_input_path, complex_int32 ***ffts, uint64_t *nSamples, 
					  uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample )
{
	FILE *caduFile;
#if DEBUG
	char fft_filename[40];
	FILE *fftFile;
	uint64_t j;
#endif

	uint64_t num_terms;

	uint16_t i;

	caduFile = fopen( cadu_input_path, "rb" );
	if ( caduFile == NULL )
	{
		fprintf( stderr, "Error opening CADU input file %s.\n", cadu_input_path );
		exit(EXIT_FAILURE);
	}

	printf("Reading header of CADU file...\n");
	fread( nSamples,      sizeof(uint64_t),   1, caduFile );
	fread( samplingRate,  sizeof(double),     1, caduFile );
	fread( bitsPerSample, sizeof(uint16_t),   1, caduFile );
	fread( channels,      sizeof(uint16_t),   1, caduFile );

#if DEBUG
	printf( "cadu_input_path = %s\n", cadu_input_path );
	printf( "nSamples = %llu\n", *nSamples );
	printf( "samplingRate = %f\n", *samplingRate );
	printf( "bitsPerSample = %u\n", *bitsPerSample );
	printf( "channels = %u\n", *channels );
#endif

	printf("Reading body of CADU file...\n");
	*ffts = (complex_int32**) malloc( sizeof(complex_int32*) * (*channels) );
	for ( i = 0; i < *channels; ++i )
	{
		(*ffts)[i] = (complex_int32*) malloc( sizeof(complex_int32) * ((*nSamples)/2 + 1) );
		memset( (*ffts)[i], 0, sizeof(complex_int32) * ((*nSamples)/2 + 1) );

		fread( &num_terms,  sizeof(uint64_t),              1, caduFile );
		fread( (*ffts)[i],  sizeof(complex_int32), num_terms, caduFile );

	#if DEBUG
		sprintf( fft_filename, "test_fft_uncompressed_ch_%u.dat", i);
		printf("Writing FFT to %s...\n", fft_filename);
		fftFile = fopen( fft_filename, "w" );
		for ( j = 0; j < ((*nSamples)/2 + 1); ++j )
		{
			fprintf( fftFile, "%u\n",
					 (uint32_t) sqrt( (*ffts)[i][j].a * (*ffts)[i][j].a + 
					 				  (*ffts)[i][j].b * (*ffts)[i][j].b ) );
		}
		fclose(fftFile);
	#endif

	}

	printf("\n");

	fclose(caduFile);
}

void inverse_fft ( complex_int32 **ffts, float ***samples, uint64_t nSamples, uint16_t channels )
{
	fftwf_complex *in;
	fftwf_plan p;

	uint16_t i;
	uint64_t j;

#if DEBUG
	char sample_filename[40];
	FILE *sampleFile;
#endif

	printf("Computing inverse FFTs...\n");
	*samples = (float**) malloc( sizeof(float*) * channels );
	for ( i = 0; i < channels; ++i )
	{
		in = (fftwf_complex*) fftwf_malloc( sizeof(fftwf_complex) * (nSamples/2 + 1) );

		for ( j = 0; j < (nSamples/2 + 1); ++j )
		{
			in[j][0] = (float) ffts[i][j].a;
			in[j][1] = (float) ffts[i][j].b;
		}

		(*samples)[i] = (float*) fftwf_malloc( sizeof(float) * nSamples );
		p = fftwf_plan_dft_c2r_1d( nSamples, in, (*samples)[i], FFTW_ESTIMATE );

		fftwf_execute(p);

		fftwf_destroy_plan(p);
		fftwf_free(in);

	#if DEBUG
		sprintf( sample_filename, "test_sample_uncompressed_ch_%u.dat", i );
		printf("Writing samples to %s\n", sample_filename);
		sampleFile = fopen( sample_filename, "w" );
	#endif

		for ( j = 0; j < nSamples; ++j ) 
		{
			(*samples)[i][j] /= nSamples;
		#if DEBUG
			fprintf( sampleFile, "%f\n", (*samples)[i][j] );
		#endif
		}

	#if DEBUG
		fclose(sampleFile);
	#endif

	}

	printf("\n");
}

void convert_to_write_format ( float **samples, int32_t **samples_final, 
							  uint64_t nSamples, uint16_t channels )
{
	uint16_t i;
	uint64_t j;

	printf("Converting samples to writing format...\n");
	*samples_final = (int32_t*) malloc( sizeof(int32_t) * nSamples * channels );

	for ( i = 0; i < channels; ++i )
		for ( j = 0; j < nSamples; ++j )
			(*samples_final)[ j*channels + i ] = (int32_t) (INT32_MAX * samples[i][j]);
}

void write_aiff_file ( char *aiff_output_path, int32_t *samples_final, uint64_t nSamples,
					   uint16_t channels, double samplingRate, uint16_t bitsPerSample )
{
	AIFF_Ref ref;

	printf("Opening output AIFF file...\n");
	ref = AIFF_OpenFile(aiff_output_path, F_WRONLY);
	if ( !ref )
	{
		fprintf( stderr, "Error opening AIFF output file.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}

	printf("Setting output AIFF format...\n");
	if ( AIFF_SetAudioFormat( ref, (int) channels, samplingRate, (int) bitsPerSample ) < 1 )
	{
		fprintf( stderr, "Error setting AIFF output format.\n");
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}

	printf("Writing samples to AIFF output file...\n");
	if ( AIFF_StartWritingSamples(ref) < 1 )
	{
		fprintf( stderr, "Error starting writing samples in output.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}

	if ( AIFF_WriteSamples32Bit(ref, samples_final, nSamples * channels) < 1 )
	{
		fprintf( stderr, "Error writing samples to output.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}

	if( AIFF_EndWritingSamples(ref) < 1 )
	{
		fprintf( stderr, "Error ending writing samples in output.\n" );
	}

	printf("\n");

	AIFF_CloseFile(ref);
}

void free_arrays( complex_int32 **ffts, float **samples, int32_t *samples_final, uint64_t nSamples,
				  uint16_t channels )
{
	uint16_t i;
	uint64_t j;

	for ( i = 0; i < channels; ++i )
	{
		free( ffts[i] );
		fftwf_free( samples[i] );
	}

	free(ffts);
	free(samples);
	free(samples_final);
}