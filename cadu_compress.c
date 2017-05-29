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

// Functions to compress the AIFF file
void read_aiff_file ( char *aiff_input_path, float ***samples, uint64_t *nSamples, 
					  uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample );

void calc_channels_fft ( float **samples, complex_int32 ***ffts, 
						 uint64_t nSamples, uint16_t channels );

void crop_ffts ( complex_int32 **ffts, uint64_t **num_terms, uint16_t percentage,
				 uint64_t nSamples, uint16_t uchannels );

void write_cadu_file ( char *cadu_output_path, complex_int32 **ffts, uint64_t *num_terms, uint64_t nSamples, 
					   uint16_t channels, double samplingRate, uint16_t bitsPerSample );

void free_arrays ( float **samples, complex_int32 **ffts, uint64_t *num_terms,
				   uint64_t nSamples, uint16_t channels );

void calc_compression_rate ( char *aiff_input_path, char *cadu_output_path );

// Testing functions
void test_read_aiff_file ( char *aiff_input_path, float ***samples, uint64_t *nSamples, 
					       uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample );

void test_calc_channels_fft ( float **samples, complex_int32 ***ffts,
							  uint64_t nSamples, uint16_t channels );

void test_crop_ffts ( complex_int32 **ffts, uint64_t **num_terms, uint16_t percentage,
					  uint64_t nSamples, uint16_t uchannels );

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

#if DEBUG
	// Testing place
	test_read_aiff_file( aiff_input_path, &samples, &nSamples, 
						 &channels, &samplingRate, &bitsPerSample );

	test_calc_channels_fft( samples, &ffts, nSamples, channels );

	test_crop_ffts( ffts, &num_terms, percentage, nSamples, channels );
#else
	// Read <AIFF input file>
	read_aiff_file( aiff_input_path, &samples, &nSamples, 
					&channels, &samplingRate, &bitsPerSample );

	// Calculate FFT of channels
	calc_channels_fft( samples, &ffts, nSamples, channels );

	// Crop FFTs to keep <Percentage of information to keep>
	crop_ffts( ffts, &num_terms, percentage, nSamples, channels );
#endif

	// Write <CADU output file>
	write_cadu_file( cadu_output_path, ffts, num_terms, nSamples, channels,
	 				 samplingRate, bitsPerSample );

	free_arrays( samples, ffts, num_terms, nSamples, channels );

	calc_compression_rate( aiff_input_path, cadu_output_path );

	printf("All done xD \n\n");

	return 0;
}

void read_aiff_file ( char *aiff_input_path, float ***samples, uint64_t *nSamples, 
					  uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample )
{
	AIFF_Ref ref;

	int channels_int;
	int bitsPerSample_int;
	int segmentSize;

	float *samples_total;

	uint64_t i;
	uint16_t j;

	printf( "Opening AIFF input file...\n" );
	ref = AIFF_OpenFile( aiff_input_path, F_RDONLY );
	if ( !ref )
	{
		fprintf( stderr, "Error opening AIFF input file.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}

	printf( "Getting Audio Format info...\n" );
	if ( AIFF_GetAudioFormat( ref, nSamples, &channels_int,
							  samplingRate, &bitsPerSample_int,
							  &segmentSize ) < 1 )
	{
		fprintf( stderr, "Error getting audio format info.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}
	*channels = (uint16_t) channels_int;
	*bitsPerSample = (uint16_t) bitsPerSample_int;

	printf( "aiff_input_path = %s\n", aiff_input_path );
	printf( "nSamples = %llu\n", *nSamples );
	printf( "channels = %u\n", *channels );
	printf( "samplingRate = %f\n", *samplingRate );
	printf( "bitsPerSample = %u\n", *bitsPerSample );

	printf( "Reading samples from input audio file...\n" );
	samples_total = (float*) malloc( sizeof(float) * (*nSamples) * (*channels) );
	if( AIFF_ReadSamplesFloat( ref, samples_total, (*nSamples) * (*channels) ) < 1 )
	{
		fprintf( stderr, "Error getting audio samples.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}

	printf( "Transfering channels to different arrays...\n" );
	*samples = (float**) malloc( sizeof(float*) * (*channels) );
	for ( j = 0; j < *channels; ++j )
		(*samples)[j] = (float*) fftwf_malloc( sizeof(float) * (*nSamples) );

	for ( i = 0; i < *nSamples; ++i )
		for ( j = 0; j < *channels; ++j )
		{
			(*samples)[j][i] = samples_total[ i*(*channels) + j ];
			// printf("%f\n", samples[j][i]);
		}

	printf("\n");
	
	AIFF_CloseFile(ref);
	free(samples_total);
}

void test_read_aiff_file ( char *aiff_input_path, float ***samples, uint64_t *nSamples, 
					  	   uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample )
{
	char sample_name[30];

	FILE *testFile;

	uint16_t i;
	uint64_t j;

	printf("\
############################\n\
  Testing read_aiff_file()\n\
############################\n\n"
	);

	read_aiff_file( aiff_input_path, samples, nSamples, 
					channels, samplingRate, bitsPerSample );

	for ( i = 0; i < *channels; ++i )
	{	
		sprintf( sample_name, "test_sample_channel_%u.dat", i );
		printf( "Writing channel %u to file %s...\n", i, sample_name );
		testFile = fopen( sample_name, "w" );
		if ( testFile == NULL )
		{
			fprintf( stderr, "Error opening %s\n", sample_name );
			exit(EXIT_FAILURE);
		}

		for ( j = 0; j < *nSamples; ++j )
			fprintf( testFile, "%f\n", (*samples)[i][j] );

		fclose(testFile);
	}

	printf("\nTest read_aiff_file() OK\n\n");
}

void calc_channels_fft ( float **samples, complex_int32 ***ffts,
						 uint64_t nSamples, uint16_t channels )
{
	fftwf_complex *out;
	fftwf_plan p;

	uint16_t i;
	uint64_t j;

	printf( "Calculating FFTs...\n" );
	out = (fftwf_complex*) fftwf_malloc( sizeof(fftwf_complex) * (nSamples/2 + 1) );
	*ffts = (complex_int32**) malloc( sizeof(complex_int32*) *  channels);
	for (i = 0; i < channels; ++i)
	{
		(*ffts)[i] = (complex_int32*) malloc( sizeof(complex_int32) * (nSamples/2 + 1) );
		p = fftwf_plan_dft_r2c_1d( nSamples, samples[i], out, FFTW_ESTIMATE );

		fftwf_execute(p);

		for ( j = 0; j < (nSamples/2 + 1); ++j )
		{
			(*ffts)[i][j].a = (int32_t) out[j][0];
			(*ffts)[i][j].b = (int32_t) out[j][1];
		}

		fftwf_destroy_plan(p);
	}

	printf("\n");
}

void test_calc_channels_fft ( float **samples, complex_int32 ***ffts,
							  uint64_t nSamples, uint16_t channels )
{
	char sample_name[30];

	FILE *testFile;

	uint16_t i;
	uint64_t j;

	printf("\
###############################\n\
  Testing calc_channels_fft()\n\
###############################\n\n"
	);

	calc_channels_fft( samples, ffts, nSamples, channels );

	for ( i = 0; i < channels; ++i )
	{
		sprintf( sample_name, "test_fft_channel_%u.dat", i );
		printf( "Writing FFT magnitude of channel %u to %s...\n", i, sample_name );
		testFile = fopen(sample_name, "w");
		if ( testFile == NULL )
		{
			fprintf(stderr, "Error opening file %s.\n", sample_name);
			exit(EXIT_FAILURE);
		}

		for ( j = 0; j < (nSamples/2 + 1); ++j )
		{
			fprintf( testFile, "%u\n", 
					 (uint32_t) sqrt( (*ffts)[i][j].a * (*ffts)[i][j].a + 
					 				  (*ffts)[i][j].b * (*ffts)[i][j].b ) );
		}

		fclose(testFile);
	}

	printf("\nTest calc_channels_fft() OK\n\n");
}

void crop_ffts ( complex_int32 **ffts, uint64_t **num_terms, uint16_t percentage,
				 uint64_t nSamples, uint16_t channels )
{
	uint64_t sum, j, sum_percentage;

	uint16_t i;

	*num_terms = malloc( sizeof(uint64_t) * channels );

	printf("Computing number of FFT terms for preservation of %u%% original area...\n", percentage);
	for ( i = 0; i < channels; ++i )
	{
		sum = 0;
		for ( j = 0; j < (nSamples/2 + 1); ++j )
			sum += (uint32_t) sqrt( ffts[i][j].a * ffts[i][j].a + 
					 				ffts[i][j].b * ffts[i][j].b );

		sum_percentage = (sum/100) * percentage;

	#if DEBUG
		printf("sum for channel %u = %llu\n", i, sum);
		printf("sum_percentage for channel %u = %llu\n", i, sum_percentage);
	#endif

		sum = 0;
		for ( j = 0; j < (nSamples/2 + 1); ++j )
		{
			sum += (uint32_t) sqrt( ffts[i][j].a * ffts[i][j].a + 
					 				ffts[i][j].b * ffts[i][j].b );
			if ( sum > sum_percentage  )
			{
			#if DEBUG
				printf("sum that reaches sum_percentage = %llu\n", sum);
			#endif
				(*num_terms)[i] = j;
				break;
			}	
		}
	}
	printf("\n");
}

void test_crop_ffts ( complex_int32 **ffts, uint64_t **num_terms, uint16_t percentage,
				 	  uint64_t nSamples, uint16_t channels )
{
	uint16_t i;

	printf("\
###############################\n\
     Testing crop_ffts()\n\
###############################\n\n"
	);	

	crop_ffts( ffts, num_terms, percentage, nSamples, channels );

	for ( i = 0; i < channels; ++i )
	{
		printf( "Number of terms for channel %u: %llu / %llu (%.2f%% of original terms)\n",
				i, (*num_terms)[i], (nSamples/2 + 1), 100*((double)(*num_terms)[i])/(nSamples/2 + 1) );
	}

	printf("\nTest crop_ffts() OK\n\n");
}

void write_cadu_file ( char *cadu_output_path, complex_int32 **ffts, uint64_t *num_terms, uint64_t nSamples, 
					   uint16_t channels, double samplingRate, uint16_t bitsPerSample )
{
	uint16_t i;

	FILE *caduFile;

	caduFile = fopen( cadu_output_path, "wb" );
	if ( caduFile == NULL )
	{
		fprintf( stderr, "Error opening CADU output file %s.\n", cadu_output_path );
		exit(EXIT_FAILURE);
	}

	printf("Writing header to CADU file...\n");
	fwrite( &nSamples,      sizeof(uint64_t),   1, caduFile );
	fwrite( &samplingRate,  sizeof(double),     1, caduFile );
	fwrite( &bitsPerSample, sizeof(uint16_t),   1, caduFile );
	fwrite( &channels,      sizeof(uint16_t),   1, caduFile );

	printf("Writing body to CADU file...\n");
	for ( i = 0; i < channels; ++i )
	{
		fwrite( num_terms+i, sizeof(uint64_t),      1,            caduFile );
		fwrite( ffts[i],     sizeof(complex_int32), num_terms[i], caduFile );
	}

	fclose(caduFile);
	printf("\n");
}

void free_arrays ( float **samples, complex_int32 **ffts, uint64_t *num_terms,
				   uint64_t nSamples, uint16_t channels )
{
	uint16_t i;

	printf("Freeing memory...\n");
	
	for ( i = 0; i < channels; ++i )
	{
		free( samples[i] );
		free( ffts[i] );
	}

	free(samples);
	free(ffts);
	free(num_terms);

	printf("\n");
}

void calc_compression_rate ( char *aiff_input_path, char *cadu_output_path )
{
	FILE *aiffFile, *caduFile;
	long size_aiff, size_cadu;
	double rate;

	printf("Calculating compression rate...\n");

	aiffFile = fopen(aiff_input_path, "r");
	if ( aiffFile == NULL )
	{
		fprintf( stderr, "Error opening AIFF input file %s.\n", aiff_input_path );
		exit(EXIT_FAILURE);
	}

	caduFile = fopen( cadu_output_path, "r" );
	if ( caduFile == NULL )
	{
		fprintf( stderr, "Error opening CADU output file %s.\n", cadu_output_path );
		exit(EXIT_FAILURE);
	}

	fseek(aiffFile, 0, SEEK_END);
	size_aiff = ftell(aiffFile);
	fclose(aiffFile);

	fseek(caduFile, 0, SEEK_END);
	size_cadu = ftell(caduFile);
	fclose(caduFile);

	rate = 100 * ((double)size_cadu) / size_aiff;

	printf("%s is %.2f%% of %s file size\n", cadu_output_path, rate, aiff_input_path);
	printf("\n");
}