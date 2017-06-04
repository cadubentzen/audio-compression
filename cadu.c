#include "cadu.h"


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

#if DEBUG
	char sample_name[30];
	FILE *testFile;
#endif

	printf( "Opening AIFF input file... " );
	ref = AIFF_OpenFile( aiff_input_path, F_RDONLY );
	if ( !ref )
	{
		fprintf( stderr, "Error opening AIFF input file.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}
	printf( "DONE\n" );

	printf( "Getting Audio Format info... " );
	if ( AIFF_GetAudioFormat( ref, nSamples, &channels_int,
							  samplingRate, &bitsPerSample_int,
							  &segmentSize ) < 1 )
	{
		fprintf( stderr, "Error getting audio format info.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}
	printf( "DONE\n\n" );

	*channels = (uint16_t) channels_int;
	*bitsPerSample = (uint16_t) bitsPerSample_int;

	printf( "aiff_input_path = %s\n", aiff_input_path );
	printf("nSamples = %" PRIu64 "\n", *nSamples);
	printf( "channels = %" PRIu16 "\n", *channels );
	printf( "samplingRate = %.2f\n", *samplingRate );
	printf( "bitsPerSample = %" PRIu16 "\n\n", *bitsPerSample );

	printf( "Reading samples from input audio file... " );
	samples_total = (float*) malloc( sizeof(float) * (*nSamples) * (*channels) );
	if( AIFF_ReadSamplesFloat( ref, samples_total, (*nSamples) * (*channels) ) < 1 )
	{
		fprintf( stderr, "Error getting audio samples.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}
	printf( "DONE\n" );

	printf( "Transfering channels to different arrays... " );
	*samples = (float**) malloc( sizeof(float*) * (*channels) );
	for ( j = 0; j < *channels; ++j )
		(*samples)[j] = (float*) fftwf_malloc( sizeof(float) * (*nSamples) );

	for ( i = 0; i < *nSamples; ++i )
		for ( j = 0; j < *channels; ++j )
		{
			(*samples)[j][i] = samples_total[ i*(*channels) + j ];
			// printf("%f\n", samples[j][i]);
		}

	printf( "DONE\n" );
	printf("\n");

#if DEBUG
	printf("\
############################\n\
  Testing read_aiff_file()\n\
############################\n\n"
	);

	for ( j = 0; j < *channels; ++j )
	{	
		sprintf( sample_name, "test_sample_channel_%" PRIu16 ".dat", j );
		printf( "Writing channel %" PRIu16 " to file %s... ", j, sample_name );
		testFile = fopen( sample_name, "w" );
		if ( testFile == NULL )
		{
			fprintf( stderr, "Error opening %s\n", sample_name );
			exit(EXIT_FAILURE);
		}

		for ( i = 0; i < *nSamples; ++i )
			fprintf( testFile, "%f\n", (*samples)[j][i] );

		printf( "DONE\n" );
		fclose(testFile);
	}

	printf("\nTest read_aiff_file() OK\n\n");
#endif
	
	AIFF_CloseFile(ref);
	free(samples_total);
}


void calc_channels_fft ( float **samples, complex_int32 ***ffts,
						 uint64_t nSamples, uint16_t channels )
{
	fftwf_complex *out;
	fftwf_plan p;

	uint16_t i;
	uint64_t j;

#if DEBUG
	char sample_name[30];
	FILE *testFile;
#endif

	printf( "Calculating FFTs... " );
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
	printf( "DONE\n" );

#if DEBUG
	printf("\
###############################\n\
  Testing calc_channels_fft()\n\
###############################\n\n"
	);

	for ( i = 0; i < channels; ++i )
	{
		sprintf( sample_name, "test_fft_channel_%" PRIu16 ".dat", i );
		printf( "Writing FFT magnitude of channel %" PRIu16 " to %s... ", i, sample_name );
		testFile = fopen(sample_name, "w");
		if ( testFile == NULL )
		{
			fprintf(stderr, "Error opening file %s.\n", sample_name);
			exit(EXIT_FAILURE);
		}

		for ( j = 0; j < (nSamples/2 + 1); ++j )
		{
			fprintf( testFile, "%" PRIu32 "\n", 
					 (uint32_t) sqrt( (*ffts)[i][j].a * (*ffts)[i][j].a + 
					 				  (*ffts)[i][j].b * (*ffts)[i][j].b ) );
		}
		printf( "DONE\n" );

		fclose(testFile);
	}
#endif

	printf("\n");
}


void crop_ffts ( complex_int32 **ffts, uint64_t **num_terms, uint16_t percentage,
				 uint64_t nSamples, uint16_t channels )
{
	uint64_t sum, j, sum_percentage;

	uint16_t i;

	*num_terms = malloc( sizeof(uint64_t) * channels );

#if DEBUG
	printf("Computing number of FFT terms for preservation of %u%% original area...\n", percentage);
#else
	printf("Computing number of FFT terms for preservation of %u%% original area... ", percentage);
#endif
	for ( i = 0; i < channels; ++i )
	{
		sum = 0;
		for ( j = 0; j < (nSamples/2 + 1); ++j )
			sum += (uint32_t) sqrt( ffts[i][j].a * ffts[i][j].a + 
					 				ffts[i][j].b * ffts[i][j].b );

		sum_percentage = (sum/100) * percentage;

	#if DEBUG
		printf("sum for channel %" PRIu16 " = %" PRIu64 "\n", i, sum);
		printf("sum_percentage for channel %" PRIu16 " = %" PRIu64 "\n", i, sum_percentage);
	#endif

		sum = 0;
		for ( j = 0; j < (nSamples/2 + 1); ++j )
		{
			sum += (uint32_t) sqrt( ffts[i][j].a * ffts[i][j].a + 
					 				ffts[i][j].b * ffts[i][j].b );
			if ( sum > sum_percentage  )
			{
			#if DEBUG
				printf("sum that reaches sum_percentage = %" PRIu64 "\n", sum);
			#endif
				(*num_terms)[i] = j;
				break;
			}	
		}
	}
#if !DEBUG
	printf( "DONE\n" );
#endif
	printf("\n");

#if DEBUG
	printf("\
###############################\n\
     Testing crop_ffts()\n\
###############################\n\n"
	);	

	for ( i = 0; i < channels; ++i )
	{
		printf( "Number of terms for channel %" PRIu16 ": %" PRIu64 " / %" PRIu64 " (%.2f%% of original terms)\n",
				i, (*num_terms)[i], (nSamples/2 + 1), 100*((double)(*num_terms)[i])/(nSamples/2 + 1) );
	}

	printf("\nTest crop_ffts() OK\n\n");
#endif
}

void discretize_ffts_int16 ( complex_int32 **ffts, complex_int16 ***ffts16, int32_t **max_terms,
							 uint64_t *num_terms, uint64_t nSamples, uint16_t channels )
{
	uint16_t i;
	uint64_t j;
	int32_t max_aux, aux;

#if DEBUG
	char sample_name[40];
	FILE *testFile;
#endif

	printf("Discretizing FFTs to 8 bits... ");
	// Allocate memory
	*ffts16 = (complex_int16**) malloc( sizeof(complex_int16*) *  channels );
	*max_terms = (int32_t*) malloc( sizeof(int32_t) * channels );
	for ( i = 0; i < channels; ++i )
	{
		(*ffts16)[i] = (complex_int16*) malloc( sizeof(complex_int16) * num_terms[i] );
		max_aux = 0;
		for ( j = 0; j < num_terms[i]; ++j)
		{
			aux = magnitude_complex_int32( ffts[i][j] );
			if ( aux > max_aux ) max_aux = aux;
		}
		(*max_terms)[i] = max_aux; 
	}

	// Populate ffts8
	for ( i = 0; i < channels; ++i )
	{
		for ( j = 0; j < num_terms[i]; ++j)
		{
			(*ffts16)[i][j].a = (int16_t) ( INT16_MAX * ( ( (float) ffts[i][j].a ) / (*max_terms)[i] ) );
			(*ffts16)[i][j].b = (int16_t) ( INT16_MAX * ( ( (float) ffts[i][j].b ) / (*max_terms)[i] ) );
		}
	}
	printf("DONE\n");

#if DEBUG
	printf("\
###############################\n\
  Testing discretize_ffts()\n\
###############################\n\n"
	);

	for ( i = 0; i < channels; ++i )
	{
		sprintf( sample_name, "test_fft_discrete_channel_%" PRIu16 ".dat", i );
		printf( "Writing FFT magnitude of channel %" PRIu16 " to %s... ", i, sample_name );
		testFile = fopen(sample_name, "w");
		if ( testFile == NULL )
		{
			fprintf(stderr, "Error opening file %s.\n", sample_name);
			exit(EXIT_FAILURE);
		}

		for ( j = 0; j < (nSamples/2 + 1); ++j )
		{
			if(j < num_terms[i])
			{
				fprintf( testFile, "%" PRIu16 "\n", 
					 (uint32_t) sqrt( (*ffts16)[i][j].a * (*ffts16)[i][j].a + 
					 				  (*ffts16)[i][j].b * (*ffts16)[i][j].b ) );
			}
			else {
				fprintf( testFile, "0\n");
			}
		}
		printf( "DONE\n" );

		fclose(testFile);
	}
#endif
	printf("\n");
}

void discretize_ffts_int8 ( complex_int32 **ffts, complex_int8 ***ffts8, int32_t **max_terms,
							uint64_t *num_terms, uint64_t nSamples, uint16_t channels )
{
	uint16_t i;
	uint64_t j;
	int32_t max_aux, aux;

#if DEBUG
	char sample_name[40];
	FILE *testFile;
#endif

	printf("Discretizing FFTs to 8 bits... ");
	// Allocate memory
	*ffts8 = (complex_int8**) malloc( sizeof(complex_int8*) *  channels );
	*max_terms = (int32_t*) malloc( sizeof(int32_t) * channels );
	for ( i = 0; i < channels; ++i )
	{
		(*ffts8)[i] = (complex_int8*) malloc( sizeof(complex_int8) * num_terms[i] );
		max_aux = 0;
		for ( j = 0; j < num_terms[i]; ++j)
		{
			aux = magnitude_complex_int32( ffts[i][j] );
			if ( aux > max_aux ) max_aux = aux;
		}
		(*max_terms)[i] = max_aux; 
	}

	// Populate ffts8
	for ( i = 0; i < channels; ++i )
	{
		for ( j = 0; j < num_terms[i]; ++j)
		{
			(*ffts8)[i][j].a = (int8_t) ( INT8_MAX * ( ( (float) ffts[i][j].a ) / (*max_terms)[i] ) );
			(*ffts8)[i][j].b = (int8_t) ( INT8_MAX * ( ( (float) ffts[i][j].b ) / (*max_terms)[i] ) );
		}
	}
	printf("DONE\n");

#if DEBUG
	printf("\
###############################\n\
  Testing discretize_ffts()\n\
###############################\n\n"
	);

	for ( i = 0; i < channels; ++i )
	{
		sprintf( sample_name, "test_fft_discrete_channel_%" PRIu16 ".dat", i );
		printf( "Writing FFT magnitude of channel %" PRIu16 " to %s... ", i, sample_name );
		testFile = fopen(sample_name, "w");
		if ( testFile == NULL )
		{
			fprintf(stderr, "Error opening file %s.\n", sample_name);
			exit(EXIT_FAILURE);
		}

		for ( j = 0; j < (nSamples/2 + 1); ++j )
		{
			if(j < num_terms[i])
			{
				fprintf( testFile, "%" PRIu32 "\n", 
					 (uint32_t) sqrt( (*ffts8)[i][j].a * (*ffts8)[i][j].a + 
					 				  (*ffts8)[i][j].b * (*ffts8)[i][j].b ) );
			}
			else {
				fprintf( testFile, "0\n");
			}
		}
		printf( "DONE\n" );

		fclose(testFile);
	}
#endif
	printf("\n");
}


void write_cadu_file_int32 ( char *cadu_output_path, complex_int32 **ffts, uint64_t *num_terms, uint64_t nSamples, 
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

	printf("Writing header to CADU file... ");
	fwrite( &nSamples,      sizeof(uint64_t),   1, caduFile );
	fwrite( &samplingRate,  sizeof(double),     1, caduFile );
	fwrite( &bitsPerSample, sizeof(uint16_t),   1, caduFile );
	fwrite( &channels,      sizeof(uint16_t),   1, caduFile );
	printf( "DONE\n" );

	printf("Writing body to CADU file... ");
	for ( i = 0; i < channels; ++i )
	{
		fwrite( num_terms+i, sizeof(uint64_t),      1,            caduFile );
		fwrite( ffts[i],     sizeof(complex_int32), num_terms[i], caduFile );
	}
	printf( "DONE\n" );

	fclose(caduFile);
	printf("\n");
}

void write_cadu_file_int16 ( char *cadu_output_path, complex_int16 **ffts16, uint64_t *num_terms, int32_t *max_terms,
							 uint64_t nSamples, uint16_t channels, double samplingRate, uint16_t bitsPerSample )
{
	uint16_t i;

	FILE *caduFile;

	caduFile = fopen( cadu_output_path, "wb" );
	if ( caduFile == NULL )
	{
		fprintf( stderr, "Error opening CADU output file %s.\n", cadu_output_path );
		exit(EXIT_FAILURE);
	}

	printf("Writing header to CADU file... ");
	fwrite( &nSamples,      sizeof(uint64_t),   1, caduFile );
	fwrite( &samplingRate,  sizeof(double),     1, caduFile );
	fwrite( &bitsPerSample, sizeof(uint16_t),   1, caduFile );
	fwrite( &channels,      sizeof(uint16_t),   1, caduFile );
	printf( "DONE\n" );

	printf("Writing body to CADU file... ");
	for ( i = 0; i < channels; ++i )
	{
		fwrite( num_terms+i, sizeof(uint64_t),      1,            caduFile );
		fwrite( max_terms+i, sizeof(int32_t),       1,            caduFile );
		fwrite( ffts16[i],   sizeof(complex_int16),  num_terms[i], caduFile );
	}
	printf( "DONE\n" );

	fclose(caduFile);
	printf("\n");
}

void write_cadu_file_int8 ( char *cadu_output_path, complex_int8 **ffts8, uint64_t *num_terms, int32_t *max_terms,
							uint64_t nSamples, uint16_t channels, double samplingRate, uint16_t bitsPerSample )
{
	uint16_t i;

	FILE *caduFile;

	caduFile = fopen( cadu_output_path, "wb" );
	if ( caduFile == NULL )
	{
		fprintf( stderr, "Error opening CADU output file %s.\n", cadu_output_path );
		exit(EXIT_FAILURE);
	}

	printf("Writing header to CADU file... ");
	fwrite( &nSamples,      sizeof(uint64_t),   1, caduFile );
	fwrite( &samplingRate,  sizeof(double),     1, caduFile );
	fwrite( &bitsPerSample, sizeof(uint16_t),   1, caduFile );
	fwrite( &channels,      sizeof(uint16_t),   1, caduFile );
	printf( "DONE\n" );

	printf("Writing body to CADU file... ");
	for ( i = 0; i < channels; ++i )
	{
		fwrite( num_terms+i, sizeof(uint64_t),      1,            caduFile );
		fwrite( max_terms+i, sizeof(int32_t),       1,            caduFile );
		fwrite( ffts8[i],    sizeof(complex_int8),  num_terms[i], caduFile );
	}
	printf( "DONE\n" );

	fclose(caduFile);
	printf("\n");
}


void read_cadu_file_int32 ( char *cadu_input_path, complex_int32 ***ffts, uint64_t *nSamples, 
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

	printf("Reading header of CADU file... ");
	fread( nSamples,      sizeof(uint64_t),   1, caduFile );
	fread( samplingRate,  sizeof(double),     1, caduFile );
	fread( bitsPerSample, sizeof(uint16_t),   1, caduFile );
	fread( channels,      sizeof(uint16_t),   1, caduFile );
	printf( "DONE\n" );

#if DEBUG
	printf( "\ncadu_input_path = %s\n", cadu_input_path );
	printf( "nSamples = %llu\n", *nSamples );
	printf( "samplingRate = %f\n", *samplingRate );
	printf( "bitsPerSample = %u\n", *bitsPerSample );
	printf( "channels = %u\n", *channels );
#endif

	printf("Reading body of CADU file... ");
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
	printf( "DONE\n" );
	printf( "\n" );

	fclose(caduFile);
}

void read_cadu_file_int16 ( char *cadu_input_path, complex_int32 ***ffts, uint64_t *nSamples, 
					  		uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample )
{
	FILE *caduFile;
	complex_int16 **ffts16;
	int32_t max_term;

#if DEBUG
	char fft_filename[40];
	FILE *fftFile;
#endif

	uint64_t num_terms, j;

	uint16_t i;

	caduFile = fopen( cadu_input_path, "rb" );
	if ( caduFile == NULL )
	{
		fprintf( stderr, "Error opening CADU input file %s.\n", cadu_input_path );
		exit(EXIT_FAILURE);
	}

	printf("Reading header of CADU file... ");
	fread( nSamples,      sizeof(uint64_t),   1, caduFile );
	fread( samplingRate,  sizeof(double),     1, caduFile );
	fread( bitsPerSample, sizeof(uint16_t),   1, caduFile );
	fread( channels,      sizeof(uint16_t),   1, caduFile );
	printf( "DONE\n" );

	ffts16 = (complex_int16**) malloc( sizeof(complex_int16*) * (*channels) );

#if DEBUG
	printf( "\ncadu_input_path = %s\n", cadu_input_path );
	printf( "nSamples = %llu\n", *nSamples );
	printf( "samplingRate = %f\n", *samplingRate );
	printf( "bitsPerSample = %u\n", *bitsPerSample );
	printf( "channels = %u\n", *channels );

	printf("Reading body of CADU file...\n");
#else
	printf("Reading body of CADU file... ");
#endif
	
	*ffts = (complex_int32**) malloc( sizeof(complex_int32*) * (*channels) );
	for ( i = 0; i < *channels; ++i )
	{
		(*ffts)[i] = (complex_int32*) malloc( sizeof(complex_int32) * ((*nSamples)/2 + 1) );
		memset( (*ffts)[i], 0, sizeof(complex_int32) * ((*nSamples)/2 + 1) );

		fread( &num_terms, sizeof(uint64_t), 1, caduFile );
		ffts16[i] = (complex_int16*) malloc( sizeof(complex_int16) * num_terms );

		fread( &max_term, sizeof(int32_t), 1, caduFile );

		fread( ffts16[i], sizeof(complex_int16), num_terms, caduFile );

		for ( j = 0; j < num_terms; ++j )
		{
			(*ffts)[i][j].a = (int32_t) ( ffts16[i][j].a * ( ( (float) max_term ) / INT16_MAX ) );
			(*ffts)[i][j].b = (int32_t) ( ffts16[i][j].b * ( ( (float) max_term ) / INT16_MAX ) );
		}

		free( ffts16[i] );


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

	printf( "DONE\n" );
	printf( "\n" );

	free(ffts16);

	fclose(caduFile);
}

void read_cadu_file_int8 ( char *cadu_input_path, complex_int32 ***ffts, uint64_t *nSamples, 
					  	   uint16_t *channels, double *samplingRate, uint16_t *bitsPerSample )
{
	FILE *caduFile;
	complex_int8 **ffts8;
	int32_t max_term;

#if DEBUG
	char fft_filename[40];
	FILE *fftFile;
#endif

	uint64_t num_terms, j;

	uint16_t i;

	caduFile = fopen( cadu_input_path, "rb" );
	if ( caduFile == NULL )
	{
		fprintf( stderr, "Error opening CADU input file %s.\n", cadu_input_path );
		exit(EXIT_FAILURE);
	}

	printf("Reading header of CADU file... ");
	fread( nSamples,      sizeof(uint64_t),   1, caduFile );
	fread( samplingRate,  sizeof(double),     1, caduFile );
	fread( bitsPerSample, sizeof(uint16_t),   1, caduFile );
	fread( channels,      sizeof(uint16_t),   1, caduFile );
	printf( "DONE\n" );

	ffts8 = (complex_int8**) malloc( sizeof(complex_int8*) * (*channels) );

#if DEBUG
	printf( "\ncadu_input_path = %s\n", cadu_input_path );
	printf( "nSamples = %llu\n", *nSamples );
	printf( "samplingRate = %f\n", *samplingRate );
	printf( "bitsPerSample = %u\n", *bitsPerSample );
	printf( "channels = %u\n", *channels );

	printf("Reading body of CADU file...\n");
#else
	printf("Reading body of CADU file... ");
#endif
	
	*ffts = (complex_int32**) malloc( sizeof(complex_int32*) * (*channels) );
	for ( i = 0; i < *channels; ++i )
	{
		(*ffts)[i] = (complex_int32*) malloc( sizeof(complex_int32) * ((*nSamples)/2 + 1) );
		memset( (*ffts)[i], 0, sizeof(complex_int32) * ((*nSamples)/2 + 1) );

		fread( &num_terms, sizeof(uint64_t), 1, caduFile );
		ffts8[i] = (complex_int8*) malloc( sizeof(complex_int8) * num_terms );

		fread( &max_term, sizeof(int32_t), 1, caduFile );

		fread( ffts8[i], sizeof(complex_int8), num_terms, caduFile );

		for ( j = 0; j < num_terms; ++j )
		{
			(*ffts)[i][j].a = (int32_t) ( ffts8[i][j].a * ( ( (float) max_term ) / INT8_MAX ) );
			(*ffts)[i][j].b = (int32_t) ( ffts8[i][j].b * ( ( (float) max_term ) / INT8_MAX ) );
		}

		free( ffts8[i] );


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

	printf( "DONE\n" );
	printf( "\n" );

	free(ffts8);

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

	printf("Computing inverse FFTs...\n");
#else
	printf("Computing inverse FFTs... ");
#endif

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
	printf( "DONE\n" );

	printf("\n");
}


void convert_to_write_format ( float **samples, int32_t **samples_final, 
							  uint64_t nSamples, uint16_t channels )
{
	uint16_t i;
	uint64_t j;

	printf("Converting samples to writing format... ");
	*samples_final = (int32_t*) malloc( sizeof(int32_t) * nSamples * channels );

	for ( i = 0; i < channels; ++i )
		for ( j = 0; j < nSamples; ++j )
			(*samples_final)[ j*channels + i ] = (int32_t) (INT32_MAX * samples[i][j]);

	printf( "DONE\n" );
}


void write_aiff_file ( char *aiff_output_path, int32_t *samples_final, uint64_t nSamples,
					   uint16_t channels, double samplingRate, uint16_t bitsPerSample )
{
	AIFF_Ref ref;

	printf("Opening output AIFF file... ");
	ref = AIFF_OpenFile(aiff_output_path, F_WRONLY);
	if ( !ref )
	{
		fprintf( stderr, "Error opening AIFF output file.\n" );
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}
	printf( "DONE\n" );

	printf("Setting output AIFF format... ");
	if ( AIFF_SetAudioFormat( ref, (int) channels, samplingRate, (int) bitsPerSample ) < 1 )
	{
		fprintf( stderr, "Error setting AIFF output format.\n");
		AIFF_CloseFile(ref);
		exit(EXIT_FAILURE);
	}
	printf( "DONE\n" );

	printf("Writing samples to AIFF output file... ");
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
	printf( "DONE\n" );

	printf("\n");

	AIFF_CloseFile(ref);
}

int32_t magnitude_complex_int32 ( complex_int32 complex_number )
{
	int32_t magnitude = (int32_t) sqrt( complex_number.a * complex_number.a +
										  complex_number.b * complex_number.b );

	return magnitude;
}


void free_arrays ( float **samples, complex_int32 **ffts, complex_int16 **ffts16, complex_int8 **ffts8, 
				   int32_t *samples_final, uint64_t *num_terms, uint64_t nSamples, 
				   uint16_t channels )
{
	uint16_t i;

	printf("Freeing memory... ");
	
	for ( i = 0; i < channels; ++i )
	{
		free( samples[i] );
		free( ffts[i] );
		if(ffts16 != NULL) free( ffts16[i] ); 
		if(ffts8 != NULL) free( ffts8[i] ); 
	}

	free(samples);
	free(ffts);
	if(ffts16 != NULL) free( ffts16 ); 
	if(ffts8 != NULL) free( ffts8 ); 
	if(samples_final != NULL) free(samples_final);
	if(num_terms != NULL) free(num_terms);

	printf( "DONE\n" );
	printf("\n");
}


void calc_compression_rate ( char *aiff_input_path, char *cadu_output_path )
{
	FILE *aiffFile, *caduFile;
	long size_aiff, size_cadu;
	double rate;

	printf("Calculating compression rate... ");

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

	printf( "DONE\n" );
	printf("%s is %.2f%% of %s file size\n", cadu_output_path, rate, aiff_input_path);
	printf("\n");
}
