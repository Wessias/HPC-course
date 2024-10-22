#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <threads.h>
#include <string.h>

//Pre-defined values and stuff I want to be global
#define MAX_ITER 128 
#define TOLERANCE 0.001
#define TOLERANCE_SQUARED 0.000001
#define RADIUS 0.998001
#define ESCAPE_RADIUS 10000000000.0
#define MAX_RE_AND_IM 2.0
#define MIN_RE_AND_IM -2.0
int num_threads = 1;
int size = 100;
int degree = 1;
double complex roots[10] = {0}; 
double step_size = 0.0;
unsigned char **attractors;
unsigned char **convergences;


//------------------------------------------------------------------------------------------------------------------
//Pretty colors which were deliberately chosen to optimize the pleasure of your eyes
const char *rgb_values[12] = {
  "255 99 71 ",    // Tomato (Warm Red) Close to 0
  "100 149 237 ",  // Cornflower Blue (Cool Blue) Close to root: 1
  "255 215 0 ",    // Gold (Warm Yellow) 			 2
  "60 179 113 ",   // Medium Sea Green (Cool Green)		 3
  "255 105 180 ",  // Hot Pink (Warm Pink)			 4
  "75 0 130 ",     // Indigo (Deep Blue)			 5
  "255 160 122 ",  // Light Salmon (Warm Light Orange)	         6
  "173 216 230 ",  // Light Blue (Soft Cool Blue)		 7
  "240 230 140 ",  // Khaki (Warm Beige)			 8
  "50 205 50 ",    // Lime Green (Vibrant Green)		 9
  "238 130 238 ",  // Violet (Bright Purple)		10, I thought the degree was allowed to be 10 at first but when I eventually realized this was not the case I could not be bothered changing the indexing 
  "0 0 0 "   // Black despair                                    diverges
};

// Define greyscale values as strings with a space at the end
const char *greyscale[129] = {
  "0 0 0 ", "1 1 1 ", "3 3 3 ", "5 5 5 ", "7 7 7 ",
  "9 9 9 ", "11 11 11 ", "13 13 13 ", "15 15 15 ", "17 17 17 ",
  "19 19 19 ", "21 21 21 ", "23 23 23 ", "25 25 25 ", "27 27 27 ",
  "29 29 29 ", "31 31 31 ", "33 33 33 ", "35 35 35 ", "37 37 37 ",
  "39 39 39 ", "41 41 41 ", "43 43 43 ", "45 45 45 ", "47 47 47 ",
  "49 49 49 ", "51 51 51 ", "53 53 53 ", "55 55 55 ", "57 57 57 ",
  "59 59 59 ", "61 61 61 ", "63 63 63 ", "65 65 65 ", "67 67 67 ",
  "69 69 69 ", "71 71 71 ", "73 73 73 ", "75 75 75 ", "77 77 77 ",
  "79 79 79 ", "81 81 81 ", "83 83 83 ", "85 85 85 ", "87 87 87 ",
  "89 89 89 ", "91 91 91 ", "93 93 93 ", "95 95 95 ", "97 97 97 ",
  "99 99 99 ", "101 101 101 ", "103 103 103 ", "105 105 105 ", "107 107 107 ",
  "109 109 109 ", "111 111 111 ", "113 113 113 ", "115 115 115 ", "117 117 117 ",
  "119 119 119 ", "121 121 121 ", "123 123 123 ", "125 125 125 ", "127 127 127 ",
  "129 129 129 ", "131 131 131 ", "133 133 133 ", "135 135 135 ", "137 137 137 ",
  "139 139 139 ", "141 141 141 ", "143 143 143 ", "145 145 145 ", "147 147 147 ",
  "149 149 149 ", "151 151 151 ", "153 153 153 ", "155 155 155 ", "157 157 157 ",
  "159 159 159 ", "161 161 161 ", "163 163 163 ", "165 165 165 ", "167 167 167 ",
  "169 169 169 ", "171 171 171 ", "173 173 173 ", "175 175 175 ", "177 177 177 ",
  "179 179 179 ", "181 181 181 ", "183 183 183 ", "185 185 185 ", "187 187 187 ",
  "189 189 189 ", "191 191 191 ", "193 193 193 ", "195 195 195 ", "197 197 197 ",
  "199 199 199 ", "201 201 201 ", "203 203 203 ", "205 205 205 ", "207 207 207 ",
  "209 209 209 ", "211 211 211 ", "213 213 213 ", "215 215 215 ", "217 217 217 ",
  "219 219 219 ", "221 221 221 ", "223 223 223 ", "225 225 225 ", "227 227 227 ",
  "229 229 229 ", "231 231 231 ", "233 233 233 ", "235 235 235 ", "237 237 237 ",
  "239 239 239 ", "241 241 241 ", "243 243 243 ", "245 245 245 ", "247 247 247 ",
  "249 249 249 ", "251 251 251 ", "253 253 253 ", "255 255 255 "
};
//------------------------------------------------------------------------------------------------------------------

typedef struct {
  int val;
  char pad[60];
} int_padded;


typedef struct {
  int id;
  unsigned char *attractor; //These track attractor index and iteration count 
  unsigned char *convergence;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thread_data_t;


typedef struct{
  char filenameAttr[50];
  char filenameConv[50];
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} write_data_t;

//------------------------------------------------------------------------------------------------------------------


int write_ppm_files(void* arg) {
    write_data_t *data = (write_data_t *)arg;

    // Open attractor file
    FILE *fp_attractor = fopen(data->filenameAttr, "w");
    if (fp_attractor == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n", data->filenameAttr);
        return thrd_error;
    }

    // Open convergence file
    FILE *fp_convergence = fopen(data->filenameConv, "w");
    if (fp_convergence == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n", data->filenameConv);
        fclose(fp_attractor);
        return thrd_error;
    }

    // Write the PPM headers for both files
    char header[50];
    int header_size = snprintf(header, sizeof(header), "P3\n%d %d\n255\n", size, size);
    fwrite(header, sizeof(char), header_size, fp_attractor);
    fwrite(header, sizeof(char), header_size, fp_convergence);

    // Allocate buffer for a single row for both files
    char *buffer_attractor = malloc(size * 12); // max "255 255 255 " for each pixel (12 chars)
    char *buffer_convergence = malloc(size * 12);

    // Write pixel data in chunks, buffering each row for both files
    for (int row = 0, ibnd; row < size;) {
        // If no new lines are available, we wait.
        for (mtx_lock(data->mtx);;) {
            // We extract the minimum of all status variables.
            ibnd = size;
            for (int t = 0; t < num_threads; ++t) {
                if (ibnd > data->status[t].val) {
                    ibnd = data->status[t].val;
                }
            }

            if (ibnd <= row) {
                cnd_wait(data->cnd, data->mtx);
            } else {
                mtx_unlock(data->mtx);
                break;
            }
        }

        // Buffer the pixel data for both files
        for (; row < ibnd; ++row) {
            char *ptr_attractor = buffer_attractor;
            char *ptr_convergence = buffer_convergence;

            for (int col = 0; col < size; col++) {
                // Writing attractor data
                const char *rgb = rgb_values[attractors[row][col]];
                memcpy(ptr_attractor, rgb, strlen(rgb));  // Copy RGB string directly
                ptr_attractor += strlen(rgb);  // Move pointer forward

                // Writing convergence data (greyscale)
                const char *gray = greyscale[convergences[row][col]];
                memcpy(ptr_convergence, gray, strlen(gray));  // Copy grayscale string directly
                ptr_convergence += strlen(gray);  // Move pointer forward
            }

            // Write the buffered rows at once
            *ptr_attractor = '\0';  // Null terminate the buffer
            *ptr_convergence = '\0';  // Null terminate the buffer
            fwrite(buffer_attractor, 1, ptr_attractor - buffer_attractor, fp_attractor);
            fwrite("\n", 1, 1, fp_attractor);  // Newline after each row

            fwrite(buffer_convergence, 1, ptr_convergence - buffer_convergence, fp_convergence);
            fwrite("\n", 1, 1, fp_convergence);  // Newline after each row

            // Free attr/conv row
            free(attractors[row]);
            free(convergences[row]);
        }
    }

    // Free buffers
    free(buffer_attractor);
    free(buffer_convergence);

    // Close files
    fclose(fp_attractor);
    fclose(fp_convergence);

    return thrd_success;
}


//------------------------------------------------------------------------------------------------------------------


int newton_algo(void *arg){

  thread_data_t *data = (thread_data_t *)arg;

  //pre-compute coefficients for the "next z" step at the bottom (newton method)
  double coef = (double)((degree - 1.0)/degree);
  double coef2 = (double)1.0/(degree);

  for (int row = data->id; row < size; row += num_threads) {
    data->convergence = (unsigned char*) malloc(sizeof(unsigned char)*size);
    data->attractor = (unsigned char*) malloc(sizeof(unsigned char)*size);
    for (int col = 0; col < size; col++) {

      double real = MIN_RE_AND_IM + col * step_size;
      double imag = MAX_RE_AND_IM - row * step_size;
      double complex z = real + imag * I;

      int iterations;
      int attractor_index = -1;


      double complex dfz, dfz1, dfz2;
      // Newton's method
      for (iterations = 0; iterations < MAX_ITER; iterations++) {


	// Check if the absolute value of its real or imaginary part exceeds the escape radius
	double z_re = creal(z);
	double z_im = cimag(z);
	if (fabs(z_re) > ESCAPE_RADIUS || fabs(z_im) > ESCAPE_RADIUS) {
	  attractor_index = 11; // Mark as diverged
	  break;
	}

	//Check if close to 0
	double z_abs_sqrd = z_re * z_re + z_im * z_im;
	if (z_abs_sqrd < TOLERANCE_SQUARED){
	  attractor_index = 0;
	  break;
	}

	// Check if the current value is close to unit circle then check same for all roots
        // (1-TOLERANCE)^2 < z_abs_sqrd
	if ( RADIUS < z_abs_sqrd){
	  for (int k = 0; k < degree; k++) {
	    double delta_real = z_re - creal(roots[k]);
	    double delta_imag = z_im - cimag(roots[k]);
	    if ( delta_real*delta_real + delta_imag*delta_imag < TOLERANCE_SQUARED) {

	      attractor_index = k + 1;
	      break;
	    }
	  }
	}


	if (attractor_index != -1) {
	  break;}  // Converged to a root
		   
	// Manual power computation for z^(degree-1)
	switch(degree) {
	  case 1:
	    dfz = 1.0;    // 1
	    break;
	  case 2:
	    dfz =  z;     // z
	    break;
	  case 3:
	    dfz =  z * z;     // z^2
	    break;
	  case 4:
	    dfz1 = z*z; //z^2
	    dfz =  z * dfz1;     // z^3
	    break;
	  case 5:
	    dfz1 = z*z; //z^2
	    dfz =  dfz1 *dfz1;     // z^4
	    break;
	  case 6:
	    dfz1 = z*z; //z^2
	    dfz2 = dfz1 * dfz1; //z^4
	    dfz =  z * dfz2;     // z^5
	    break;
	  case 7:
	    dfz1 = z*z; //z^2
	    dfz2 = dfz1 * dfz1; //z^4
	    dfz =  dfz1 * dfz2;     // z^6
	    break;
	  case 8:
	    dfz1 = z*z; //z^2
	    dfz2 = dfz1 * dfz1; //z^4
	    dfz =  z * dfz1 * dfz2;     // z^7
	    break;
	  case 9:
	    dfz1 = z*z; //z^2
	    dfz2 = dfz1 * dfz1; //z^4
	    dfz = dfz2 * dfz2;     // z^8
	    break;
	  default:
	    fprintf(stderr, "Degree not supported: %d\n", degree);
	    exit(EXIT_FAILURE);
	}
	// Newton's iteration

	z = z * coef + (coef2 / dfz); //z*(degree-1)/degree + ((1/d)/dfz)

      }

      // Check if it took more than MAX_ITER iterations
      if (attractor_index == -1){
	attractor_index = 11;  // Diverged
      }

      // Store attractor and convergence values in local arrays
      data->attractor[col] = (unsigned char)attractor_index;
      data->convergence[col] = (unsigned char)iterations;
    }



    //Update global arrays and status for write thread
    mtx_lock(data->mtx);
    data->status[data->id].val = row + num_threads;
    attractors[row] = data->attractor;
    convergences[row] = data->convergence;
    mtx_unlock(data->mtx);
    cnd_signal(data->cnd);
  }

  return thrd_success;
}



//------------------------------------------------------------------------------------------------------------------


int main(int argc, char *argv[]) {
  if (argc < 4) {
    fprintf(stderr, "Usage: %s -t<num_threads> -l<size> <degree>\n", argv[0]);
    return 1;
  }


  // Parsing command line argum[ents
  for (int i = 1; i < argc - 1; i++) {
    if (argv[i][0] == '-') {
      if (argv[i][1] == 't') {
	num_threads = atoi(&argv[i][2]);
      } else if (argv[i][1] == 'l') {
	size = atoi(&argv[i][2]);
      }
    }
  }

  degree = atoi(argv[argc - 1]); // Last argument is the degree

  //Allocate global arrays and pre-compute stepsize
  step_size = (double)4.0/size;
  attractors = (unsigned char**) malloc(sizeof(unsigned char*)*size);
  convergences =(unsigned char**) malloc(sizeof(unsigned char*)*size);
  

  //Compute roots of unity
  for (int i = 0; i < degree; i++){
    double theta = (double)(2 * M_PI * i) / degree;
    roots[i] = cos(theta) + I * sin(theta);
  }

 //Init stuff needed for threads
  thrd_t write_thread;
  thrd_t threads[num_threads];
  write_data_t write_data;
  thread_data_t thread_data[num_threads];
  int_padded status[num_threads];

  mtx_t mtx;
  cnd_t cnd;
  mtx_init(&mtx, mtx_plain);
  cnd_init(&cnd);


  //Init computing threads
  for (int t = 0; t < num_threads; t++){
    thread_data[t].id = t;
    thread_data[t].mtx = &mtx;
    thread_data[t].cnd = &cnd;
    thread_data[t].status = status;
    status[t].val = -1;



    // Check if thread is created successfully
    int rc = thrd_create(&threads[t], newton_algo, &thread_data[t]);
    if (rc != thrd_success) {
      fprintf(stderr, "Error: Unable to create thread %d\n", t);
      exit(EXIT_FAILURE);  // Exit the program if thread creation fails
    }


  }





  //Init write thread
  snprintf(write_data.filenameAttr, sizeof(write_data.filenameAttr), "newton_attractors_x%d.ppm", degree);
  snprintf(write_data.filenameConv, sizeof(write_data.filenameConv), "newton_convergence_x%d.ppm", degree);
  write_data.mtx = &mtx;
  write_data.cnd = &cnd;
  write_data.status = status;

  int rc = thrd_create(&write_thread, write_ppm_files, (void *)&write_data);
  if (rc != thrd_success){
    fprintf(stderr, "Error: Unable to create write thread\n");
    exit(EXIT_FAILURE);
  }



  // Wait for all threads to complete and combine results. Probably unnecessary since I
  // am already waiting for the write thread (which depends on the compute threads)
  for (int t = 0; t < num_threads; t++) {
    thrd_join(threads[t], NULL);

  }
  // Wait for the write thread to finish before exiting
  int r;
  thrd_join(write_thread, &r);
  if(r != thrd_success){
    fprintf(stderr, "Error: failed to join write thread\n");
    exit(EXIT_FAILURE);
  }

  //Free stuff
  mtx_destroy(&mtx);
  cnd_destroy(&cnd);
  free(convergences);
  free(attractors);

  return 0;
}






