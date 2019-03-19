#define _DEFAULT_SOURCE
#define PI 3.141592653589793238462643383279
#define TIMELIMIT 2

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <fftw3.h>
#include <stdbool.h>
#include <string.h>


void generate_cosine_data(double *cosine, double fs, int rank, int *n, int matrix_size);
void fill_row(double *cosine, double fs, int row_length, int start_idx, int n_sum, int matrix_size);
void plot1D(double *cosine, int dim, int rank, int *n, char *title);

int main(int argc, char* argv[]){

    // Loop variables
    int i,j;

    // Parse inputs
    bool plot; //to plot or not to plot -- that is the question!
    char *plot_opt; //user passes either "plot" or "noplot" for plotting
    double fs; //sampling frequency (double values)
    int nthreads, niters, rank;
    int n[100]; //will hold all of the rank data... max of 100 dims
    char *pEnd;
    if (argc == 1){
        printf("No arguments were passed! Please enter: (1.) \"noplot\" or \"plot\" for plotting, (2.) number of threads to use, (3.) number of iterations to execute, (4.) the sampling frequency \"fs\" for the cosine, (5.) the rank of the cosine, and (6.) the size of each dimension.\n");
        exit(0);
    }
    else if (argc < 5){
        printf("Only %d argument(s) given. Minimum number of arguments is 6. Please enter: (1.) \"noplot\" or \"plot\" for plotting, (2.) number of threads to use, (3.) number of iterations to execute, (4.) the sampling frequency \"fs\" for the cosine, (5.) the rank of the cosine, and (6.) the size of each dimension", argc-1);
        exit(0);
    }
    else{
        plot_opt = argv[1]; //set to "plot" to plot or "noplot" to not plot

        if (strcmp(plot_opt, "plot") == 0)
            plot=true;
        else if (strcmp(plot_opt, "noplot") == 0)
            plot=false;
        else{
            printf("Invalid input '%s' for plotting. Please use \"plot\" or \"noplot\"\n", plot_opt);
            exit(0);
        }

        nthreads = (int)strtol(argv[2], &pEnd, 10);
        niters = (int)strtol(argv[3], &pEnd, 10);
        fs = atof(argv[4]);
        rank = (int)strtol(argv[5], &pEnd, 10);

        if (argc <= rank+5){
            printf("Rank is set to %d, but %d dimensions were passed. The number of dimensions passed must equal the rank. Exiting now.\n", rank, (rank+5-argc));
            exit(0);
        }

        for (i=6; i<rank+6; i++){
            n[i-6] = (int)strtol(argv[i], &pEnd, 10);
        }

        if (nthreads < 1){
            printf("Number of threads must be greater than or equal to 1.\n");
            exit(0);
        }
        if (niters < 1){
            printf("Number of iterations must be greater than or equal to 1.\n");
            exit(0);
        }
        if (fs <= 10e-9){
            printf("Sampling frequency must be greater than 0.0. You entered: %0.2e\n", fs);
            exit(0);
        }
        if (rank < 1){
            printf("The rank must be greater than or equal to 1.\n");
            exit(0);
        }
    }

    // Cosine variables
    int n_total = 1;

    // Plot variables
    char *title = "Resulting cosine Curve After Forward and Backward DFTs";

    // FFTW variables
    unsigned flags = FFTW_ESTIMATE;

    // Performance variables
    struct timeval forward_dft_start, forward_dft_stop;
    struct timeval backward_dft_start, backward_dft_stop;
    double forward_dft_execution_time = 0.0;
    double backward_dft_execution_time = 0.0;
    double total_f_dft_exec_time = 0.0;
    double total_b_dft_exec_time = 0.0;

    // Get the total number of indices
    for (i=0; i<rank; i++){
        n_total *= n[i];
    }

    // Set threading
    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads);

    // Allocate memory for cosine data
    double *cosine = (double*)malloc(n_total * sizeof(double));

    // Fill N-dimensional cosine matrix
    generate_cosine_data(cosine, fs, rank, n, n_total);

    // Set time limit so that FFTW doesn't spend too much time trying to figure out the "best" algorithm.
    fftw_set_timelimit(TIMELIMIT);

    // Initialize real-to-complex cosine input and output
    double *cosine_original = (double*)fftw_malloc(n_total * sizeof(double));
    fftw_complex *cosine_complex = (fftw_complex*)fftw_malloc(n_total * sizeof(fftw_complex));

    // Initialize the cosine that will be returned from the complex DFT
    double *cosine_back = (double*)fftw_malloc(n_total * sizeof(fftw_complex));

    // Iterate
    for (j=0; j<niters; j++){
        // Create FFTW plans
        fftw_plan forward_cos_dft_plan = fftw_plan_dft_r2c(rank, n, cosine_original, cosine_complex, flags);
        fftw_plan backward_cos_dft_plan = fftw_plan_dft_c2r(rank, n, cosine_complex, cosine_back, flags);

        // Fill input cosine array (this MUST be done after the fftw plans are created)
        for (i=0; i<n_total; i++)
            cosine_original[i] = cosine[i];

        // Execute Forward DFT and capture performance time
        gettimeofday(&forward_dft_start, NULL); //start clock
        fftw_execute(forward_cos_dft_plan);
        gettimeofday(&forward_dft_stop, NULL); //stop clock
        forward_dft_execution_time = (forward_dft_stop.tv_sec - forward_dft_start.tv_sec) * 1000.0;// sec to ms
        forward_dft_execution_time += (forward_dft_stop.tv_usec - forward_dft_start.tv_usec)/ 1000.0;// us to ms
        forward_dft_execution_time *= (1.0e-3); 
        total_f_dft_exec_time += forward_dft_execution_time;

        // Execute Backward DFT and capture performance time
        gettimeofday(&backward_dft_start, NULL); //start clock
        fftw_execute(backward_cos_dft_plan);
        gettimeofday(&backward_dft_stop, NULL); //stop clock
        backward_dft_execution_time = (backward_dft_stop.tv_sec - backward_dft_start.tv_sec) * 1000.0;// sec to ms
        backward_dft_execution_time += (backward_dft_stop.tv_usec - backward_dft_start.tv_usec)/ 1000.0;// us to ms
        backward_dft_execution_time *= (1.0e-3); 
        total_b_dft_exec_time += backward_dft_execution_time;

        // Destroy FFTW plans
        fftw_destroy_plan(forward_cos_dft_plan);
        fftw_destroy_plan(backward_cos_dft_plan);
    }

    //for (i=n[0]; i<n[0]+n[1]; i++){
    //    printf("cosine[%d] = %0.2f\n", i, cosine[i]);
    //}

    // Handle threading
    fftw_cleanup_threads();

    // Get average times
    double average_forward_dft_exec_time = total_f_dft_exec_time / niters;
    double average_backward_dft_exec_time = total_b_dft_exec_time / niters;

    // Compute teraflops (see here for info on how to calculate mflops: http://www.fftw.org/speed/)
    long double forward_dft_mflops_approx = 5 * n_total * log2l(n_total) / average_forward_dft_exec_time;
    long double backward_dft_mflops_approx = 5 * n_total * log2l(n_total) / average_backward_dft_exec_time;

    long double forward_dft_tflops_approx = forward_dft_mflops_approx * (1e-6);
    long double backward_dft_tflops_approx = backward_dft_mflops_approx * (1e-6);

    // Plot result to ensure we get back what we put in!
    if (plot == true)
        plot1D(cosine_back, 1, rank, n, title);

    // Create file to save results to
    FILE *results_file = fopen("fftw_cosine_performance_results.json", "w");

    // Get timestamp
    time_t raw_time = time(NULL);
    struct tm *timeinfo;
    timeinfo = localtime(&raw_time);

    // Save as JSON
    fprintf(results_file, "{\n");
    fprintf(results_file, "    \"timestamp\": \"%d-%d-%d %d:%d:%d\",\n", timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
    fprintf(results_file, "    \"performance_results\": {\n");
    fprintf(results_file, "        \"inputs\": {\n");
    fprintf(results_file, "            \"rank\": %d,\n", rank);
    fprintf(results_file, "            \"dims\": [");
    for (i=0; i<rank-1; i++){
        fprintf(results_file, " %d,", n[i]);
    }
    fprintf(results_file, " %d],\n", n[rank-1]);
    fprintf(results_file, "            \"fs_Hz\": %0.2e,\n", fs);
    fprintf(results_file, "            \"iterations\": %d,\n", niters);
    fprintf(results_file, "            \"threads\": %d\n", nthreads);
    fprintf(results_file, "        },\n");
    fprintf(results_file, "        \"forward_dft_results\": {\n");
    fprintf(results_file, "            \"average_execution_time_seconds\": %0.5f,\n", average_forward_dft_exec_time);
    fprintf(results_file, "            \"average_tflops\": %0.5Lf\n", forward_dft_tflops_approx);
    fprintf(results_file, "        },\n");
    fprintf(results_file, "        \"backward_dft_results\": {\n");
    fprintf(results_file, "            \"average_execution_time_seconds\": %0.5f,\n", average_backward_dft_exec_time);
    fprintf(results_file, "            \"average_tflops\": %0.5Lf\n", backward_dft_tflops_approx);
    fprintf(results_file, "        }\n");
    fprintf(results_file, "    }\n");
    fprintf(results_file, "}\n");

    fclose(results_file);

    printf("\nPERFORMANCE RESULTS\n");
    printf("===================\n");
    printf("Input Info:\n");
    printf("    One %dD cosine: %d", rank, n[0]);
    for (i=1; i<rank; i++)
        printf(" x %d", n[i]);
    printf(" samples\n");
    printf("    fs = %0.2e Hz\n", fs);
    printf("    %d iterations\n", niters);
    printf("    %d threads used\n", nthreads);
    printf("DFT Results\n");
    printf("    Forward DFT execution time: %0.3f sec\n", average_forward_dft_exec_time);
    printf("    Forward DFT TFlops: %0.3Lf\n", forward_dft_tflops_approx);
    printf("    Backward DFT execution time: %0.3f sec\n", average_backward_dft_exec_time);
    printf("    Backward DFT TFlops: %0.3Lf\n", backward_dft_tflops_approx);

    return 0;
}

void fill_row(double *cosine, double fs, int row_length, int start_idx, int n_sum, int matrix_size){
/* Helper function to fill a row of data in an N-dimensional cosine matrix
 *
 * Inputs
 * ======
 *   double *cosine
 *       The array to be filled with cosine data
 *
 *   int fs
 *       Sample size for the cosine
 *
 *   int row_length
 *       Length of the dimension (AKA the length of the row to fill)
 *
 *   int start_idx
 *       Index of the cosine array to start filling
 *
 *   int n_sum
 *       Sum of all the values in n
 *
 *   int matrix_size
 *       Total number of samples in the "cosine" matrix across all dimensions (i.e., n0*n1*n2*...*nK)
 */
    if (start_idx < matrix_size){
        int i; //iterative var
        for (i=0; i<row_length; i++){
            cosine[i+start_idx] = cos(i*fs*180/PI);
        }

        fill_row(cosine, fs, row_length, start_idx+n_sum, n_sum, matrix_size);
    }
}

void generate_cosine_data(double *cosine, double fs, int rank, int *n, int matrix_size){
/* Generates data for a forward FFT
 *
 * Inputs
 * ======
 *   double *cosine
 *       The array to be filled with cosine data
 *
 *   int fs
 *       Sample size for the cosine
 *
 *   int rank
 *       Number of dimensions in the data array
 *
 *   int *n
 *       An array which contains the dimensions of the data array
 *
 *   int matrix_size
 *       Total number of samples in the "cosine" matrix across all dimensions (i.e., n0*n1*n2*...*nK)
 */

    // Init values
    int i; //iterative value
    int start = 0; //start of the next dimension
    int prev = 0; //start of the previous dimension
    int dim; //current dimension
    int n_sum = 0; //sum of all the dimensions

    // Get the sum of all N values. This will help us with indexing the cosine matrix when we go to
    // fill it with data
    for (i=0; i<rank; i++){
        n_sum += n[i];
    }

    // Fill cosine array
    for (i=0; i<rank; i++){

        // Get current dimension
        dim = n[i];

        // Get starting point
        start += prev;

        // Fill array
        fill_row(cosine, fs, dim, start, n_sum, matrix_size);

        // Update 'prev'
        prev = dim;
    }
}

void plot1D(double *cosine, int dim_to_plot, int rank, int *n, char *title){
/* Plot cosine data for a specific dimension
 *
 * Inputs
 * ======
 *   double *cosine
 *       The array to be filled with cosine data
 *
 *   int dim_to_plot
 *       Dimension of the nD cosine to plot
 *
 *   int rank
 *       Number of dimensions in the data array
 *
 *   int *n
 *       An array which contains the dimensions of the data array
 *
 *   char *title
 *       Title of the plot
 */
    // Get the length of the dimension
    int N = n[dim_to_plot];

    // Init values
    int i; //iterative value
    int n_sum = 0; //sum of all the values in n
    int row_start_idx = 0; //Start and end of the row for the given dimension
    double *xvals = (double*)malloc(N * sizeof(double)); //allocate memory for x values
    double *yvals = (double*)malloc(N * sizeof(double)); //allocate memory for y values

    // Get the starting point and the sum of all values in n
    for (i=1; i<rank-1; i++){
        n_sum += n[i];

        if (dim_to_plot <= rank)
            row_start_idx += n[i-1];
    }

    // Now gather data and store in x- and y-value arrays
    for (i=0; i<N; i++){

        // Get x-values (simply just store index i)
        xvals[i] = i;

        // Get y-values
        yvals[i] = cosine[i+row_start_idx];
        row_start_idx += n_sum;
    }

    // File to save data in
    FILE *cosine_data_file = fopen("cosine_data.txt", "w");
    for (i=0; i<N; i++){
        fprintf(cosine_data_file, "%lf %lf \n", xvals[i], yvals[i]);
    }

    // Open gnuplot
    FILE *gnuplot_pipe = popen("gnuplot -persistent", "w");

    // Prepare title (no more than 100 chars)
    char plot_title[100];
    sprintf(plot_title, "set title \"%s\"", title);

    // Prepare gnuplot commands
    char *gnuplot_cmds[] = {plot_title, "plot 'cosine_data.txt' with line"};

    // Now plot
    for (i=0; i<2; i++){
        fprintf(gnuplot_pipe, "%s \n", gnuplot_cmds[i]);
    }
}
