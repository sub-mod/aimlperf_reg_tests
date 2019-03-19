# FFTW Regression Tests

This folder contains regression tests to assess the performance of FFTW's forward and backward FFTs by performing an image blurring task.

## Building FFTW

To build FFTW, either use Podman to build the Dockerfile under the "Dockerfiles" folder, or run the scripts in the "installation" folder yourself. Both methods utilize rpmbuild to build the library.

### Podman

If you decide to use Podman, make note of the `FFTW_CFLAGS` and `FFTW_BUILD_FLAGS` environment variables in the Dockerfile. These variables determine which build flags FFTW will be configured to use. Once you've chosen your flags, build the image via:

```
$ podman build -f Dockerfiles/Dockerfile.[name]
```

where `[name]` is either `rpmbuild_only` or `rpmbuild_regression_tests`. The first, as the name says, only builds the FFTW rpms. The second also builds the FFTW rpms, but then adds and runs the regression tests.

### Installation Scripts

If you would like to use the installation scripts, you can run either `build.sh` or `install.sh` to build the library, but of course `install.sh` will also install the generated rpms for you. See

```
$ sh build.sh -h
```

for help on how to use the build script. Similarly, use

```
$ sh install.sh -h
```

for help on how to use the install script.

Note that if you wish to build FFTW while using root, you will need to use the `-s` flag in both scripts to override the automatic exit. (The exit is for protection.)

## How to Build the Tests

To build the tests,

```
$ . ./compile_benchmark_code.sh /path/to/main/fftw/folder
```

This command will generate two executables: `2d_fft` and `nd_cosine_ffts`.

## How to Run the Tests

### Using the Existing Shell Script

To run the tests with the `run_benchmarks.sh` script,

```
$ ./run_benchmarks.sh -i <number_of_iterations> -e <executable> [options]
```

where `<executable>` is either `2d_fft` or `nd_cosine_ffts`.

Note that for `nd_cosine_ffts`, there are additional required flags. See `./run_benchmarks.sh -h` for help and more info on the options available.

This tool is useful because it allows you to run a series of benchmark tests with one command, rather than running everything by hand multiple times. It also saves all of the runs to a log file, defined by the user with the `-l` option.

Example for `2d_fft`:

```
$ ./run_benchmarks.sh -i 1000 -e 2d_fft -v "2 4 9"
```

This command will run the `2d_fft` executable on 2, 4, and 9 threads, and each time it will process 1000 images.

Example for `nd_cosine_ffts`:

```
$ ./run_benchmarks.sh -e nd_cosine_ffts -i 14 -r 2 -f 0.001 -d "300 300" -j "fftw_cosine_performance_results.json"
```

This command will run the `nd_cosine_ffts` on fourteen 2D cosine matrices--each dimension being 300 in size--with a sampling frequency `Fs=0.001`, and it will save the performance results to `fftw_cosine_performance_results.json`. Since no thread values were specified, the thread values will be powers of 2, up to the maximum number of real cores on your system.

To optimize performance, you can use the `-n` option to enable numactl.

### Running by Hand

To run the image blurring test by hand,

```
$ ./2d_fft <number-of-threads> <number-of-executions>
```

e.g.,

```
$ ./2d_fft 24 2
```

will execute the tests two times spread across twenty four threads.

To run the cosine FFT tests by hand,

```
$ ./nd_cosine_ffts "noplot" "test.json" 24 10 0.00001 2 30000 30000 
```

The `"noplot"` parameter tells the executable not to plot the results. If you want to plot the results, however, change `"noplot"` to `"plot"`--but make sure you have gnuplot installed! Look for `test.json` to see performance results.

If you want a quick rundown of parameter info, simply run

```
$ ./nd_cosine_ffts
```

This will throw an error, but the error will tell you all the parameters that are required and in what order.


## Sample Outputs

Below are sample outputs from each FFTW test set.

### Image Blurring

```
PERFORMANCE RESULTS
===================
Operations:
    2 images of size 2392x2500 analyzed
    24 threads used
FFT Performance Results
    26.112 FFT performance GFlops
    0.077 sec FFT execution time
Inverse FFT Performance Results
    23.066 IFFT performance GFlops
    0.087 sec IFFT execution time
FFT + IFFT Setup time
    Took 0.274 sec to setup 2 images
    Took 0.137 sec to setup a single image
Blur time (non-FFTW computations)
    Took 0.098 sec to blur 2 images
    Took 0.049 sec to setup a single image
Wall time
    Took 0.535 sec to blur 2 images (includes non-FFTW computations)
    Took 0.268 sec to blur single image (includes non-FFTW computations)
Wall time (excluding blur time)
    Took 0.437 sec to blur 2 images (only FFTW computations)
    Took 0.219 sec to blur single image (only FFTW computations)
```

### Cosine FFTs

```
PERFORMANCE RESULTS
===================
Input Info:
    One 2D cosine: 30000 x 30000 samples
    fs = 1.00e-03 Hz
    14 iterations
    24 threads used
DFT Results
    Forward DFT execution time: 1.698 sec 
    Forward DFT TFlops: 3591.624
    Backward DFT execution time: 1.784 sec 
    Backward DFT TFlops: 3419.213

```
