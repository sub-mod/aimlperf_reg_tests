# FFTW Regression Tests

This folder contains regression tests to assess the performance of FFTW's forward and backward FFTs by performing an image blurring task.

# How to Build the Tests

To build the tests,

```
. ./run_me.sh /path/to/main/fftw/folder
```

This command will generate an executable called `2d_fft`.

# How to Run the Tests

To run the tests,

```
./2d_fft <number-of-threads> <number-of-executions>
```

e.g.,

```
./2d_fft 24 2
```

will execute the tests two times spread across twenty four threads.


# Sample Output

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
