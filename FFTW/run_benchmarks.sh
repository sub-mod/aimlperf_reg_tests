#!/bin/bash

usage() {
    echo "Usage: $0 [-e] [-t] [-l log_filename] [-v thread_values] [-h]"
    echo "  REQUIRED:"
    echo "  -e  Number of executions. Use this value to emulate the number of images processed."
    echo ""
    echo "  OPTIONAL:"
    echo "  -t  Max number of threads to use. Omit this option if you want to use the max number of (real) cores on your system."
    echo "  -l  The resulting log of all the runs will be saved to a file with this name. (Default: fftw_runs.log)"
    echo "  -v  Values of the threads to use. For example, \"2 4 6 8\" will tell this script to run the tests on 2, 4, 6, and 8 threads."
    echo "  -h  Help. Display the usage."
    exit
}

# Set default values
run_log="fftw_runs.log"
max_threads=$(lscpu | awk '/^Core\(s\) per socket:/ {cores=$NF}; /^Socket\(s\):/ {sockets=$NF}; END{print cores*sockets}') #from https://stackoverflow.com/a/31646165
thread_values="-1"

options=":he:t:l:v:"
while getopts "$options" x
do
    case "$x" in
      h)  
          usage
          ;;
      e)
          num_executions=${OPTARG}  
          ;;
      t)  
          max_threads=${OPTARG}
          ;;
      l)  
          run_log=${OPTARG}
          ;;
      v)  
          thread_values=${OPTARG}
          ;;
      *)  
          usage
          ;;
    esac
done
shift $((OPTIND-1))

# Check if the 2d_fft program exists
if [ ! -f 2d_fft ]; then
    echo "The program 2d_fft does not exist! Please compile it by running `. ./compile_benchmark_code.sh /path/to/fftw/lib`"
    exit
fi


# If no thread values were supplied, then iterate in powers of two
if [ "$thread_values" == -1 ]; then
    echo "Using default thread values."
    for (( k=1; k<$max_threads; k*=2 ))
    do
        echo "Executing ./2d_fft $k $num_executions"
	numactl -C 0-$num_threads -i 0,1 ./2d_fft $k $num_executions >> $run_log
    done
    if [ $max_threads > $k ]; then
        echo "Executing ./2d_fft $max_threads $num_executions"
        numactl -C 0-$num_threads -i 0,1 ./2d_fft $max_threads $num_executions >> $run_log
    fi
# Else, use the thread values the user specified
else
    echo "Using custom thread values."
    for k in $thread_values; do
        echo "Executing ./2d_fft $k $num_executions"
        numactl -C 0-$num_threads -i 0,1 ./2d_fft $k $num_executions >> $run_log
    done
fi
