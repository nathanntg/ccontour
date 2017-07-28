# Consensus Contours

This is a C implementation of the consensus contour algorithm developed by Dr. Tim Gardner. The [basic algorithm is documented on his website](http://people.bu.edu/timothyg/blog/index.html) as well as in the IEEE [Sparse Contour Representations of Sound](http://ieeexplore.ieee.org/document/6256698/) publication. This implementation was inspired by the [acontour](https://github.com/jmarkow/acontour) repository that contains a MATLAB implementation of the algorithm.

This implementation relies heavily on the macOS Accelerate framework (which uses vectorized instructions to optimize common calculations, such as the FFT and the complex ratio required for reassignment).

Currently, this code has a few limitations:

* Parameters are hard coded in C (timescales, angles, etc). Eventually a configuration structure will allow overriding these defaults.
* The original consensus contour work suggests eliminating all but the longest contours. In order to reduce memory usage and extend this method to realtime audio processing, this code preserves all contours.

Advantages to this code:

* The code supports both single and double precision.
* The code includes a **mex** file implementation to allow usage in MATLAB.
* For a 5 second audio file, the `ccontour` implementation is 68x faster than the `acontour` (MATLAB) implementation. On most modern computers, it is sufficiently fast to calculate contours in realtime.


## Usage: C

The API is modeled off the Accelerate framework, specifically using an opaque pointer to a structure to hold all the memory and pre-calculated values required for the algorithm.

To perform the setup:

```c
CCCSetup createCCCSetup(unsigned long fft_length, unsigned long fft_overlap, float fs, bool pow_weight);
```

The returned `CCCSetup` is an opaque structure that contains memory for all temporary buffers, all timescales and angles. The create setup function takes an `fft_length` (must be a power of 2), the `fft_overlap` (only applies to calculating spectrograms), the sample rate (`fs`) in Hz and whether or not to weight the returned contours by the spectral power.

You can then use the `CCCSetup` variable to calculate contours for a single column (the setup will hold information from the previous column, so you will need one setup per signal that will be processed simultaneously). To do so, pass a pointer to the start of the signal (where the signal has at least `fft_length` elements) and a pointer to the output vector (that has `fft_length / 2` elements):

```c
void cccColumn(const CCCSetup setup, const float *signal, float *contour);
```

You can also calculate the full spectrogram for a signal in one call. The program provides a `struct ConsensusContourSize` type that contains size information to help in pre-allocating the output memory. To get the required dimensions, call:

```c
struct ConsensusContourSize dim = cccSize(const CCCSetup setup, const unsigned long signal_len);
```

The returned `struct ConsensusContourSize` contains fields `dim.bytes` telling you how many bytes are required for the full spectrogram, as well as `dim.rows` (the number of rows based on the `fft_length` and `fft_overlap`) and the `dim.cols` (equal to `fft_length / 2`).

Once you allocate an output vector with sufficient space (see `bytes` above), you can calculate the spectrogram all at once:

```C
void cccSpectrogram(const CCCSetup setup, const struct ConsensusContourSize dim, const float *signal, float *consensus_contours);
```

Behind the scenes, this calls `cccColumn` while moving through the signal. 

Note that both `cccColumn` and `cccSpectrogram` will zero out the output vector before proceeding.

## Usage: MATLAB

To use the MATLAB mex function, you must first compile it. Again, because of the dependency on the Accelerate framework, this is only supported for macOS. 

To compile, run:

```matlab
compile_ccontour_mex
```

Once compiled, you can run:

```matlab
consensus_contours = ccontour(signal, fs);
```

Where `signal` is a 1D vector of type `single` or `double`, and `fs` is the sample rate in Hz.