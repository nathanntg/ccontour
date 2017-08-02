# Consensus Contours

This is a C implementation of the consensus contour algorithm described by Yoonseob Lim,  Barbara Shinn-Cunningham and Tim Gardner. The [basic algorithm is documented on his website](http://people.bu.edu/timothyg/blog/index.html) as well as in the IEEE [Sparse Contour Representations of Sound](http://ieeexplore.ieee.org/document/6256698/) publication. This implementation was inspired by the [acontour](https://github.com/jmarkow/acontour) repository that contains a MATLAB implementation of the algorithm.

This implementation relies heavily on the macOS Accelerate framework (which uses vectorized instructions to optimize common calculations, such as the FFT and the complex ratio required for reassignment).

Currently, this code has a few limitations:

* The original consensus contour work suggests eliminating all but the longest contours. In order to reduce memory usage and extend this method to realtime audio processing, this code preserves all contours.
* This just returns a consensus contour spectrogram, but does not provide programmatic access to the identified contours.

Advantages to this code:

* The code supports both single and double precision.
* The code includes a **mex** file implementation to allow usage in MATLAB.
* For a 5 second audio file, the `ccontour` implementation is 68x faster than the `acontour` (MATLAB) implementation. On most modern computers, it is sufficiently fast to calculate contours in realtime.


## Usage: C

The API is modeled off the Accelerate framework, specifically using an opaque pointer to a structure to hold all the memory and pre-calculated values required for the algorithm. During setup, there is a configuration phase where a mutable opaque pointer is updated with configuration information used by the algorithm. Next, the configuration is turned into a setup, which acts as an immutable opaque pointer that has memory and pre-calculated values necessary for implementing the algorithm.

To get a default configuration, calll:

```c
CCCConfig createCCCConfig();
```

This configuration is a good, out-of-the box setup and can be used immediately. Alternatively, you can specify any of the following values:

```c
void cccConfigSetFFTLength(CCCConfig config, unsigned long fft_length);
// default: 1,024
// must be a power of 2

void cccConfigSetFFTOverlap(CCCConfig config, unsigned long fft_overlap);
// default: 1,005

void cccConfigSetWeightByPower(CCCConfig config, bool pow_weight);
// default: true

void cccConfigSetSampleRate(CCCConfig config, float fs);
// default: 44,100 Hz

void cccConfigSetTimescales(CCCConfig config, unsigned long num_timescales, const float timescales[]);
// default: 0.5ms, 0.7ms, 0.9ms, 1.1ms, 1.3ms, 1.5ms, 1.7ms, 1.9ms, 2.1ms

void cccConfigSetAngles(CCCConfig config, unsigned long num_angles, const float angles[]);
// default: [eight angles in radians composing full circle]
```

Once you have completed configuration, you can turn the configuration into an immutable setup that allocates all memory needed for the algorithm and pre-calculates a number of values.

To perform the setup:

```c
CCCSetup createCCCSetup(const CCCConfig config);
```

During the setup process, any necessary values are copied and the configuration is no longer required. You can de-allocate the configuration by calling:

```c
void destroyCCCConfig(CCCConfig config);
```

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

Finally, when you are done with the setup, you can release all allocated memory by calling:

```c
void destroyCCCSetupD(CCCSetup setup);
```

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

Where `signal` is a 1D vector of type `single` or `double`, and `fs` is the sample rate in Hz. You can optionally pass additional configuration parameters as character value pairs. For help on all options available via the MATLAB interface, use the help function:

```matlab
help ccontour
```

