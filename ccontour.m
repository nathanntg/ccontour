%CCONTOUR Help file for the ccontour.c MEX file
%   This function calculates consensus contours for an input signal.
%   Consensus contours are described in a paper cited below by Yoonseob 
%   Lim,  Barbara Shinn-Cunningham and Tim Gardner. The algorithm
%   calculates reassigned spectrograms across a range of timescales and
%   then looks for consensus across a range of angular contours.
%
%   This function returns a spectrogram style summary of the consensus
%   contour representation, which provides a more precise representation.
%
%   SPECT = CCONTOUR(SIGNAL, FS) calculates the default consensus contour
%   spectrogram for a vector SIGNAL (e.g., audio) recorded with a sampling
%   rate of FS samples per second (Hz). The returned matrix SPECT will
%   contain rows for each frequency bin and columns in accordance with the
%   length of the signal.
%
%   The function has two additional optional outputs for the frequency and
%   time bins that correspond with the spectrogram rows and columns
%   respectively:
%
%   [SPECT, F, T] = CCONTOUR(...);
%
%   Optional parameters can be passed as string value pairs. The parameters
%   are (in order of likelihood of use):
%
%   CCONTOUR(..., 'fft_length', FFT_LENGTH) is the length of the window
%   used when calculating the reassigned spectrograms. The default is 1024.
%   This determines the number of frequency bins (rows) in the returned
%   spectogram. This must be a power of two.
%
%   CCONTOUR(..., 'fft_shift', FFT_OVERLAP) is the shift in the signal
%   between columns in the returned spectrogram. The default is 19.
%
%   CCONTOUR(..., 'pow_weight', POW_WEIGHT) allows you to determine
%   whether the returned image should have contours weighted by their
%   spectral power. This defaults to true.
%
%   CCONTOUR(..., 'timescales', TIMESCALES) specify the different
%   timescales when calculating the reassigned spectograms. TIMESCALES 
%   should be a vector of the same type as SIGNAL. This defaults to 
%   0.5ms, 0.7ms, 0.9ms, 1.1ms, 1.3ms, 1.5ms, 1.7ms, 1.9ms, 2.1ms.
%
%   CCONTOUR(..., 'angles', ANGLES) specify the different
%   angles when identifying contours in the reassigned spectrograms. 
%   ANGLES should be a vector of the same type as SIGNAL. This defaults to
%   eight angles spread over a full circle.
%
%   Citation: http://ieeexplore.ieee.org/document/6256698/
%   C library: https://github.com/nathanntg/ccontour
% 
