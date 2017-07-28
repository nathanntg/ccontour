function compile_ccontour_mex(varargin)
%COMPILE_DT_MEX Compile MEX versions of dynamic time warping functions
%   In addition to MATLAB versions of the dynamic time warping code, this
%   repository offers MEX versions, which are much preferred. The MEX
%   versions offer close to a 100x speed benefit and are designed to use
%   substantially less memory.
%
%   The repository includes Mac versions of the MEX files already compiled,
%   but this function allows compiling the MEX files for other platforms.
%
%   If you see warnings about non-MEX implementations when running code
%   from this repository, call this function once to compile the MEX
%   functions.
%
%   No parameters are required, but you can optionally provide parameters:
%
%   COMPILE_DT_MEX('debug', true) turns off optimizations, adds debugging
%   symbols to the binaries and turns on verbose compilation, all to help
%   debug issues with the MEX files.
%
%   COMPILE_DT_MEX('warnings', true) turns on all warnings during compile
%   time (only tested with Clang compiler).

%% parameters
debug = false;
warnings = false;

% load custom parameters
nparams = length(varargin);
if 0 < mod(nparams, 2)
    error('Parameters must be specified as parameter/value pairs');
end
for i = 1:2:nparams
    nm = lower(varargin{i});
    if ~exist(nm, 'var')
        error('Invalid parameter: %s.', nm);
    end
    eval([nm ' = varargin{i+1};']);
end

% check platform
if ~ismac
    error('Only supports the macOS platform (requires Accelerate framework).');
end

% print nice message
fprintf('Compiling ccontour functions...\n');

c = {};
cf = {'-DPRECISION=2'};

% show warnings
if warnings
    cf{end + 1} = '-Weverything"';
end

% enable debugging
if debug
    c{end + 1} = '-g';
    c{end + 1} = '-v';
else
    c{end + 1} = '-silent';
    cf{end + 1} = '-O3';
end

% build cflags
if numel(cf) > 0
    c{end + 1} = ['CFLAGS="\$CFLAGS ' strjoin(cf) '"'];
end

% include Accelerate framework
if ismac
    c{end + 1} = 'LDFLAGS="\$LDFLAGS -framework Accelerate"';
end

% call mex functions
functions = {{'ccontour.c', 'consensus_contour.c'}};
for j = 1:length(functions)
    if iscell(functions{j})
        fprintf('%s\n', functions{j}{1});
        d = [c functions{j}];
    else
        fprintf('%s\n', functions{j});
        d = [c [functions{j} '.c']];
    end
    mex(d{:});
end
