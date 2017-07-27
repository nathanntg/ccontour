% compile
%mex('CFLAGS="\$CFLAGS -Weverything"', 'LDFLAGS="\$LDFLAGS -framework Accelerate"', '-g', '-v', 'ccc.c', 'consensus_contour.c');
mex('CFLAGS="\$CFLAGS -Weverything -O3"', 'LDFLAGS="\$LDFLAGS -framework Accelerate"', 'ccc.c', 'consensus_contour.c');
