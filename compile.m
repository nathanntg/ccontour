% compile
mex('CFLAGS="\$CFLAGS -Weverything"', 'LDFLAGS="\$LDFLAGS -framework Accelerate"', '-g', '-v', 'ccc.c');
