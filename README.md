# GPengine
```
Down load source files. Download chaotic time series mg_int_128.dat
Compile with: make -f GPengine.make
Run command line: 
GPengine mg_int_128.dat mg.code0 i2 i3 i4 i5 i6 i7 i8 i9 o10 s02317 t19

mg_int_128.dat training data
i2 ... i9      columns of input training data
o10            column of target output data
mg.code0       small output file
s02317         seed initial population
t19            nineteen threads

Linear GP run parameters are code in GPengine.h

GPengine came orginally from Peter Nordin and was used
in work with his "Elvis" humanoid robot
http://gpbib.cs.ucl.ac.uk/gp-html/langdon_2001_elvis.html
and later as the GP system investigating repeated patterns
created by crossover in evolved linear programs
http://gpbib.cs.ucl.ac.uk/gp-html/langdon_2005_CS.html

Most recently GPengine has been used to power
long term evolution experiments (LTEE) in linear genetic programming
http://gpbib.cs.ucl.ac.uk/gp-html/Langdon_2025_IMOL.html
and http://gpbib.cs.ucl.ac.uk/gp-html/Langdon_2026_raLGP.html

For speed, the 2026 version assumes your Linux computer has Intel's AVX512 vector instructions
and so GPengine.make compiles with g++ -march=skylake-avx512

If your computer does not support avx512 GPengine may terminate with exit status 132, eg
```
./GPengine mg_int_128.dat i2 i3 i4 i5 i6 i7 i8 i9 o10 s92317
#GPengine $Revision: 1.25 $ rev=1.129x $ AVX512 re2=1.14 $ WBL December 2025 ./GPengine mg_int_128.dat i2 i3 i4 i5 i6 i7 i8 i9 o10 s92317 seed=92317 threads=8 GenerateLimit=249500 Fri Jan 30 10:47:21 2026

Illegal instruction (core dumped)
echo $status
132
```
