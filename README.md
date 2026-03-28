# GPengine: Linear Genetic Programming
GPengine has been used in a number of experiments.
It was originally written by Peter Nordin at the end of the last century.

At the end of last year it was upgraded to take advantage of modern Intel X86
Linux pthreads to allow parallel fitness evaluation for use with Long Term Evolution Experiments to 100000 generations
and bloated programs of up to 4 million instructions.

In some cases, the use of multi-core parallism was still not sufficient and
so GPengine was upgraded to make use of Intel Vextor eXtentions.
First to use SSE 256 bit instructions
and (the code here) secondly to use AVX 512 bit instrctions.

# January 2026 Linux parallel pthreads and Intel AVX-512 vector instructions
# To use
Down load source files. Download chaotic time series mg_int_128.dat

Compile with: make -f GPengine.make

Run command line, eg: 

GPengine mg_int_128.dat mg.code0 i2 i3 i4 i5 i6 i7 i8 i9 o10 s02317 t19

* mg_int_128.dat training data
* i2 ... i9      columns of input training data
* o10            column of target output data
* mg.code0       small output file
* s02317         pseudo random number generator (PRNG) seed
* t19            nineteen threads

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

The Genetic Improvement experiments are documented in:
* SSE 256 bits http://gpbib.cs.ucl.ac.uk/gp-html/langdon_2026_GI.html
* AVX 512 bits https://arxiv.org/abs/2512.09157

# AVX512 not available, error: exit status 132

For speed, the 2026 version of GPengine assumes your Linux computer has Intel's AVX512 vector instructions
and so GPengine.make compiles with g++ -march=skylake-avx512

If your computer does not support avx512 GPengine may terminate with exit status 132, eg
```
./GPengine mg_int_128.dat i2 i3 i4 i5 i6 i7 i8 i9 o10 s92317
#GPengine $Revision: 1.25 $ rev=1.129x $ AVX512 re2=1.14 $ WBL December 2025 ./GPengine mg_int_128.dat i2 i3 i4 i5 i6 i7 i8 i9 o10 s92317 seed=92317 threads=8 GenerateLimit=249500 Fri Jan 30 10:47:21 2026

Illegal instruction (core dumped)
echo $status
132
```
