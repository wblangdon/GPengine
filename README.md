# GPengine
```
Down load source files. Download chaotic time series mg_int_128.dat
Compile with: make -f GPengine.make
Run command line: 
GPengine mg_int_128.dat mg.code0 i2 i3 i4 i5 i6 i7 i8 i9 o10 s02317

mg_int_128.dat training data
i2 ... i9      columns of input training data
o10            column of target output data
mg.code0       small output file
s02317         seed initial population

Linear GP run parameters are code in GPengine.h
