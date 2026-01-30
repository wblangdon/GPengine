#WBL 23 June 2025 Run default GPengine for Mackey-Glass

#Modifications:
#WBL 13 Dec 2025 (cf ten8 r1.5) Add best.n remove stats.awk pipe
#WBL 12 Aug 2025 Remove mg.code? Use all CPU cores
#WBL  4 Aug 2025 pipe into stats.awk, restore mg.code?
#WBL 23 Jul 2025 Add time
#WBL 22 Jul 2025 Remove mg.code? Use all CPU cores
#WBL 28 Jun 2025 Add $1

echo $0 '$Revision: 1.14 $' "start" `date` `pwd` $HOST

g++ --version
rcsdiff ../main.cpp ../GPengine.h ../GPengine.cpp ../entropy.cpp ../GPengine.make ../stats.awk

if( $1 == 0 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.0 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s02317 >& mg.0
if($status) exit $status;
endif

if( $1 == 1 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.1 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s12317 >& mg.1
if($status) exit $status;
endif

if( $1 == 2 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.2 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s22317 >& mg.2
if($status) exit $status;
endif

if( $1 == 3 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.3 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s32317 >& mg.3
if($status) exit $status;
endif

if( $1 == 4 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.4 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s42317 >& mg.4
if($status) exit $status;
endif

if( $1 == 5 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.5 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s52317 >& mg.5
if($status) exit $status;
endif

if( $1 == 6 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.6 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s62317 >& mg.6
if($status) exit $status;
endif

if( $1 == 7 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.7 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s72317 >& mg.7
if($status) exit $status;
endif

if( $1 == 8 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.8 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s82317 >& mg.8
if($status) exit $status;
endif

if( $1 == 9 ) then
time GPengine ~/mackey_glass/mg_int_128.dat best.9 \
  i2 i3 i4 i5 i6 i7 i8 i9 o10 s92317 >& mg.9
if($status) exit $status;
endif


echo "$0 $1 done" `date`
