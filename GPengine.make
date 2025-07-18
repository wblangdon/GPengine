############################################################################
#
#       Project		Linux GPengine $Revision: 1.6 $ 
#
#       Author		W.B.Langdon
#
#	Created		19 June 2025 from gp/ei2022/lung.make r1.2
#
#Modifications:
#WBL 5 Jul 2025 add AVX etc based on avx.make r1.9
############################################################################

NAME 	= GPengine
OBJS 	= main.o GPengine.o
#LIBS 	= -lm
CC 	= g++
COMPILE = $(CC) -c $(CFLAGS)
CFLAGS 	= -O3
#CFLAGS += -Delvis
CFLAGS 	+= -g
CFLAGS += -DNDEBUG
#CFLAGS += -g3 -Wall -check -pg
#CFLAGS += -pthread 
CFLAGS 	+= -fmax-errors=5
CFLAGS 	+= -fpermissive 
#CFLAGS += -Wno-deprecated #remove old style headers warning
#https://gcc.gnu.org/onlinedocs/gcc-10.1.0/gcc/x86-Options.html
#we have -mavx512f only on cluster
CFLAGS 	+= -march=skylake -DAVX
#CFLAGS 	+= -mavx512bw -mavx -DAVX avx512bw not here
LINK 	= $(CC) -o 
LDFLAGS =  -g
#LDFLAGS =  -g -rdynamic -pthread -lpthread


$(NAME):$(OBJS) 
	$(LINK) $(NAME) $(LDFLAGS) $(OBJS) $(LIBS)

#C++ source files

main.o:	main.cpp GPengine.h
	$(COMPILE) main.cpp

GPengine.o:	GPengine.cpp GPengine.h
	$(COMPILE) GPengine.cpp 

