############################################################################
#
#       Project		Linux GPengine $Revision: 1.14 $ 
#
#       Author		W.B.Langdon
#
#	Created		19 June 2025 from gp/ei2022/lung.make r1.2
#
#Modifications:
#WBL 30 Jan 2026 comment out stats
#WBL 30 Dec 2025 restore stats (cf r1.9)
#WBL 11 Dec 2025 Add skylake-avx512
#WBL 5 Jul 2025 add AVX etc based on avx.make r1.9
############################################################################

NAME 	= GPengine
OBJS 	= main.o GPengine.o entropy.o
#LIBS 	= -lm
CC 	= g++
COMPILE = $(CC) -c $(CFLAGS)
CFLAGS 	= -O3
#CFLAGS += -Delvis
CFLAGS 	+= -g
#CFLAGS += -Dstats
CFLAGS += -DNDEBUG
#CFLAGS += -g3 -Wall -check -pg
CFLAGS  += -pthread 
CFLAGS 	+= -fmax-errors=5
CFLAGS 	+= -fpermissive 
#CFLAGS += -Wno-deprecated #remove old style headers warning
#https://gcc.gnu.org/onlinedocs/gcc-10.1.0/gcc/x86-Options.html
CFLAGS 	+= -march=skylake-avx512
LINK 	= $(CC) -o 
LDFLAGS =  -g  -pthread -lpthread
#LDFLAGS =  -g -rdynamic -pthread -lpthread


$(NAME):$(OBJS) 
	$(LINK) $(NAME) $(LDFLAGS) $(OBJS) $(LIBS)

#C++ source files

main.o:	main.cpp GPengine.h
	$(COMPILE) main.cpp

GPengine.o:	GPengine.cpp GPengine.h
	$(COMPILE) GPengine.cpp 

entropy.o:	entropy.cpp GPengine.h
	$(COMPILE) entropy.cpp
