BASEFLAGS = -Wall -g -Wshadow -Wstrict-prototypes -Wmissing-prototypes -static
# Edit the following to match your system:
BINDIR = $(HOME)/bin          # where to put the executables
CC = gcc                      # specify your c compiler

# Uncomment one of the following:
#E = .exe                     # suffix for executables under MSDOS
E =                           # (empty) suffix for executables under Unix

# Uncomment one of the following:
CFLAGS = $(BASEFLAGS)             # use this for debugging
#CFLAGS = $(BASEFLAGS) -DNDEBUG     # use this for production

# Leave the rest alone
.c.o:
	$(CC) $(CFLAGS) -c $<
TARGETS = mmest$(E) mmci$(E) mmgen$(E) normdist$(E) ci2ptx$(E) \
 ci2txt$(E) seqstat$(E) mmstr$(E) mmcistr$(E)

SRC = alloc2d.c bye.c charplot.c ci2ptx.c ci2txt.c estimate.c getcic.c \
getdata.c header.c identity.c iscoales.c lblaxis.c \
seqstat.c maxlike.c memstack.c misc.c mmci.c mmcistr.c mmest.c mmgen.c \
mmstr.c normdist.c readstr.c samedist.c simulate.c testh0.c \
unirand.c

all : $(TARGETS) 
	touch all

MMCIOBJ = mmci.o alloc2d.o bye.o estimate.o getcic.o getdata.o \
header.o identity.o iscoales.o lblaxis.o memstack.o \
misc.o simulate.o testh0.o unirand.o chol.o savetree.o eprintf.o
mmci$(E) : $(MMCIOBJ)
	$(CC) $(CFLAGS) -o $@ $(MMCIOBJ) -lm

MMESTOBJ = mmest.o alloc2d.o getdata.o misc.o estimate.o identity.o \
header.o memstack.o getcic.o bye.o
mmest$(E) : $(MMESTOBJ) 
	$(CC) $(CFLAGS) -o $@ $(MMESTOBJ) -lm

MMGENOBJ = mmgen.o iscoales.o unirand.o memstack.o charplot.o \
header.o identity.o misc.o simulate.o alloc2d.o savetree.o \
estimate.o getcic.o bye.o
mmgen$(E) : $(MMGENOBJ) 
	$(CC) $(CFLAGS) -o $@ $(MMGENOBJ) -lm

CI2TXTOBJ = ci2txt.o getdata.o misc.o lblaxis.o alloc2d.o \
getcic.o bye.o maxlike.o header.o
ci2txt$(E) : $(CI2TXTOBJ) 
	$(CC) $(CFLAGS) -o $@ $(CI2TXTOBJ) -lm

CI2PTXOBJ = ci2ptx.o getdata.o misc.o lblaxis.o alloc2d.o \
getcic.o bye.o maxlike.o
ci2ptx$(E) : $(CI2PTXOBJ) 
	$(CC) $(CFLAGS) -o $@ $(CI2PTXOBJ) -lm

NORMDISTOBJ = normdist.o
normdist$(E) : $(NORMDISTOBJ) 
	$(CC) $(CFLAGS) -o $@ $(NORMDISTOBJ) -lm

SEQSTATOBJ = seqstat.o misc.o bye.o getcic.o header.o
seqstat$(E) : $(SEQSTATOBJ) 
	$(CC) $(CFLAGS) -o $@ $(SEQSTATOBJ) -lm

MMSTROBJ = mmstr.o readstr.o alloc2d.o bye.o header.o unirand.o str_misc.o 
mmstr$(E) : $(MMSTROBJ) 
	$(CC) $(CFLAGS) -o $@ $(MMSTROBJ) -lm

MMCISTROBJ = mmcistr.o samedist.o alloc2d.o unirand.o iscoales.o memstack.o \
header.o bye.o misc.o readstr.o getcic.o str_misc.o
mmcistr$(E) : $(MMCISTROBJ) 
	$(CC) $(CFLAGS) -o $@ $(MMCISTROBJ) -lm

install : $(TARGETS)
	cp $(TARGETS) $(BINDIR)
	touch install

clean : 
	rm -f core junk* *.o *.zip *.tar *.tgz *~ $(TARGETS) all depend

zip : 
	make clean
	zip -k -l -j mmsrc.zip Makefile *.c *.h 

tgz :
	make clean
	tar zcf mmsrc.tgz Makefile *.c *.h

depend :  *.c *.h
	echo '#Automatically generated dependency info' > depend
	$(CC) -M *.c >> depend

include depend
