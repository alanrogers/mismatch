.c.o:
	$(CC) $(CFLAGS) -c -o ${<D}/${@F}  $<
#MMDIR is the top directory of the mismatch package on your system
MMDIR = $(HOME)/src/mismatch3
VPATH = $(MMDIR)/Misc
INCL = -I$(MMDIR)/Misc
CFLAGS = -DNDEBUG -O2 -Wall -g -DVERSION=$(VERSION) -DDATE=$(DATE) $(INCL)

OBJ = ci2ptx.o \
      $(MMDIR)/Misc/getdata.o \
      $(MMDIR)/Misc/misc.o \
      $(MMDIR)/Misc/lblaxis.o \
      $(MMDIR)/Misc/alloc2d.o \
      $(MMDIR)/Misc/ndx.o \
      $(MMDIR)/Misc/getcic.o \
      $(MMDIR)/Misc/bye.o \
      $(MMDIR)/Misc/maxlike.o
TARGET = ci2ptx

all : $(TARGET)

$(TARGET) : $(OBJ) depend
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

clean :
	rm -f core junk* *.o $(TARGET) depend

depend : ${OBJ:.o=.c}
	echo '#Automatically generated dependency info' > depend
	$(CC) $(INCL) -M ${OBJ:.o=.c} >> depend

include depend
