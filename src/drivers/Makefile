vpath %.c ..
INCL = -I..
.c.o:
	$(CC) $(CFLAGS) $(INCL) -c -o ${@F}  $<
CFLAGS = -g -Wall -Wshadow -Wconversion -Wstrict-prototypes \
  -Wmissing-prototypes
CC = gcc
########################################################################
# xmmtheory
XMMTHEORY = xmmtheory.o mmtheory.o eprintf.o
xmmtheory : $(XMMTHEORY)
	$(CC) $(CFLAGS) -o $@ $(XMMTHEORY) -lm
########################################################################
# xeprintf
XEPRINTF = xeprintf.o eprintf.o
xeprintf : $(XEPRINTF)
	$(CC) $(CFLAGS) -o $@ $(XEPRINTF) -lm
########################################################################
# xgetf
XGETF = xgetf.o identity.o bye.o alloc2d.o misc.o
xgetf : $(XGETF)
	$(CC) $(CFLAGS) -o $@ $(XGETF) -lm
########################################################################
# xget_mse
XGET_MSE = xget_mse.o estimate.o misc.o ndx.o identity.o bye.o alloc2d.o
xget_mse : $(XGET_MSE)
	$(CC) $(CFLAGS) -o $@ $(XGET_MSE) -lm
########################################################################
clean :
	rm -f *.o core depend
########################################################################
depend :
	echo '#Automatically generated dependency info' > depend
	$(CC) -M $(INCL) *.c >> depend
########################################################################
include depend

