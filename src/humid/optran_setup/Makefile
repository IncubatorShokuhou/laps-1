# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../

include $(SRCROOT)/src/include/makefile.inc

#FFLAGS+=-pg

DBFFLAGS+=-pg

LAPSLIBS = 
#$(SRCROOT)/src/lib/liblaps.a


DBLAPSLIBS = 
#$(SRCROOT)/src/lib/libdebuglaps.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = binary_coeff.x

SRC= binary_coeff.f

MAIN =

LIB= libbinary.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)



$(EXE): $(ARCHIVE) $(MAIN)
	$(RANLIB) $(LIB)
	$(FC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i


distclean: clean

%.f : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f

%.o : %.F
























