# Makefile for IDA under LINUX:
#
# abbreviations :
#
FF = /opt/absoft10.2/bin/f77
LL = oba/l1.o oba/l2.o oba/l3.o oba/l4.o oba/l5.o \
oba/l6.o oba/l0x11.o
GG = oba/g1.o oba/g2.o
II = oba/i00.o oba/i01.o oba/i10.o oba/i20.o oba/i23.o oba/i25.o \
oba/i30.o oba/i32.o oba/i40.o oba/i41.o oba/i42.o oba/i43.o \
oba/i50.o oba/i60.o oba/i66.o oba/i67.o \
oba/i70.o oba/i71.o oba/i72.o oba/i73.o oba/i74.o oba/i75.o \
oba/i76.o oba/i80.o oba/i87.o oba/i95.o oba/i99.o \
oba/i77.o oba/i81.o oba/napi.o oba/napif.o
#
LD = l_dim.f
ID = i_dim.f
#
#CFLAGS = -O4 -ffixed-line-length-132 -assume bac -check nounderflow 
#CFLAGS = -O4 -132 -dusty -w -save
#CFLAGS = -O0 -g -132 -dusty -w -save
CFLAGL = -O1 -lm -Bstatic
# -g debugging symbols (excl. -O)
# -O optimize
# -A suppress alignment warnings
# -W wide format
# -s all storage static and initialized to zero
# -w suppress warnings
# -f fold symbols to lower case (for compatibility with libnag and libV77)
# -B108 append underscore to symbols (dito)
# missing: no underflow check

CB   = -check bounds 
# Compiling with libnag.a binary provided by LRZ
CFLAGS = -O  -A  -W -s -w -f -B108
LNAG = -static -m32 -L/usr/local/lib -lnag -lgfortran -lm 
LGET = -static -L/usr/local/lib -lget
#LNAG = -static -L/usr/local/lib -lnag -lgfortran -lm   
#LNEX is used for the Nexus read in
LNEX = -static -L/usr/local/lib -lz -ljpeg -L/usr/local/lib/hdf/lib -lhdf5 -lmfhdf -ldf 
LV77 = -L/opt/absoft10.2/lib -lV77 -lU77 -lNeXus
# -X options to be evaluated by ld
# -Bstatic option for ld to link all libraries statically
# libnag: Nag libraries
# libg2c: G77 routines references by libnag
# libm: math libraries
# libV77.a provides VAX-compatible DATE, IDATE and TIME
#
#
# Compile via C-stage
#   gcc -O -c -I/usr/local/lib/NAGWare -Wuninitialized
# link :
#
exa/frida1 : $(II) $(GG) $(LL)
	$(FF) -o exa/frida1 $(II) $(GG) $(LL) $(LNAG) $(LGET) $(LNEX) $(LV77) 
#
# compile :
#
#i_dim.f	    : i_dim.alp	    ; cp i_dim.alp i_dim.f
oba/i00.o   : i00.f $(ID)   ; $(FF) -c $(CFLAGS) i00.f -o oba/i00.o
oba/i01.o   : i01.f $(ID)   ; $(FF) -c $(CFLAGS) i01.f -o oba/i01.o
oba/i10.o   : i10.f $(ID)   ; $(FF) -c $(CFLAGS) i10.f -o oba/i10.o
oba/i20.o   : i20.f $(ID)   ; $(FF) -c $(CFLAGS) i20.f -o oba/i20.o
oba/i23.o   : i23.f $(ID)   ; $(FF) -c $(CFLAGS) i23.f -o oba/i23.o
oba/i25.o   : i25.f $(ID)   ; $(FF) -c $(CFLAGS) i25.f -o oba/i25.o
oba/i30.o   : i30.f $(ID)   ; $(FF) -c $(CFLAGS) i30.f -o oba/i30.o
oba/i32.o   : i32.f $(ID)   ; $(FF) -c $(CFLAGS) i32.f -o oba/i32.o
oba/i40.o   : i40.f $(ID)   ; $(FF) -c $(CFLAGS) i40.f -o oba/i40.o
oba/i41.o   : i41.f $(ID)   ; $(FF) -c $(CFLAGS) i41.f -o oba/i41.o
oba/i42.o   : i42.f $(ID)   ; $(FF) -c $(CFLAGS) i42.f -o oba/i42.o
oba/i43.o   : i43.f $(ID)   ; $(FF) -c $(CFLAGS) i43.f -o oba/i43.o
oba/i50.o   : i50.f $(ID)   ; $(FF) -c $(CFLAGS) i50.f -o oba/i50.o
oba/i60.o   : i60.f $(ID)   ; $(FF) -c $(CFLAGS) i60.f -o oba/i60.o
oba/i66.o   : i66.f $(ID)   ; $(FF) -c $(CFLAGS) i66.f -o oba/i66.o
oba/i67.o   : i67.f $(ID)   ; $(FF) -c $(CFLAGS) i67.f -o oba/i67.o
oba/i70.o   : i70.f $(ID)   ; $(FF) -c $(CFLAGS) i70.f -o oba/i70.o
oba/i71.o   : i71.f $(ID)   ; $(FF) -c $(CFLAGS) i71.f -o oba/i71.o
oba/i72.o   : i72.f $(ID)   ; $(FF) -c $(CFLAGS) i72.f -o oba/i72.o
oba/i73.o   : i73.f $(ID)   ; $(FF) -c $(CFLAGS) i73.f -o oba/i73.o
oba/i74.o   : i74.f $(ID)   ; $(FF) -c $(CFLAGS) i74.f -o oba/i74.o
oba/i75.o   : i75.f $(ID)   ; $(FF) -c $(CFLAGS) i75.f -o oba/i75.o
oba/i76.o   : i76.f $(ID)   ; $(FF) -c $(CFLAGS) i76.f -o oba/i76.o
oba/i77.o   : i77.f $(ID)   ; $(FF) -c $(CFLAGS) i77.f -o oba/i77.o
oba/i80.o   : i80.f $(ID)   ; $(FF) -c $(CFLAGS) i80.f -o oba/i80.o
oba/i81.o   : i81.f $(ID)   ; $(FF) -c $(CFLAGS) i81.f -o oba/i81.o
oba/i82.o   : i81.f $(ID)   ; $(FF) -c $(CFLAGS) i81.f -o oba/i82.o
oba/i87.o   : i87.f $(ID)   ; $(FF) -c $(CFLAGS) i87.f -o oba/i87.o
oba/i95.o   : i95.f $(ID)   ; $(FF) -c $(CFLAGS) i95.f -o oba/i95.o
oba/i99.o   : i99.f $(ID)   ; $(FF) -c $(CFLAGS) i99.f -o oba/i99.o
oba/l0x11.o : l0x11.f ; $(FF) -c $(CFLAGS) l0x11.f -o oba/l0x11.o
oba/l1.o    : l1.f    ; $(FF) -c $(CFLAGS) l1.f -o oba/l1.o
oba/l2.o    : l2.f    ; $(FF) -c $(CFLAGS) l2.f -o oba/l2.o
oba/l3.o    : l3.f $(LD)   ; $(FF) -c $(CFLAGS) l3.f -o oba/l3.o
oba/l4.o    : l4.f    ; $(FF) -c $(CFLAGS) l4.f -o oba/l4.o
oba/l5.o    : l5.f    ; $(FF) -c $(CFLAGS) l5.f -o oba/l5.o
oba/l6.o    : l6.f    ; $(FF) -c $(CFLAGS) l6.f -o oba/l6.o
oba/g1.o    : g1.f    ; $(FF) -c $(CFLAGS) g1.f -o oba/g1.o
oba/g2.o    : g2.f g_dim.f   ; $(FF) -c $(CFLAGS) g2.f -o oba/g2.o
#remove the following two lines if nexus read in is not used.
oba/napif.o : napif.f napif.inc; /opt/absoft10.2/bin/f77 -m32 -g -f -X -static -c napif.f -o oba/napif.o
#oba/napif.o : napif.f napif.inc; /opt/absoft10.2/bin/f77 -m32 -g -f -B108 -X -static -c napif.f -o oba/napif.o
oba/napi.o  : napi.c ; cc -m32  -static -c napi.c -o oba/napi.o

