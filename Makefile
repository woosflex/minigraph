CC=			gcc
CXX ?=			g++
CFLAGS=		-g -Wall -Wc++-compat -std=c99 -msse4 -O3
CPPFLAGS=
INCLUDES=
OBJS=		kalloc.o kthread.o algo.o sys.o gfa-base.o gfa-io.o gfa-aug.o gfa-bbl.o gfa-ed.o \
            sketch.o misc.o bseq.o options.o shortk.o miniwfa.o \
			index.o lchain.o gchain1.o galign.o gcmisc.o map-algo.o cal_cov.o \
			format.o gmap.o ggsimple.o ggen.o asm-call.o
PROG=		minigraph
LIBS=		-lz -lpthread -lm

# --- TracEon backend (optional; default off) ---
# TRACEON=1 backs minigraph's per-bucket khashl mg_hidx minimizer table with
# TracEon's traceon_kmer C API (libtraceon_kmer.a, see
# https://github.com/woosflex/TracEon). The final executable is linked via
# $(CXX) in TRACEON builds, which resolves libstdc++/libtraceon_kmer symbols
# (minigraph compiles as C; only mg_traceon_crc.cpp is C++). TRACEON_INC/
# TRACEON_LIB point at the TracEon include dir and build dir and are
# overridable on the command line. `make TRACEON=1` builds the traceon binary;
# plain `make` stays byte-for-byte on the stock khashl path.
#
# The same TRACEON_BACKEND build additionally gains the TG2 ".tgcache"
# zero-rebuild index cache (mg_traceon_cache.c + the CRC32C shim
# mg_traceon_crc.cpp): `minigraph -i ref.tgcache ref.gfa reads.fq` builds and
# dumps; `minigraph ref.tgcache reads.fq` mmap()s the index and points the
# buckets at the mapped slot arrays — zero table rebuild (see the format
# documentation at the top of mg_traceon_cache.c).
TRACEON ?= 0
TRACEON_INC ?= $(HOME)/agent_workspace/TracEon/include
TRACEON_LIB ?= $(HOME)/agent_workspace/TracEon/build
ifeq ($(TRACEON),1)
CPPFLAGS += -DTRACEON_BACKEND -I$(TRACEON_INC)
LIBS += -L$(TRACEON_LIB) -ltraceon_kmer
OBJS += mg_traceon_cache.o mg_traceon_crc.o
LINKER=$(CXX)
else
LINKER=$(CC)
endif

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .o
.PHONY:all clean depend FORCE

all:$(PROG)

# Rebuild every object whenever the compile flags change (TRACEON=1 vs stock,
# asan, ...): minigraph has no dependency tracking, so without this stamp a
# plain `make` after `make TRACEON=1` would silently relink stale objects.
BUILD_FLAGS=$(CFLAGS) $(CPPFLAGS) $(INCLUDES)
.build_flags: FORCE
	@echo '$(BUILD_FLAGS)' | cmp -s - $@ || echo '$(BUILD_FLAGS)' > $@
$(OBJS): .build_flags

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c -g -Wall -O2 -DTRACEON_HAS_AVX2 $(CPPFLAGS) $< -o $@

minigraph:$(OBJS) main.o
		$(LINKER) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

algo.o: kalloc.h algo.h miniwfa.h kvec-km.h ksort.h
asm-call.o: mgpriv.h minigraph.h gfa.h ggen.h bseq.h gfa-priv.h algo.h
bseq.o: bseq.h kvec-km.h kalloc.h kseq.h
cal_cov.o: mgpriv.h minigraph.h gfa.h gfa-priv.h algo.h kalloc.h
format.o: kalloc.h mgpriv.h minigraph.h gfa.h
galign.o: mgpriv.h minigraph.h gfa.h kalloc.h miniwfa.h
gchain1.o: mgpriv.h minigraph.h gfa.h ksort.h khashl.h kalloc.h gfa-priv.h
gcmisc.o: mgpriv.h minigraph.h gfa.h kalloc.h
gfa-aug.o: gfa-priv.h gfa.h ksort.h
gfa-base.o: gfa-priv.h gfa.h kstring.h khashl.h kalloc.h ksort.h
gfa-bbl.o: gfa-priv.h gfa.h kalloc.h ksort.h kvec.h
gfa-ed.o: gfa-priv.h gfa.h kalloc.h ksort.h khashl.h kdq.h kvec-km.h
gfa-io.o: kstring.h gfa-priv.h gfa.h kseq.h
ggen.o: kthread.h kalloc.h sys.h bseq.h ggen.h minigraph.h gfa.h mgpriv.h
ggen.o: gfa-priv.h
ggsimple.o: mgpriv.h minigraph.h gfa.h gfa-priv.h kalloc.h bseq.h algo.h
ggsimple.o: sys.h ggen.h kvec-km.h
gmap.o: kthread.h kalloc.h bseq.h sys.h mgpriv.h minigraph.h gfa.h gfa-priv.h
index.o: mgpriv.h minigraph.h gfa.h khashl.h kalloc.h kthread.h kvec-km.h
index.o: sys.h
kalloc.o: kalloc.h
kthread.o: kthread.h
lchain.o: mgpriv.h minigraph.h gfa.h kalloc.h krmq.h
main.o: mgpriv.h minigraph.h gfa.h gfa-priv.h sys.h ketopt.h
map-algo.o: kalloc.h mgpriv.h minigraph.h gfa.h khashl.h ksort.h
miniwfa.o: miniwfa.h kalloc.h
misc.o: mgpriv.h minigraph.h gfa.h ksort.h
options.o: mgpriv.h minigraph.h gfa.h sys.h
shortk.o: mgpriv.h minigraph.h gfa.h ksort.h kavl.h algo.h khashl.h kalloc.h
sketch.o: kvec-km.h kalloc.h mgpriv.h minigraph.h gfa.h
sys.o: sys.h
