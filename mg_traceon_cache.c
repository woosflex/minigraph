/*
 * mg_traceon_cache.c — "TG2" open-addressing index cache (TracEon backend)
 *
 * OPT-IN, additive cache mode for RECURRING runs, ported from minimap2's
 * TRC2 ".tcache" (mm_traceon_cache.c) to minigraph's khashl-based index.
 * minigraph has no .mmi-style dump/load path at all (it indexes GFA graphs,
 * not linear references), so this file adds one:
 *
 *   ./minigraph -i ref.tgcache ref.gfa reads.fq   # one-time build (slow path)
 *   ./minigraph ref.tgcache reads.fq              # recurring runs: ~zero-rebuild
 *
 * A .tgcache file stores (a) each bucket's minimizer entries as a khashl-style
 * open-addressing slot array + occupancy bitmap — mmap()ed and pointed at with
 * ZERO inserts on load — and (b) the graph reference block (segment names,
 * sequences, topology, and the derived edseq) so the graph mapper has the same
 * state it would after gfa_read()+mg_index(). Loading is: mmap + one
 * whole-file CRC32C check + pointer fixups + a heap copy of the (small) graph
 * reference block. No minimizer sketch, no radix sort, no table rebuild.
 *
 * Lookups are linear probes identical to khashl's: idx_hash(a)=(a)>>1 (bit 0
 * is the single-occurrence flag that hashing and equality ignore), probe start
 * __kh_h2b(hash, bits) = hash * 2654435769U >> (32-bits), advance while the
 * occupancy bit is set. This reproduces the exact khashl semantics of
 * index.c's mg_hidx table — including probe ORDER — so a .tgcache lookup and
 * a fresh khashl build return byte-identical (key, value) results.
 *
 * Only the mg_hidx table is replaced (by kmerindex_t* in the in-memory
 * TRACEON build, or by the mmap'd open-addressing array in tcache mode); the
 * per-bucket position array (mg_idx_bucket_t::p) and every other bucket field
 * stay exactly as upstream. minigraph's OTHER khashl maps (gfa-base's h_s2i,
 * shortk's sp/sp2, gchain1's kc, ...) are left untouched — stock khashl.
 *
 * ---------------------------------------------------------------------------
 * TG2 FILE FORMAT (little-endian, host byte order — same convention as TRC2)
 * ---------------------------------------------------------------------------
 *   Header (128 bytes):
 *     u32  magic     "TG2C"
 *     u16  version   1
 *     u16  flags     bit0 has_link_aux (always 1 for graphs dumped from the
 *                    full gfa_t; 0 would mean the link_aux blob is absent)
 *     i32  k, w, b   index parameters (b = bucket bits; n_bucket == 1<<b)
 *     u32  n_bucket  1<<b
 *     u32  n_seg     number of graph segments
 *     u32  n_sseq    number of persistent (stable) names
 *     i32  max_rank  max segment rank
 *     u64  n_arc     number of arcs (both directions present)
 *     u64  n_entry   total (key,value) entries across ALL buckets
 *     u64  p_total   total u64 count across ALL per-bucket position (p) arrays
 *     u64  sum_len   total graph sequence length in bases
 *     u64  names_len segment-names blob size ([u8 len][len bytes], no NUL)
 *     u64  seqs_len  segment-sequences blob size (concatenated, no NULs)
 *     u64  sseq_names_len stable-names blob size ([u8 len][len bytes])
 *     u64  edseq_len total bytes of the reverse-complement edseq blob
 *     u64  aux_len   total bytes of the link_aux blob (see below)
 *     u64  reserved  0
 *
 *   Graph reference block (restores gfa_t + gfa_edseq_t WITHOUT parsing):
 *     names blob   names_len bytes: per segment [u8 len][len bytes]
 *     (pad to 8)
 *     lens array   n_seg * u32 segment lengths
 *     (pad to 8)
 *     meta array   n_seg * 16 bytes: [u32 del:16,circ:16][i32 snid][i32 soff]
 *                  [i32 rank]
 *     (pad to 8)
 *     seqs blob    seqs_len bytes: per segment the raw sequence (no NULs);
 *                  offsets follow from the cumulative lens
 *     (pad to 8)
 *     sseq names   sseq_names_len bytes: per stable name [u8 len][len bytes]
 *     (pad to 8)
 *     sseq array   n_sseq * 12 bytes: [i32 min][i32 max][i32 rank]
 *     (pad to 8)
 *     arc array    n_arc * 32 raw bytes of gfa_arc_t (host layout, memcpy'd;
 *                  both directions present, exactly as gfa_finalize left it)
 *     (pad to 8)
 *     idx array    n_vtx * u64 (gfa_t::idx, vertex arc index; n_vtx==2*n_seg)
 *     link_aux     aux_len bytes: per link k: [u32 l_aux][l_aux aux bytes]
 *                  (l_aux==0 entries contribute 4 zero bytes; needed by
 *                  --cov/gfa_print; the mapping path never reads it)
 *     (pad to 8)
 *     edseq blob   edseq_len bytes: per segment the reverse complement;
 *                  lengths == segment lens, offsets follow from the lens
 *
 *   Table block (zero-rebuild, mmap-pointable):
 *     p_off   (n_bucket+1) * u64 cumulative counts of p-array entries
 *             (p_off[0]=0; bucket i's p array = p_blob[p_off[i]..p_off[i+1]))
 *     b_off   (n_bucket+1) * u64 cumulative BYTE offsets into the bucket blob
 *             (b_off[0]=0; bucket i's table = b_blob[b_off[i]..b_off[i+1]));
 *             the total bucket-blob size is b_off[n_bucket]
 *     b_blob  per-bucket open-addressing tables, concatenated; empty buckets
 *             contribute ZERO bytes (b_off[i+1] == b_off[i]). Each non-empty
 *             table is, in order (8-byte aligned):
 *               u32  capacity  power of two, >= 4; load factor (n/capacity)
 *                              <= 0.75 at dump time (same target as khashl's
 *                              putp rehash threshold)
 *               u32  bits      log2(capacity)
 *               slots  capacity * 16 bytes: (u64 key, u64 value) pairs,
 *                              placed by khashl's linear probe; empty slots
 *                              are ANY bytes (written as 0); occupancy lives
 *                              ONLY in the bitmap
 *               used   ceil(capacity/32) * 4 bytes: khashl-style occupancy
 *                              bitmap (bit i in word i>>5, position i&31)
 *               (pad to 8)
 *     p_blob  p_total * u64 concatenated per-bucket position arrays (each
 *             bucket's already sorted by position, as in the khashl build)
 *
 *   Trailer:
 *     u32  crc32c  CRC-32C (Castagnoli, init 0xFFFFFFFF, final XOR) over the
 *             WHOLE file from byte 0 up to (excluding) this trailer — the
 *             TracEon .traceon v4 whole-payload integrity pattern. Computed
 *             with TracEon's include/Crc32c.h (see mg_traceon_crc.cpp).
 *
 * The layout is computed by tg_calc_layout() — the ONLY place the offset math
 * lives — used by both save and load so they can never drift apart.
 *
 * COVERAGE: the format covers the plain mapping flow (GAF/PAF output), the
 * --cov coverage path, and gfa_print. It does NOT cover the incremental graph
 * generation (--ggen) flow, which augments the graph and rebuilds the index
 * in-process; that flow never loads a .tgcache. Re-dumping from an already
 * loaded .tgcache (minigraph -i out.tgcache in.tgcache reads.fq) is supported.
 *
 * MIT License, copyright (c) 2026 — TracEon integration (see LICENSE.txt).
 */

#ifdef TRACEON_BACKEND

#define _POSIX_C_SOURCE 200809L /* fileno/ftello/fseeko for the tcache loader */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <assert.h>

#include "minigraph.h"
#include "mgpriv.h"
#include "gfa-priv.h"
#include "kalloc.h" // kcalloc/kmalloc (NULL km falls back to malloc/calloc)
#include "sys.h"    // realtime()
#include "kmerindex_c_api.h"

/* mg_idx_init lives in index.c without a header declaration (it is an internal
 * entry point); declare it here for the loader. */
extern mg_idx_t *mg_idx_init(int k, int w, int b);

/* ---- CRC-32C (Castagnoli) — C shim over TracEon's include/Crc32c.h, which is
 * ---- header-only C++; the thin wrapper lives in mg_traceon_crc.cpp and is
 * ---- linked via $(CXX) in TRACEON builds (same link model as libtraceon_kmer). */
typedef struct mg_crc32c_s mg_crc32c_s;
mg_crc32c_s *mg_crc32c_new(void);
void mg_crc32c_free(mg_crc32c_s *s);
void mg_crc32c_update(mg_crc32c_s *s, const void *data, size_t len);
uint32_t mg_crc32c_final(mg_crc32c_s *s);
uint32_t mg_crc32c(const void *data, size_t len);

#define TG_MAGIC       "TG2C"
#define TG_VERSION     1
#define TG_HDR_SIZE    128

#define TG_F_HAS_LINK_AUX 0x1

#define TG_OFF_MAGIC    0
#define TG_OFF_VERSION  4
#define TG_OFF_FLAGS    6
#define TG_OFF_K        8
#define TG_OFF_W        12
#define TG_OFF_B        16
#define TG_OFF_NBUCKET  20
#define TG_OFF_NSEG     24
#define TG_OFF_NSSEQ    28
#define TG_OFF_MAXRANK  32
#define TG_OFF_NARC     40
#define TG_OFF_NENTRY   48
#define TG_OFF_PTOTAL   56
#define TG_OFF_SUMLEN   64
#define TG_OFF_NAMESLEN 72
#define TG_OFF_SEQSLEN  80
#define TG_OFF_SSEQNAMESLEN 88
#define TG_OFF_EDSEQLEN 96
#define TG_OFF_AUXLEN   104
#define TG_OFF_RSV      112

/* Byte layout of everything after the header (offsets relative to file base). */
typedef struct {
	uint64_t names_off;      // segment names blob
	uint64_t lens_off;       // n_seg * u32 lengths (8-aligned)
	uint64_t meta_off;       // n_seg * 16 (del/circ, snid, soff, rank) (8-aligned)
	uint64_t seqs_off;       // segment sequences blob (8-aligned)
	uint64_t sseq_names_off; // stable names blob (8-aligned)
	uint64_t sseq_off;       // n_sseq * 12 (min, max, rank) (8-aligned)
	uint64_t arc_off;        // n_arc * 32 raw gfa_arc_t (8-aligned)
	uint64_t idx_off;        // n_vtx * 8 (gfa_t::idx)
	uint64_t aux_off;        // link_aux blob (8-aligned)
	uint64_t edseq_off;      // reverse-complement edseq blob (8-aligned)
	uint64_t poff_off;       // (n_bucket+1) * u64 cumulative p counts
	uint64_t boff_off;       // (n_bucket+1) * u64 cumulative bucket-blob byte offsets
	uint64_t bblob_off;      // bucket tables blob
	uint64_t bblob_size;     // bytes of the bucket tables blob
	uint64_t pblob_off;      // p_total * 8 bytes
	uint64_t payload_size;   // == file size minus the 4-byte CRC trailer
} tg_layout_t;

static uint64_t tg_round8(uint64_t x) { return (x + 7) & ~(uint64_t)7; }

static void tg_calc_layout(uint64_t sum_len, uint32_t n_seg, uint32_t n_sseq,
                           uint32_t n_bucket, uint64_t n_arc, uint64_t names_len,
                           uint64_t seqs_len, uint64_t sseq_names_len,
                           uint64_t edseq_len, uint64_t aux_len,
                           uint64_t bblob_size, uint64_t p_total, tg_layout_t *L)
{
	uint64_t n_vtx = (uint64_t)n_seg << 1;
	L->names_off = TG_HDR_SIZE;
	L->lens_off = tg_round8(L->names_off + names_len);
	L->meta_off = tg_round8(L->lens_off + (uint64_t)n_seg * 4);
	L->seqs_off = tg_round8(L->meta_off + (uint64_t)n_seg * 16);
	L->sseq_names_off = tg_round8(L->seqs_off + seqs_len);
	L->sseq_off = tg_round8(L->sseq_names_off + sseq_names_len);
	L->arc_off = tg_round8(L->sseq_off + (uint64_t)n_sseq * 12);
	L->idx_off = tg_round8(L->arc_off + n_arc * 32);
	L->aux_off = L->idx_off + n_vtx * 8;
	L->edseq_off = tg_round8(L->aux_off + aux_len);
	L->poff_off = tg_round8(L->edseq_off + edseq_len);
	L->boff_off = L->poff_off + (uint64_t)(n_bucket + 1) * 8;
	L->bblob_off = L->boff_off + (uint64_t)(n_bucket + 1) * 8;
	L->bblob_size = bblob_size;
	L->pblob_off = L->bblob_off + bblob_size;
	L->payload_size = L->pblob_off + p_total * 8;
}

static uint16_t tg_rd16(const uint8_t *p) { uint16_t v; memcpy(&v, p, 2); return v; }
static uint32_t tg_rd32(const uint8_t *p) { uint32_t v; memcpy(&v, p, 4); return v; }
static uint64_t tg_rd64(const uint8_t *p) { uint64_t v; memcpy(&v, p, 8); return v; }
static void tg_wr16(uint8_t *p, uint16_t v) { memcpy(p, &v, 2); }
static void tg_wr32(uint8_t *p, uint32_t v) { memcpy(p, &v, 4); }
static void tg_wr64(uint8_t *p, uint64_t v) { memcpy(p, &v, 8); }

/* Open-addressing table capacity for n entries: smallest power of two
 * >= ceil(n / 0.75), minimum 4 — mirrors khashl's putp rehash threshold
 * (count >= (n_buckets>>1)+(n_buckets>>2) triggers resize), so the load factor
 * at dump time is <= 0.75 and a probe never scans a full table. */
static uint32_t tg_table_capacity(uint64_t n)
{
	uint64_t need = (n * 4 + 2) / 3;
	uint32_t cap = 4;
	while (cap < need && cap <= 0x40000000u) cap <<= 1;
	assert(cap >= need);
	return cap;
}

/* khashl's __kh_h2b: probe start = (uint32)hash * 2654435769U >> (32-bits),
 * where hash = idx_hash(key) = key>>1 (truncated to 32 bits, as khashl does). */
static inline uint32_t tg_h2b(uint64_t key, uint32_t bits)
{
	return (uint32_t)((uint32_t)(key >> 1) * 2654435769U) >> (32 - bits);
}

/* Probe-insert (key,value) into an open-addressing table using khashl's exact
 * probe order. The stored key keeps its bit0 (the singleton flag); only key>>1
 * participates in hashing and comparison — identical to khashl's idx_hash /
 * idx_eq semantics. */
static void tg_probe_insert(uint64_t *slots, uint32_t *used, uint32_t cap, uint32_t bits,
                            uint64_t k, uint64_t v)
{
	uint32_t mask = cap - 1;
	uint32_t idx = tg_h2b(k, bits);
	while (used[idx >> 5] & (1u << (idx & 31))) idx = (idx + 1) & mask;
	slots[(size_t)idx << 1] = k;
	slots[((size_t)idx << 1) + 1] = v;
	used[idx >> 5] |= (1u << (idx & 31));
}

/* Popcount the occupancy bitmap of a table with `cap` slots (khashl layout:
 * u32 words, bit i in word i>>5). Only bits [0, cap) belong to this table. */
#ifdef _MSC_VER
#include <intrin.h>
#define TG_POPCNT64(x) __popcnt64(x)
#else
#define TG_POPCNT64(x) __builtin_popcountll(x)
#endif
static uint32_t tg_bitmap_popcount(const uint32_t *used, uint32_t cap)
{
	uint32_t n = 0, nw = (cap + 31) >> 5, i;
	for (i = 0; i < nw; ++i)
		n += (uint32_t)TG_POPCNT64((uint64_t)used[i]);
	{
		uint32_t rem = cap & 31;
		if (rem) n -= (uint32_t)TG_POPCNT64((uint64_t)used[nw - 1] & ~((1u << rem) - 1));
	}
	return n;
}

/* Write `n` bytes from `buf` and feed them to the running CRC. */
static int tg_write(FILE *fp, mg_crc32c_s *crc, const void *buf, size_t n)
{
	if (fwrite(buf, 1, n, fp) != n) return -1;
	mg_crc32c_update(crc, buf, n);
	return 0;
}

int mg_tcache_dump(FILE *fp, const mg_idx_t *gi)
{
	uint64_t i, j, k, sum_len = 0, names_len = 0, seqs_len = 0, sseq_names_len = 0;
	uint64_t edseq_len = 0, aux_len = 0, n_entry = 0, p_total = 0, bblob_size = 0;
	uint32_t n_bucket = 1U << gi->b, has_link_aux = 1;
	uint64_t *ecount, *pcount, *bsize;
	uint32_t *cap, *bits;
	tg_layout_t L;
	mg_crc32c_s *crc = 0;
	uint8_t hdr[TG_HDR_SIZE];
	const gfa_t *g = gi->g;
	uint64_t n_vtx;

	if (fp == 0 || gi == 0 || g == 0) return -1;

	n_vtx = (uint64_t)g->n_seg << 1;
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		sum_len += s->len;
		seqs_len += s->len;
		if (s->name) names_len += 1 + strlen(s->name);
		else names_len += 1;
	}
	for (i = 0; i < g->n_sseq; ++i)
		if (g->sseq[i].name) sseq_names_len += 1 + strlen(g->sseq[i].name);
	for (i = 0; i < gi->n_seg; ++i) // edseq: reverse-complement strings
		edseq_len += gi->es[i<<1|1].len;


	/* pass A: sizes; per-bucket table capacity + byte size follow from the
	 * entry count (load factor target <= 0.75) */
	ecount = (uint64_t*)calloc(n_bucket, 8);
	pcount = (uint64_t*)calloc(n_bucket, 8);
	bsize  = (uint64_t*)calloc(n_bucket, 8);
	cap    = (uint32_t*)calloc(n_bucket, 4);
	bits   = (uint32_t*)calloc(n_bucket, 4);
	if (ecount == 0 || pcount == 0 || bsize == 0 || cap == 0 || bits == 0) goto fail0;
	for (i = 0; i < n_bucket; ++i) {
		const mg_idx_bucket_t *b = &gi->B[i];
		ecount[i] = gi->is_tcache? (uint64_t)b->ne
		                        : (b->h? kmerindex_size((const kmerindex_t*)b->h) : 0);
		pcount[i] = (uint64_t)b->n;
		n_entry += ecount[i];
		p_total += pcount[i];
		if (ecount[i]) {
			cap[i] = tg_table_capacity(ecount[i]);
			bits[i] = 31 - (uint32_t)__builtin_clz(cap[i]);
			bsize[i] = tg_round8(8 + (uint64_t)cap[i] * 16 + ((uint64_t)cap[i] + 31) / 32 * 4);
			bblob_size += bsize[i];
		}
	}
	/* link_aux blob: per link k: u32 l_aux + l_aux bytes */
	if (g->link_aux) {
		for (k = 0; k < g->n_arc; ++k)
			aux_len += 4 + g->link_aux[k].l_aux;
	} else has_link_aux = 0;
	tg_calc_layout(sum_len, g->n_seg, g->n_sseq, n_bucket, g->n_arc, names_len,
	               seqs_len, sseq_names_len, edseq_len, aux_len, bblob_size, p_total, &L);

	/* header */
	memset(hdr, 0, sizeof(hdr));
	memcpy(hdr + TG_OFF_MAGIC, TG_MAGIC, 4);
	tg_wr16(hdr + TG_OFF_VERSION, TG_VERSION);
	tg_wr16(hdr + TG_OFF_FLAGS, has_link_aux? TG_F_HAS_LINK_AUX : 0);
	tg_wr32(hdr + TG_OFF_K, (uint32_t)gi->k);
	tg_wr32(hdr + TG_OFF_W, (uint32_t)gi->w);
	tg_wr32(hdr + TG_OFF_B, (uint32_t)gi->b);
	tg_wr32(hdr + TG_OFF_NBUCKET, n_bucket);
	tg_wr32(hdr + TG_OFF_NSEG, g->n_seg);
	tg_wr32(hdr + TG_OFF_NSSEQ, g->n_sseq);
	tg_wr32(hdr + TG_OFF_MAXRANK, g->max_rank);
	tg_wr64(hdr + TG_OFF_NARC, g->n_arc);
	tg_wr64(hdr + TG_OFF_NENTRY, n_entry);
	tg_wr64(hdr + TG_OFF_PTOTAL, p_total);
	tg_wr64(hdr + TG_OFF_SUMLEN, sum_len);
	tg_wr64(hdr + TG_OFF_NAMESLEN, names_len);
	tg_wr64(hdr + TG_OFF_SEQSLEN, seqs_len);
	tg_wr64(hdr + TG_OFF_SSEQNAMESLEN, sseq_names_len);
	tg_wr64(hdr + TG_OFF_EDSEQLEN, edseq_len);
	tg_wr64(hdr + TG_OFF_AUXLEN, aux_len);

	crc = mg_crc32c_new();
	if (crc == 0) goto fail0;
	if (tg_write(fp, crc, hdr, sizeof(hdr)) != 0) goto fail;

	/* names blob: [u8 len][len bytes] per segment */
	for (i = 0; i < g->n_seg; ++i) {
		const char *name = g->seg[i].name;
		uint8_t nb[256];
		size_t l = name? strlen(name) : 0;
		assert(l < 256);
		nb[0] = (uint8_t)l;
		if (l) memcpy(nb + 1, name, l);
		if (tg_write(fp, crc, nb, l + 1) != 0) goto fail;
	}
	/* pad names blob to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.lens_off - (TG_HDR_SIZE + names_len);
		if (pad && tg_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* lens array */
	{
		uint32_t *lens = (uint32_t*)malloc(g->n_seg * 4);
		if (lens == 0) goto fail;
		for (i = 0; i < g->n_seg; ++i) lens[i] = g->seg[i].len;
		if (tg_write(fp, crc, lens, g->n_seg * 4) != 0) { free(lens); goto fail; }
		free(lens);
	}
	/* pad lens to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.meta_off - L.lens_off - (uint64_t)g->n_seg * 4;
		if (pad && tg_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* meta array: [u32 del:16,circ:16][i32 snid][i32 soff][i32 rank] */
	{
		uint32_t *meta = (uint32_t*)malloc(g->n_seg * 16);
		if (meta == 0) goto fail;
		for (i = 0; i < g->n_seg; ++i) {
			const gfa_seg_t *s = &g->seg[i];
			meta[i*4 + 0] = (uint32_t)((uint32_t)s->del | ((uint32_t)s->circ << 16));
			meta[i*4 + 1] = (uint32_t)s->snid;
			meta[i*4 + 2] = (uint32_t)s->soff;
			meta[i*4 + 3] = (uint32_t)s->rank;
		}
		if (tg_write(fp, crc, meta, g->n_seg * 16) != 0) { free(meta); goto fail; }
		free(meta);
	}
	/* pad meta to 8 (n_seg*16 is already a multiple of 8 — no-op) */
	/* seqs blob */
	{
		char *seqs = (char*)malloc(seqs_len? seqs_len : 1);
		char *w = seqs;
		if (seqs == 0) goto fail;
		for (i = 0; i < g->n_seg; ++i) {
			const gfa_seg_t *s = &g->seg[i];
			if (s->len) { memcpy(w, s->seq, s->len); w += s->len; }
		}
		if (tg_write(fp, crc, seqs, seqs_len) != 0) { free(seqs); goto fail; }
		free(seqs);
	}
	/* pad seqs to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.sseq_names_off - (L.seqs_off + seqs_len);
		if (pad && tg_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* sseq names blob */
	for (i = 0; i < g->n_sseq; ++i) {
		const char *name = g->sseq[i].name;
		uint8_t nb[256];
		size_t l = name? strlen(name) : 0;
		assert(l < 256);
		nb[0] = (uint8_t)l;
		if (l) memcpy(nb + 1, name, l);
		if (tg_write(fp, crc, nb, l + 1) != 0) goto fail;
	}
	/* pad sseq names to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.sseq_off - (L.sseq_names_off + sseq_names_len);
		if (pad && tg_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* sseq array: [i32 min][i32 max][i32 rank] */
	if (g->n_sseq) {
		uint32_t *ss = (uint32_t*)malloc(g->n_sseq * 12);
		if (ss == 0) goto fail;
		for (i = 0; i < g->n_sseq; ++i) {
			ss[i*3 + 0] = (uint32_t)g->sseq[i].min;
			ss[i*3 + 1] = (uint32_t)g->sseq[i].max;
			ss[i*3 + 2] = (uint32_t)g->sseq[i].rank;
		}
		if (tg_write(fp, crc, ss, g->n_sseq * 12) != 0) { free(ss); goto fail; }
		free(ss);
	}
	/* pad sseq to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.arc_off - (L.sseq_off + (uint64_t)g->n_sseq * 12);
		if (pad && tg_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* arc array: raw gfa_arc_t */
	assert(sizeof(gfa_arc_t) == 32);
	if (g->n_arc && fwrite(g->arc, 32, g->n_arc, fp) != g->n_arc) goto fail;
	if (g->n_arc) mg_crc32c_update(crc, g->arc, g->n_arc * 32);
	/* pad arcs to 8 (n_arc*32 is already a multiple of 8 — no-op) */
	/* idx array */
	if (n_vtx && fwrite(g->idx, 8, n_vtx, fp) != n_vtx) goto fail;
	if (n_vtx) mg_crc32c_update(crc, g->idx, n_vtx * 8);
	/* link_aux blob */
	if (has_link_aux) {
		for (k = 0; k < g->n_arc; ++k) {
			const gfa_aux_t *a = &g->link_aux[k];
			uint32_t l = a->l_aux;
			if (tg_write(fp, crc, &l, 4) != 0) goto fail;
			if (l && tg_write(fp, crc, a->aux, l) != 0) goto fail;
		}
	}
	/* pad link_aux to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.edseq_off - (L.aux_off + aux_len);
		if (pad && tg_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* edseq blob: reverse complements */
	{
		char *buf = (char*)malloc(edseq_len? edseq_len : 1);
		char *w = buf;
		if (buf == 0) goto fail;
		for (i = 0; i < gi->n_seg; ++i) {
			const char *s = gi->es[i<<1|1].seq;
			int32_t l = gi->es[i<<1|1].len;
			if (l) { memcpy(w, s, l); w += l; }
		}
		if (tg_write(fp, crc, buf, edseq_len) != 0) { free(buf); goto fail; }
		free(buf);
	}
	/* pad edseq to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.poff_off - (L.edseq_off + edseq_len);
		if (pad && tg_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* p_off table: cumulative p-array entry counts */
	for (i = 0, j = 0; i <= n_bucket; ++i) {
		uint64_t v = j;
		if (tg_write(fp, crc, &v, 8) != 0) goto fail;
		if (i < n_bucket) j += pcount[i];
	}
	/* b_off table: cumulative byte offsets into the bucket blob */
	for (i = 0, j = 0; i <= n_bucket; ++i) {
		uint64_t v = j;
		if (tg_write(fp, crc, &v, 8) != 0) goto fail;
		if (i < n_bucket) j += bsize[i];
	}
	/* bucket blob: per bucket, an open-addressing table
	 * (u32 capacity, u32 bits, slots[capacity*16], used[(capacity+31)/32*4]) */
	for (i = 0; i < n_bucket; ++i) {
		uint64_t c = ecount[i];
		mg_idx_bucket_t *b = &gi->B[i];
		uint64_t *slots;
		uint32_t *used, capbuf[2];
		if (c == 0) continue;
		slots = (uint64_t*)calloc(cap[i], 16);
		used = (uint32_t*)calloc(((size_t)cap[i] + 31) >> 5, 4);
		if (slots == 0 || used == 0) { free(slots); free(used); goto fail; }
		if (gi->is_tcache) {
			// v2 source (re-dump): iterate the occupied slots of the source table
			uint32_t k2;
			const uint32_t *sbm = (const uint32_t*)b->bm; // khashl u32-word bitmap
			for (k2 = 0; k2 < b->cap; ++k2)
				if (sbm[k2 >> 5] & (1u << (k2 & 31)))
					tg_probe_insert(slots, used, cap[i], bits[i], b->fe[(size_t)k2 << 1], b->fe[((size_t)k2 << 1) + 1]);
		} else {
			// traceon table source: insert in the table's iteration order
			kmerindex_iter_t it;
			uint64_t kk, vv;
			uint32_t n = 0;
			kmerindex_iter_begin((const kmerindex_t*)b->h, &it);
			while (kmerindex_iter_next(&it, &kk, &vv)) {
				tg_probe_insert(slots, used, cap[i], bits[i], kk, vv);
				++n;
			}
			assert(n == c);
		}
		capbuf[0] = cap[i];
		capbuf[1] = bits[i];
		if (tg_write(fp, crc, capbuf, 8) != 0) { free(slots); free(used); goto fail; }
		if (tg_write(fp, crc, slots, (size_t)cap[i] * 16) != 0) { free(slots); free(used); goto fail; }
		if (tg_write(fp, crc, used, ((size_t)cap[i] + 31) / 32 * 4) != 0) { free(slots); free(used); goto fail; }
		{ // pad this bucket's table to 8 (bsize[] in pass A is round8'd)
			uint8_t zb[8] = {0,0,0,0,0,0,0,0};
			uint64_t pad = bsize[i] - (8 + (uint64_t)cap[i] * 16 + ((uint64_t)cap[i] + 31) / 32 * 4);
			if (pad && tg_write(fp, crc, zb, (size_t)pad) != 0) { free(slots); free(used); goto fail; }
		}
		free(slots); free(used);
	}
	/* p_blob: concatenated per-bucket position arrays (already sorted) */
	for (i = 0; i < n_bucket; ++i) {
		if (pcount[i] == 0) continue;
		if (fwrite(gi->B[i].p, 8, pcount[i], fp) != pcount[i]) goto fail;
		mg_crc32c_update(crc, gi->B[i].p, pcount[i] * 8);
	}
	/* whole-file CRC32C trailer */
	{
		uint32_t c = mg_crc32c_final(crc);
		if (fwrite(&c, 4, 1, fp) != 1) goto fail;
	}
	/* Flush and verify that the trailer is the LAST 4 bytes of the PHYSICAL
	 * file (a short/deferred write would otherwise surface as a cryptic CRC
	 * mismatch on load). */
	if (fflush(fp) != 0 || ferror(fp)) goto fail;
	{
		struct stat st;
		if (fstat(fileno(fp), &st) != 0 || (uint64_t)st.st_size != L.payload_size + 4) {
			fprintf(stderr, "[ERROR] mg_tcache_dump: write truncated (%lld of %llu bytes on disk) — the tcache file is invalid; delete it and retry with free space\n",
				(long long)st.st_size, (unsigned long long)(L.payload_size + 4));
			goto fail;
		}
	}
	mg_crc32c_free(crc);
	free(ecount); free(pcount); free(bsize); free(cap); free(bits);
	return 0;

fail:
	mg_crc32c_free(crc);
fail0:
	free(ecount); free(pcount); free(bsize); free(cap); free(bits);
	return -1;
}

static gfa_t *mg_tcache_load_graph(const uint8_t *map, const tg_layout_t *L,
                                   uint32_t n_seg, uint32_t n_sseq, uint32_t max_rank,
                                   uint64_t n_arc, uint64_t names_len, uint64_t seqs_len,
                                   uint64_t sseq_names_len, uint64_t aux_len,
                                   uint64_t edseq_len, int has_link_aux, gfa_edseq_t **es_)
{
	gfa_t *g;
	gfa_edseq_t *es;
	uint32_t i;
	uint64_t k, off;
	const uint8_t *np = map + L->names_off;
	const uint32_t *lens = (const uint32_t*)(map + L->lens_off);
	const uint32_t *meta = (const uint32_t*)(map + L->meta_off);
	const uint8_t *seqs = map + L->seqs_off;
	const uint8_t *snp = map + L->sseq_names_off;
	const uint32_t *ssarr = (const uint32_t*)(map + L->sseq_off);
	const gfa_arc_t *arcs = (const gfa_arc_t*)(map + L->arc_off);
	const uint64_t *idx = (const uint64_t*)(map + L->idx_off);
	const uint8_t *aux = map + L->aux_off;
	const uint8_t *ed = map + L->edseq_off;

	g = gfa_init();
	if (g == 0) return 0;
	g->n_seg = g->m_seg = n_seg;
	g->n_sseq = g->m_sseq = n_sseq;
	g->max_rank = max_rank;
	g->seg = (gfa_seg_t*)calloc(n_seg, sizeof(gfa_seg_t));
	g->sseq = (gfa_sseq_t*)calloc(n_sseq? n_sseq : 1, sizeof(gfa_sseq_t));
	if (g->seg == 0 || g->sseq == 0) goto fail;
	/* segments */
	off = 0;
	for (i = 0; i < n_seg; ++i) {
		gfa_seg_t *s = &g->seg[i];
		uint8_t l = np[0];
		uint32_t len = lens[i];
		if (l) {
			s->name = (char*)malloc(l + 1);
			if (s->name == 0) goto fail;
			memcpy(s->name, np + 1, l);
			s->name[l] = 0;
		} else s->name = 0;
		np += 1 + l;
		s->len = (int32_t)len;
		s->del = meta[i*4 + 0] & 0xffff;
		s->circ = (meta[i*4 + 0] >> 16) & 0xffff;
		s->snid = (int32_t)meta[i*4 + 1];
		s->soff = (int32_t)meta[i*4 + 2];
		s->rank = (int32_t)meta[i*4 + 3];
		if (len) {
			s->seq = (char*)malloc(len);
			if (s->seq == 0) goto fail;
			memcpy(s->seq, seqs + off, len);
		} else s->seq = 0;
		off += len;
	}
	assert(off == seqs_len);
	/* sseq */
	off = 0;
	for (i = 0; i < n_sseq; ++i) {
		gfa_sseq_t *s = &g->sseq[i];
		uint8_t l = snp[0];
		if (l) {
			s->name = (char*)malloc(l + 1);
			if (s->name == 0) goto fail;
			memcpy(s->name, snp + 1, l);
			s->name[l] = 0;
		} else s->name = 0;
		snp += 1 + l;
		s->min = (int32_t)ssarr[i*3 + 0];
		s->max = (int32_t)ssarr[i*3 + 1];
		s->rank = (int32_t)ssarr[i*3 + 2];
	}
	/* arcs + idx: heap copies (gfa_destroy free()s them) */
	if (n_arc) {
		g->arc = (gfa_arc_t*)malloc(n_arc * sizeof(gfa_arc_t));
		if (g->arc == 0) goto fail;
		memcpy(g->arc, arcs, n_arc * sizeof(gfa_arc_t));
	}
	g->m_arc = g->n_arc = n_arc;
	if (n_seg) {
		uint64_t n_vtx = (uint64_t)n_seg << 1;
		g->idx = (uint64_t*)malloc(n_vtx * 8);
		if (g->idx == 0) goto fail;
		memcpy(g->idx, idx, n_vtx * 8);
	}
	/* link_aux */
	if (has_link_aux && n_arc) {
		g->link_aux = (gfa_aux_t*)calloc(n_arc, sizeof(gfa_aux_t));
		if (g->link_aux == 0) goto fail;
		for (k = 0; k < n_arc; ++k) {
			uint32_t l = tg_rd32(aux);
			aux += 4;
			if (l) {
				g->link_aux[k].aux = (uint8_t*)malloc(l);
				if (g->link_aux[k].aux == 0) goto fail;
				memcpy(g->link_aux[k].aux, aux, l);
				aux += l;
			}
			g->link_aux[k].l_aux = l;
			g->link_aux[k].m_aux = l;
		}
	} else if (n_arc) {
		g->link_aux = (gfa_aux_t*)calloc(n_arc, sizeof(gfa_aux_t));
	}
	/* edseq: forward = segment sequence (heap copy in the graph), reverse =
	 * heap copies of the blob strings (gfa_edseq_destroy free()s them) */
	{
		uint32_t n_vtx = n_seg << 1;
		es = (gfa_edseq_t*)calloc(n_vtx? n_vtx : 1, sizeof(gfa_edseq_t));
		if (es == 0) goto fail;
		off = 0;
		for (i = 0; i < n_seg; ++i) {
			int32_t len = g->seg[i].len;
			es[i<<1].seq = g->seg[i].seq; // forward: points into the graph
			es[i<<1].len = len;
			if (len) {
				es[i<<1|1].seq = (const char*)malloc(len); // heap copy; freed by gfa_edseq_destroy
				if (es[i<<1|1].seq == 0) goto fail;
				memcpy((char*)es[i<<1|1].seq, ed + off, len);
				es[i<<1|1].len = len;
				off += len;
			} else {
				es[i<<1|1].seq = 0;
				es[i<<1|1].len = 0;
			}
		}
		assert(off == edseq_len);
	}
	*es_ = es;
	return g;
fail:
	gfa_destroy(g);
	return 0;
}

mg_idx_t *mg_tcache_load(FILE *fp)
{
	int fd;
	struct stat st;
	uint8_t *map;
	uint64_t size;
	tg_layout_t L;
	mg_idx_t *gi;
	gfa_t *g;
	gfa_edseq_t *es;
	uint32_t n_bucket, n_seg, n_sseq, max_rank, w, k, b, has_link_aux, i;
	uint64_t n_arc, n_entry, p_total, sum_len, names_len, seqs_len, sseq_names_len;
	uint64_t edseq_len, aux_len, bblob_size;
	double t0, t_mmap, t_crc;

	if (fp == 0) return 0;
	fd = fileno(fp);
	if (fd < 0 || fstat(fd, &st) != 0) return 0;
	size = (uint64_t)st.st_size;
	if (size < TG_HDR_SIZE + 4) return 0;
#ifdef WIN32
	if (_ftelli64(fp) >= (int64_t)size) return 0;
#else
	if (ftello(fp) >= (off_t)size) return 0; // stream already consumed
#endif
	t0 = realtime();
	map = (uint8_t*)mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (map == MAP_FAILED) return 0;
	t_mmap = realtime();

	/* magic + version */
	if (memcmp(map + TG_OFF_MAGIC, TG_MAGIC, 4) != 0) {
		fprintf(stderr, "[ERROR] mg_tcache_load: not a minigraph tcache file (bad magic)\n");
		goto fail;
	}
	if (tg_rd16(map + TG_OFF_VERSION) != TG_VERSION) {
		fprintf(stderr, "[ERROR] mg_tcache_load: unsupported tcache format version %u (this build reads v%u); regenerate the file with: minigraph -i <ref>.tgcache <ref>.gfa\n",
			(unsigned)tg_rd16(map + TG_OFF_VERSION), (unsigned)TG_VERSION);
		goto fail;
	}

	/* whole-file CRC32C trailer: covers [0, size-4) — TracEon v4 integrity pattern */
	{
		uint32_t want = tg_rd32(map + size - 4);
		uint32_t got = mg_crc32c(map, size - 4);
		if (want != got) {
			fprintf(stderr, "[ERROR] mg_tcache_load: CRC32C mismatch (%s): stored %08x, computed %08x\n",
				"file corrupted or truncated", want, got);
			goto fail;
		}
	}
	t_crc = realtime();

	k = tg_rd32(map + TG_OFF_K); w = tg_rd32(map + TG_OFF_W); b = tg_rd32(map + TG_OFF_B);
	n_bucket = tg_rd32(map + TG_OFF_NBUCKET); n_seg = tg_rd32(map + TG_OFF_NSEG);
	n_sseq = tg_rd32(map + TG_OFF_NSSEQ); max_rank = tg_rd32(map + TG_OFF_MAXRANK);
	n_arc = tg_rd64(map + TG_OFF_NARC); n_entry = tg_rd64(map + TG_OFF_NENTRY);
	p_total = tg_rd64(map + TG_OFF_PTOTAL); sum_len = tg_rd64(map + TG_OFF_SUMLEN);
	names_len = tg_rd64(map + TG_OFF_NAMESLEN); seqs_len = tg_rd64(map + TG_OFF_SEQSLEN);
	sseq_names_len = tg_rd64(map + TG_OFF_SSEQNAMESLEN);
	edseq_len = tg_rd64(map + TG_OFF_EDSEQLEN); aux_len = tg_rd64(map + TG_OFF_AUXLEN);
	has_link_aux = (tg_rd16(map + TG_OFF_FLAGS) & TG_F_HAS_LINK_AUX) != 0;

	if (b < 1 || b > 30) goto fail;
	if (n_bucket != 1U << b) goto fail;
	if (n_seg != 0 && sum_len == 0) goto fail;
	/* The bucket-blob size is stored implicitly as b_off[n_bucket]; compute the
	 * partial layout to locate the b_off table, read its last entry, then
	 * recompute the full layout with the true blob size. */
	tg_calc_layout(sum_len, n_seg, n_sseq, n_bucket, n_arc, names_len, seqs_len,
	               sseq_names_len, edseq_len, aux_len, 0, p_total, &L);
	bblob_size = tg_rd64(map + L.boff_off + (uint64_t)n_bucket * 8);
	tg_calc_layout(sum_len, n_seg, n_sseq, n_bucket, n_arc, names_len, seqs_len,
	               sseq_names_len, edseq_len, aux_len, bblob_size, p_total, &L);
	if (L.payload_size + 4 != size) goto fail;

	g = mg_tcache_load_graph(map, &L, n_seg, n_sseq, max_rank, n_arc, names_len,
	                         seqs_len, sseq_names_len, aux_len, edseq_len, has_link_aux, &es);
	if (g == 0) goto fail;

	gi = mg_idx_init(k, w, b);
	if (gi == 0) { gfa_destroy(g); goto fail; }
	gi->is_tcache = 1;
	gi->tcache_map = map;
	gi->tcache_size = (int64_t)size;
	gi->g = g;
	gi->g_own = g; // the reconstructed graph is owned by gi (freed by mg_idx_destroy)
	gi->es = es;
	gi->n_seg = n_seg;

	/* table block: point every bucket at its mmapped slot array + occupancy
	 * bitmap — ZERO inserts; the entry count comes from popcounting the bitmap */
	{
		const uint64_t *poff = (const uint64_t*)(map + L.poff_off);
		const uint64_t *boff = (const uint64_t*)(map + L.boff_off);
		const uint8_t *bblob = map + L.bblob_off;
		const uint64_t *pblob = (const uint64_t*)(map + L.pblob_off);
		for (i = 0; i < n_bucket; ++i) {
			mg_idx_bucket_t *b2 = &gi->B[i];
			uint64_t pe = poff[i + 1] - poff[i];
			uint64_t bs = boff[i + 1] - boff[i];
			b2->n = (int32_t)pe;
			if (pe) b2->p = (uint64_t*)(pblob + poff[i]);
			if (bs) {
				const uint8_t *t = bblob + boff[i];
				uint32_t cap = tg_rd32(t), bits = tg_rd32(t + 4);
				uint64_t expect = tg_round8(8 + (uint64_t)cap * 16 + ((uint64_t)cap + 31) / 32 * 4);
				if (cap < 4 || (cap & (cap - 1)) != 0 || bits != 31 - (uint32_t)__builtin_clz(cap) || expect != bs) goto fail;
				b2->cap = cap;
				b2->fe = (const uint64_t*)(t + 8);
				b2->bm = t + 8 + (uint64_t)cap * 16; // khashl u32 bitmap, viewed as bytes
				b2->ne = (int32_t)tg_bitmap_popcount((const uint32_t*)b2->bm, cap);
			}
		}
	}

	/* advance the FILE position to EOF so a hypothetical multi-part loop ends */
#ifdef WIN32
	_fseeki64(fp, (int64_t)size, SEEK_SET);
#else
	fseeko(fp, (off_t)size, SEEK_SET);
#endif
	if (mg_verbose >= 3) {
		double t1 = realtime();
		fprintf(stderr, "[M::%s] mmap %.1f ms, crc32c %.1f ms, fixup %.1f ms (total %.1f ms) — %lld entries, %lld p-entries, %lld ref bases, %u segments\n",
			__func__, (t_mmap - t0) * 1e3, (t_crc - t_mmap) * 1e3, (t1 - t_crc) * 1e3,
			(t1 - t0) * 1e3, (long long)n_entry, (long long)p_total, (long long)sum_len, n_seg);
	}
	return gi;

fail:
	munmap(map, size);
	return 0;
}

int mg_tcache_is_tcache(const char *fn)
{
	int fd, ret = 0;
	char magic[4];
	if (fn == 0 || strcmp(fn, "-") == 0) return 0;
	fd = open(fn, O_RDONLY);
	if (fd >= 0) {
		if (read(fd, magic, 4) == 4 && memcmp(magic, TG_MAGIC, 4) == 0)
			ret = 1;
		close(fd);
	}
	return ret;
}

#endif /* TRACEON_BACKEND */
