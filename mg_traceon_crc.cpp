/*
 * mg_traceon_crc.cpp — thin C shim over TracEon's header-only include/Crc32c.h
 *
 * TracEon's CRC-32C (Castagnoli) lives in a C++ header (namespace TracEon,
 * constexpr table, target("sse4.2") intrinsics). minigraph's tcache code is C
 * (mg_traceon_cache.c), so this TU is the one place that includes Crc32c.h and
 * exposes the same whole-payload checksum to C: raw accumulator init
 * 0xFFFFFFFF, final XOR 0xFFFFFFFF — identical to the `.traceon` v4 format and
 * to the minimap2/mm2-fast TRC2 tcache trailers.
 *
 * Compiled with $(CXX) (see the TRACEON block in Makefile), which resolves
 * libstdc++/libtraceon_kmer symbols at the final link. The whole file is
 * linked only into TRACEON builds.
 *
 * MIT License, copyright (c) 2026 — TracEon integration (see LICENSE.txt).
 */
#include "Crc32c.h"

#include <cstddef>
#include <cstdint>
#include <new>

extern "C" {

struct mg_crc32c_s { TracEon::Crc32c c; };

mg_crc32c_s *mg_crc32c_new(void) noexcept {
    try { return new mg_crc32c_s(); }
    catch (...) { return nullptr; }
}

void mg_crc32c_free(mg_crc32c_s *s) noexcept {
    delete s;
}

void mg_crc32c_update(mg_crc32c_s *s, const void *data, size_t len) noexcept {
    s->c.update(data, len);
}

uint32_t mg_crc32c_final(mg_crc32c_s *s) noexcept {
    return s->c.finalize();
}

uint32_t mg_crc32c(const void *data, size_t len) noexcept {
    return TracEon::crc32c(data, len);
}

} // extern "C"
