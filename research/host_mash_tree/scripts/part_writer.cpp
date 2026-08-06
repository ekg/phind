// part_writer.cpp
//
// Consumes `mash dist` stdout lines and writes a per-job part file of float32
// distances in TRIANGLE ORDER (row-major, a<b, row index = global a).
//
// Orientation contract (verified empirically):
//   mash dist REF QUERY emits lines in (query-outer, ref-inner) order.
//
// Usage:
//   part_writer <ids.txt> <mode> <start> <size>
//     mode "offdiag": caller must run `mash dist chunk_j.msh chunk_i.msh`
//                     with i<j (ref = b-side = later chunk). Then query ids
//                     are the a-side and ref ids the b-side; every emitted
//                     pair has a<b and the value stream is already
//                     (a-major, b-minor) = triangle row order. Writer appends
//                     the raw float32 value per line to stdout (binary).
//     mode "diag":    ref and query come from the same chunk; local indices
//                     are in [0,size). The writer buffers a size*size float32
//                     grid (ref=row x, query=col y), then dumps the upper
//                     triangle (x-major, y-minor, x<y) to stdout (binary).
//   <start> is the global index of the chunk's first sequence (diag only;
//   unused for offdiag).
//
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <unistd.h>

static std::unordered_map<std::string_view, uint32_t> build_index(
    const std::vector<std::string>& ids) {
  std::unordered_map<std::string_view, uint32_t> m;
  m.reserve(ids.size() * 2);
  for (uint32_t i = 0; i < ids.size(); ++i) m.emplace(std::string_view(ids[i]), i);
  return m;
}

static inline void write_val(float v) {
  // write to fd 1 (binary stdout)
  const char* p = reinterpret_cast<const char*>(&v);
  ssize_t w = write(1, p, 4);
  (void)w;
}

int main(int argc, char** argv) {
  if (argc != 5) {
    fprintf(stderr, "usage: %s <ids.txt> <offdiag|diag> <start> <size>\n", argv[0]);
    return 2;
  }
  std::vector<std::string> ids;
  {
    FILE* f = fopen(argv[1], "r");
    if (!f) { perror("ids"); return 2; }
    char* line = nullptr; size_t cap = 0; ssize_t len;
    while ((len = getline(&line, &cap, f)) > 0) {
      if (line[len - 1] == '\n') line[len - 1] = '\0';
      if (len > 1) ids.emplace_back(line);
    }
    free(line);
    fclose(f);
  }
  const bool diag = strcmp(argv[2], "diag") == 0;
  const uint64_t start = strtoull(argv[3], nullptr, 10);
  const uint64_t size = strtoull(argv[4], nullptr, 10);
  auto idx = build_index(ids);

  constexpr size_t BUF = 1 << 20;
  std::vector<char> buf(BUF);
  size_t used = 0;
  uint64_t lines = 0, kept = 0, skipped = 0;

  std::vector<float> grid;
  if (diag) grid.assign((size_t)size * size, 0.f);

  auto process_line = [&](char* s, size_t len) {
    ++lines;
    char* t1 = (char*)memchr(s, '\t', len);
    if (!t1) { ++skipped; return; }
    char* t2 = (char*)memchr(t1 + 1, '\t', len - (t1 + 1 - s));
    if (!t2) { ++skipped; return; }
    char* t3 = (char*)memchr(t2 + 1, '\t', len - (t2 + 1 - s));
    if (!t3) { ++skipped; return; }
    auto it1 = idx.find(std::string_view(s, t1 - s));        // ref id
    auto it2 = idx.find(std::string_view(t1 + 1, t2 - t1 - 1)); // query id
    if (it1 == idx.end() || it2 == idx.end()) { ++skipped; return; }
    uint32_t r = it1->second, q = it2->second;
    float d = strtof(t2 + 1, nullptr);
    if (d < 0.f) d = 0.f;
    if (d > 1.f) d = 1.f;
    if (diag) {
      // r = ref = x (row), q = query = y (col); local indices
      uint64_t x = (uint64_t)r - start, y = (uint64_t)q - start;
      if (x >= size || y >= size) { ++skipped; return; }
      grid[(size_t)x * size + y] = d;
      ++kept;
    } else {
      // ref = b-side (later chunk), query = a-side (earlier chunk)
      uint32_t a = q, b = r;
      if (a >= b) { ++skipped; return; }  // safety; should not happen
      write_val(d);
      ++kept;
    }
  };

  while (true) {
    ssize_t got = read(0, buf.data() + used, BUF - used);
    if (got <= 0) break;
    used += (size_t)got;
    char* startp = buf.data();
    while (true) {
      char* nl = (char*)memchr(startp, '\n', buf.data() + used - startp);
      if (!nl) break;
      process_line(startp, nl - startp);
      startp = nl + 1;
    }
    size_t remain = buf.data() + used - startp;
    if (remain && startp != buf.data()) memmove(buf.data(), startp, remain);
    used = remain;
  }
  if (used > 0) process_line(buf.data(), used);

  if (diag) {
    // dump upper triangle row-major: for x: for y in x+1..size-1
    uint64_t dumped = 0;
    for (uint64_t x = 0; x < size; ++x) {
      for (uint64_t y = x + 1; y < size; ++y) {
        write_val(grid[(size_t)x * size + y]);
        ++dumped;
      }
    }
    fprintf(stderr, "lines=%llu kept=%llu skipped=%llu dumped=%llu\n",
            (unsigned long long)lines, (unsigned long long)kept,
            (unsigned long long)skipped, (unsigned long long)dumped);
  } else {
    fprintf(stderr, "lines=%llu kept=%llu skipped=%llu\n",
            (unsigned long long)lines, (unsigned long long)kept,
            (unsigned long long)skipped);
  }
  return 0;
}
