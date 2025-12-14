#include "repair.h"
#include "huffman.h"
#include "encoder.h"
#include <omp.h>
#include <iostream>
#include <cassert>
#include <filesystem>
#include <vector>
namespace fs = std::filesystem;
#include "libdeflate.h"

uint64_t compress_egtb(const char* filename, int nthreads, int compression_level, uint64_t block_size) {

  // load file
  FILE *f = fopen(filename, "rb");
  if (f == NULL) {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }
  fseek(f, 0, SEEK_END);
  uint64_t file_size = ftell(f);
  rewind(f);
  
  uint64_t num_pos = file_size / sizeof(int16_t);
  int16_t *TB = (int16_t*) malloc(num_pos * sizeof(int16_t));
  
  fread(TB, sizeof(int16_t), num_pos, f);
  fclose(f);


  std::cout << "Read " << num_pos << " entries (" << (double) num_pos*2 / (1024*1024*1024) << " GiB) from " << filename << std::endl;
  
  uint64_t n_blocks = ceil((double) num_pos / block_size);

  bool large_blocks = false; 
  uint64_t final_size = 0;

  #pragma omp parallel num_threads(nthreads)
  {
    // thread-private allocations
    libdeflate_compressor* compressor = libdeflate_alloc_compressor(compression_level);
    std::vector<uint8_t> compressed(block_size*2);

    
    #pragma omp for schedule(static) reduction(+:final_size) reduction(||: large_blocks)
    for (uint64_t start_ix = 0; start_ix < num_pos; start_ix += block_size) {
      uint64_t size = std::min(block_size, num_pos - start_ix) * 2;

      uint64_t compressed_size = libdeflate_deflate_compress(compressor, (uint8_t*) &TB[start_ix], size, compressed.data(), compressed.size());

      if (compressed_size == 0) {
        std::cerr << "Compression failed\n";
        exit(1);
      }
      large_blocks = large_blocks || (compressed_size > UINT16_MAX);
      final_size += compressed_size;

    }

    libdeflate_free_compressor(compressor);
  }

  final_size += large_blocks ? 4 * n_blocks : 2 * n_blocks;

  std::cout << "final size: " << final_size << " " << (double) final_size / file_size << " (large_blocks = " << large_blocks << ")" << std::endl;

  free(TB);

  return final_size;
}


int main(int argc, char *argv[]) {
  assert(argc >= 3);
  uint64_t n_threads = atoi(argv[1]);
  uint64_t compression_level = atoi(argv[2]);
  uint64_t block_size = atoi(argv[3]); // [4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152]
  std::cout << "n_threads: " << n_threads << ", compression_level: " << compression_level << ", block_size: " << block_size << std::endl;
  
  uint64_t count = 0;
  uint64_t total_size = 0;
  uint64_t total_compressed_size = 0;
  uint64_t total_zipped_size = 0;
  std::string path = "../egtbs/5men/";
  for (const auto & entry : fs::recursive_directory_iterator(path)) {
    if (entry.path().extension() != ".egtb") continue;
    count++;
    uint64_t size = entry.file_size();
    total_size += size;
    uint64_t zipped_size = fs::file_size(fs::path(entry.path().u8string() + ".zip"));
    total_zipped_size += zipped_size;
    uint64_t compressed_size = compress_egtb(entry.path().u8string().c_str(), n_threads, compression_level, block_size);
    total_compressed_size += compressed_size;
    std::cout << entry.path() << " size: " << size;
    std::cout << " zipped size: " << zipped_size << " (" << (double) zipped_size / size << ")";
    std::cout << " compressed size: " << compressed_size << " (" << (double) compressed_size / size << ")";
    std::cout << (compressed_size <= zipped_size ? " OK" : " :(") << std::endl;
    std::cout << std::endl;
  }
  std::cout << "Compressed " << count << " files" << std::endl;
  std::cout << "TOTAL" << " size: " << total_size;
  std::cout << " zipped size: " << total_zipped_size << " (" << (double) total_zipped_size / total_size << ")";
  std::cout << " compressed size: " << total_compressed_size << " (" << (double) total_compressed_size / total_size << ")" << std::endl;
  exit(1);

  // zip size: 8709004054 (0.144411)

  // compression_level = 9 block_size = 32768 -> 8563460076 (0.141998) 3m20
  // compression_level = 10 block_size = 32768 -> 7936397012 (0.1316) 7m
  // compression_level = 12 block_size = 32768 -> 7785773584 (0.129102) 19m

  // compression_level = 9 block_size = 8192 -> 9137047977 (0.151509) 3m
  // compression_level = 10 block_size = 8192 -> 8581590257 (0.142298) 6m20
  // compression_level = 12 block_size = 8192 -> 8416370763 (0.139559) 17m

  // compression_level = 10 block_size = 1048576 ->  7627181263 (0.126473) 8m
  
  assert(argc == 5);
  const char* filename = argv[4];
  // ../egtbs/5men/1pawns/KRKRP.egtb

  compress_egtb(filename, n_threads, compression_level, block_size);

  return 0;
}