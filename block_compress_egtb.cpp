#include <omp.h>
#include <iostream>
#include <cassert>
#include <filesystem>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <cstring>
#include <math.h>
#include <unistd.h> // close
#include <fcntl.h> // open
#include <sys/mman.h> // mmap
#include <sys/stat.h> // file size
#include "misc.h"
namespace fs = std::filesystem;

// #define LIBDEFLATE
// #define ZSTD
// #define LZ4

#define EXT ".bz"

#include <zstd.h>

uint64_t compute_checksum(int16_t* TB, uint64_t num_pos, int nthreads) {
  uint64_t s = 0;
    #pragma omp parallel for num_threads(nthreads) schedule(static) reduction(^:s)
    for (uint64_t ix = 0; ix < num_pos; ix++) {
      s ^= (ix << 10) | ((uint64_t) abs(TB[ix]));
    }
    return s;
}

uint64_t compress_egtb(std::string filename, int nthreads, int compression_level, uint64_t block_size, bool write, bool verbose) {

  // load file
  FILE* f = fopen(filename.c_str(), "rb");
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


  if (verbose) std::cout << "Read " << num_pos << " entries (" << (double) num_pos*2 / (1024*1024*1024) << " GiB) from " << filename << std::endl;

  uint64_t checksum = compute_checksum(TB, num_pos, nthreads);

  uint64_t n_blocks = ceil((double) num_pos / block_size);
  uint64_t* block_sizes = (uint64_t*) malloc(n_blocks * sizeof(uint64_t));  

  bool large_blocks = block_size > uint64_t(UINT16_MAX); 
  uint64_t final_size = 0;

  #pragma omp parallel num_threads(nthreads)
  {
    // thread-private allocations
    std::vector<uint8_t> compressed(block_size*2);

    
    #pragma omp for schedule(static) reduction(+:final_size)
    for (uint64_t start_ix = 0; start_ix < num_pos; start_ix += block_size) {
      uint64_t size = std::min(block_size, num_pos - start_ix) * sizeof(int16_t);
      uint64_t block_ix = start_ix / block_size;

      uint64_t compressed_size = ZSTD_compress(compressed.data(), compressed.size(),  (uint8_t*) &TB[start_ix], size, compression_level);
      assert (!ZSTD_isError(compressed_size));


      if (compressed_size == 0) {
        std::cerr << "Compression failed\n";
        exit(1);
      }
      final_size += compressed_size;

      block_sizes[block_ix] = compressed_size;
      std::memcpy(&TB[start_ix], compressed.data(), compressed_size);
    }

  }

  final_size += large_blocks ? sizeof(uint32_t) * n_blocks : sizeof(uint16_t) * n_blocks;
  final_size += sizeof(uint64_t); // checksum
  final_size += sizeof(uint64_t); // num_pos
  final_size += sizeof(uint64_t); // block_size

  if (verbose) std::cout << "final size: " << final_size << " " << (double) final_size / file_size << " (large_blocks = " << large_blocks << ")" << std::endl;


  if (write) {
    std::string compressed_filename = filename + EXT;
    f = fopen(compressed_filename.c_str(), "wb");
    if (f == NULL) {
      std::cerr << "Error opening file: " << compressed_filename << std::endl;
      exit(1);
    }
    fwrite(&checksum, sizeof(uint64_t), 1, f);
    fwrite(&num_pos, sizeof(uint64_t), 1, f);
    fwrite(&block_size, sizeof(uint64_t), 1, f);
    if (large_blocks) {
      uint32_t* write_block_sizes = (uint32_t*) malloc(n_blocks * sizeof(uint32_t));  
      for (uint64_t i = 0; i < n_blocks; i++) write_block_sizes[i] = block_sizes[i];
      fwrite(write_block_sizes, sizeof(uint32_t), n_blocks, f);
    } else {
      uint16_t* write_block_sizes = (uint16_t*) malloc(n_blocks * sizeof(uint16_t));  
      for (uint64_t i = 0; i < n_blocks; i++) write_block_sizes[i] = block_sizes[i];
      fwrite(write_block_sizes, sizeof(uint16_t), n_blocks, f);
    }

    for (uint64_t start_ix = 0; start_ix < num_pos; start_ix += block_size) {
      uint64_t block_ix = start_ix / block_size;
      fwrite(&TB[start_ix], sizeof(uint8_t), block_sizes[block_ix], f);
      assert(block_ix * block_size == start_ix);
    }
    if (verbose) std::cout << "Wrote to " << compressed_filename << std::endl;
    fclose(f);
  }

  free(block_sizes);
  free(TB);
  return final_size;
}


int main_test(int argc, char *argv[]) {
  assert(argc == 4);
  uint64_t n_threads = atoi(argv[1]);
  uint64_t compression_level = atoi(argv[2]);
  uint64_t block_size = atoi(argv[3]); // [4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152]

  std::cout << "n_threads: " << n_threads << ", compression_level: " << compression_level << ", block_size: " << block_size << std::endl;
  bool verbose = true;
  bool write = true;
  bool test = true;

  uint64_t count = 0;
  uint64_t total_size = 0;
  uint64_t total_compressed_size = 0;
  uint64_t total_zipped_size = 0;
  std::string path = "egtbs/5men/";
  std::vector<std::string> files;
  TimePoint t0 = now();
  for (const auto & entry : fs::recursive_directory_iterator(path)) {
    if (entry.path().extension() != ".egtb") continue;
    count++;
    files.push_back(entry.path().u8string());
    uint64_t size = entry.file_size();
    total_size += size;
    uint64_t zipped_size = fs::file_size(fs::path(entry.path().u8string() + ".zip"));
    total_zipped_size += zipped_size;
    uint64_t compressed_size = compress_egtb(entry.path().u8string(), n_threads, compression_level, block_size, write, verbose);
    total_compressed_size += compressed_size;
    if (verbose) {
      std::cout << entry.path() << " size: " << size;
      std::cout << " zipped size: " << zipped_size << " (" << (double) zipped_size / size << ")";
      std::cout << " compressed size: " << compressed_size << " (" << (double) compressed_size / size << ")";
      std::cout << (compressed_size <= zipped_size ? " OK" : " :(") << std::endl;
      std::cout << std::endl;
    }
  }
  std::cout << "Compressed " << count << " files" << std::endl;
  std::cout << "TOTAL" << " size: " << total_size;
  std::cout << " zipped size: " << total_zipped_size << " (" << (double) total_zipped_size / total_size << ")";
  std::cout << " compressed size: " << total_compressed_size << " (" << (double) total_compressed_size / total_size << ")" << std::endl;

  TimePoint t1 = now();
  std::cout << "Finished in " << (double) (t1 - t0) / 1000 << " s" << std::endl;

  if (test) {
    CompressedEGTB** cegtbs = (CompressedEGTB**) malloc(count * sizeof(CompressedEGTB*));
    for (uint64_t i = 0; i < count; i++) {
      cegtbs[i] = new CompressedEGTB(files[i] + EXT);
      int16_t* decompressed_TB = cegtbs[i]->decompress_into_array(n_threads);
      free(decompressed_TB);
    }
    
    int n = 100000;
    t0 = now();
    for (int i = 0; i < n; i++) {
      int fix = rand() % count;
      CompressedEGTB* cegtb = cegtbs[fix];
      uint64_t ix = rand() % cegtb->num_pos;
      cegtb->get(ix);
    }
    t1 = now();
    std::cout << "Mean access time " << (double) (t1 - t0) /n << "ms vs filesize " << (double) total_compressed_size / (1024*1024*1024) <<  "GiB" << std::endl;
    
    for (uint64_t i = 0; i < count; i++) {
      free(cegtbs[i]);
    }
    free(cegtbs);
  }

  return 0;
}

void run_cmd(std::string cmd) {
  std::cout << "CMD " << cmd << std::endl;
  system(cmd.c_str());
}

int main_compress(int argc, char *argv[]) {
  assert(argc == 4);
  uint64_t n_threads = atoi(argv[1]);
  uint64_t compression_level = atoi(argv[2]);
  uint64_t block_size = atoi(argv[3]);

  std::cout << "n_threads: " << n_threads << ", compression_level: " << compression_level << ", block_size: " << block_size << std::endl;
  bool verbose = true;
  bool write = true;

  uint64_t count = 0;
  uint64_t total_size = 0;
  uint64_t total_compressed_size = 0;
  uint64_t total_zipped_size = 0;
  std::string path = "egtbs/";
  for (const auto & entry : fs::recursive_directory_iterator(path)) {
    if (entry.path().extension() != ".zip") continue;
    TimePoint t0 = now();
    count++;
    std::string zipped_filename = entry.path().u8string();
    std::string filename = zipped_filename.substr(0, zipped_filename.size() - 4);
    std::cout << count << ": " << zipped_filename << " " << filename << std::endl;

    run_cmd("unzip -n " + zipped_filename);
    run_cmd("md5sum " + filename);

    uint64_t zipped_size = entry.file_size();
    total_zipped_size += zipped_size;
    uint64_t size = fs::file_size(fs::path(filename));
    total_size += size;
    uint64_t compressed_size = compress_egtb(filename, n_threads, compression_level, block_size, write, verbose);
    total_compressed_size += compressed_size;

    run_cmd("rm " + zipped_filename);
    run_cmd("rm " + filename);

    TimePoint t1 = now();
    if (verbose) {
      std::cout << filename << " size: " << size;
      std::cout << " zipped size: " << zipped_size << " (" << (double) zipped_size / size << ")";
      std::cout << " compressed size: " << compressed_size << " (" << (double) compressed_size / size << ")";
      std::cout << (compressed_size <= zipped_size ? " OK" : " :(") << std::endl;
      std::cout << "Finished in " << (double) (t1 - t0) / 1000 << " s" << std::endl;
      std::cout << std::endl;
    }

  }
  std::cout << "Compressed " << count << " files" << std::endl;
  std::cout << "TOTAL" << " size: " << total_size;
  std::cout << " zipped size: " << total_zipped_size << " (" << (double) total_zipped_size / total_size << ")";
  std::cout << " compressed size: " << total_compressed_size << " (" << (double) total_compressed_size / total_size << ")" << std::endl;

  return 0;
}

int main(int argc, char *argv[]) {
  main_test(argc, argv);
  // main_compress(argc, argv); // ./compress/compress_zstd.out 32 19 1048576 |& tee -a log_compress.txt
}