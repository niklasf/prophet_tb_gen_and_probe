
#include <iostream>
#include <cassert>
#include "egtb.h"
#include "compressed_tb.h"
#include "misc.h"
#include <filesystem>
namespace fs = std::filesystem;

void load_file(std::string filename, int16_t*& TB, u_int64_t& num_pos, bool verbose) {
  // load file
  FILE* f = fopen(filename.c_str(), "rb");
  if (f == NULL) {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }
  fseek(f, 0, SEEK_END);
  uint64_t file_size = ftell(f);
  rewind(f);
  
  num_pos = file_size / sizeof(int16_t);
  TB = (int16_t*) malloc(num_pos * sizeof(int16_t));
  
  fread(TB, sizeof(int16_t), num_pos, f);
  fclose(f);
  if (verbose) std::cout << "Read " << num_pos << " entries (" << (double) num_pos*2 / (1024*1024*1024) << " GiB) from " << filename << std::endl;
}

// legacy code that was used to compress from .zip

int test_from_zip(int argc, char *argv[]) {
  assert(argc == 4);
  uint64_t n_threads = atoi(argv[1]);
  uint64_t compression_level = atoi(argv[2]);
  uint64_t block_size = atoi(argv[3]); 

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
    std::string filename = entry.path().u8string();
    files.push_back(filename);
    uint64_t size = entry.file_size();
    total_size += size;
    uint64_t zipped_size = fs::file_size(fs::path(filename + ".zip"));
    total_zipped_size += zipped_size;

    int16_t* TB;
    uint64_t num_pos;
    load_file(filename, TB, num_pos, verbose);
    uint64_t compressed_size = compress_egtb(TB, num_pos, n_threads, compression_level, block_size, filename + COMP_EXT, write, verbose);
    free(TB);

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
    CompressedTB** ctbs = (CompressedTB**) malloc(count * sizeof(CompressedTB*));
    for (uint64_t i = 0; i < count; i++) {
      ctbs[i] = new CompressedTB(i, files[i] + COMP_EXT);
      int16_t* decompressed_TB = (int16_t*) malloc(ctbs[i]->num_pos * sizeof(int16_t));
      ctbs[i]->decompress_to_array(n_threads, decompressed_TB); // performs checksum test
      free(decompressed_TB);
    }
    
    int n = 100000;
    t0 = now();
    for (int i = 0; i < n; i++) {
      int fix = rand() % count;
      CompressedTB* ctb = ctbs[fix];
      uint64_t ix = rand() % ctb->num_pos;
      ctb->get_value(ix);
    }
    t1 = now();
    std::cout << "Mean access time " << (double) (t1 - t0) /n << "ms vs filesize " << (double) total_compressed_size / (1024*1024*1024) <<  "GiB" << std::endl;
    
    for (uint64_t i = 0; i < count; i++) {
      free(ctbs[i]);
    }
    free(ctbs);
  }

  return 0;
}

void run_cmd(std::string cmd) {
  std::cout << "CMD " << cmd << std::endl;
  system(cmd.c_str());
}

// legacy code that was used to compress from .zip

int compress_from_zip(int argc, char *argv[]) {
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

    int16_t* TB;
    uint64_t num_pos;
    load_file(entry.path().u8string(), TB, num_pos, verbose);
    uint64_t compressed_size = compress_egtb(TB, num_pos, n_threads, compression_level, block_size, filename + COMP_EXT, write, verbose);
    free(TB);

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
  // [4096, 8192, 16384, 32768, , 131072, 262144, 524288, 1048576, 2097152]
  // test_from_zip(argc, argv);
  // compress_from_zip(argc, argv);

  assert(argc == 4);
  uint64_t n_threads = atoi(argv[1]);
  uint64_t compression_level = atoi(argv[2]);
  uint64_t block_size = atoi(argv[3]);

  bool verbose = true;
  bool write = true;
  bool test = true;

  std::cout << "n_threads: " << n_threads << ", compression_level: " << compression_level << ", block_size: " << block_size << std::endl;

  uint64_t count = 0;
  uint64_t total_size = 0;
  uint64_t total_compressed_size = 0;
  
  std::string path = "egtbs";
  std::vector<std::string> files;
  
  TimePoint t0 = now();
  for (const auto & entry : fs::recursive_directory_iterator(path)) {
    if (entry.path().extension() != COMP_EXT) continue;
    std::string bz_filename = entry.path().u8string();
    files.push_back(bz_filename);
    count++;
    std::string filename = bz_filename.substr(0, bz_filename.size() - 3);
    CompressedTB ctb = CompressedTB(0, bz_filename);
    uint64_t size = ctb.num_pos * sizeof(int16_t);
    total_size += size;
    uint64_t compressed_size;
    TimePoint tt0 = now();
    if (ctb.block_size != block_size) {
      int16_t* TB = (int16_t*) malloc(size);
      ctb.decompress_to_array(n_threads, TB);
      compressed_size = compress_egtb(TB, ctb.num_pos, n_threads, compression_level, block_size, filename + COMP_EXT, write, verbose);
      free(TB);
    } else {
      compressed_size = ctb.compressed_filesize;
    }
    TimePoint tt1 = now();
    total_compressed_size += compressed_size;
    if (verbose) {
      std::cout << filename << " size: " << size;
      std::cout << " compressed size: " << compressed_size << " (" << (double) compressed_size / size << ")";
      std::cout << std::endl;
      std::cout << "Finished in " << (double) (tt1 - tt0) / 1000 << " s" << std::endl;
    }
  }
  TimePoint t1 = now();

  std::cout << "Compressed " << count << " files" << std::endl;
  std::cout << "TOTAL" << " size: " << total_size;
  std::cout << " compressed size: " << total_compressed_size << " (" << (double) total_compressed_size / total_size << ")" << std::endl;
  std::cout << "Finished in " << (double) (t1 - t0) / 1000 << " s" << std::endl;


  if (test) {
    CompressedTB** ctbs = (CompressedTB**) malloc(count * sizeof(CompressedTB*));
    for (uint64_t i = 0; i < count; i++) {
      ctbs[i] = new CompressedTB(i, files[i]);
      // int16_t* decompressed_TB = (int16_t*) malloc(ctbs[i]->num_pos * sizeof(int16_t));
      // ctbs[i]->decompress_to_array(n_threads, decompressed_TB); // performs checksum test
      // free(decompressed_TB);
    }
    
    int n = 100000;
    t0 = now();
    for (int i = 0; i < n; i++) {
      int fix = rand() % count;
      CompressedTB* ctb = ctbs[fix];
      uint64_t ix = rand() % ctb->num_pos;
      ctb->get_value(ix);
    }
    t1 = now();
    std::cout << "Mean access time " << (double) (t1 - t0) / n << "ms vs filesize " << (double) total_compressed_size / (1024*1024*1024) <<  "GiB" << std::endl;
    
    for (uint64_t i = 0; i < count; i++) {
      free(ctbs[i]);
    }
    free(ctbs);
  }
}