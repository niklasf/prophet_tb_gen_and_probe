#ifndef COMPRESSED_TB_H_INCLUDED
#define COMPRESSED_TB_H_INCLUDED

#include <omp.h>
#include <iostream>
#include <cassert>
#include <zstd.h>
#include <math.h>
#include <unistd.h> // close
#include <fcntl.h> // open
#include <sys/mman.h> // mmap
#include <sys/stat.h> // file size
#include <cstdint>
#include <stdio.h>

#define COMP_EXT ".bz"

struct CompressedTB {
  std::string compressed_filename;
  uint64_t compressed_filesize;
  uint64_t checksum;
  uint64_t num_pos;
  uint64_t block_size;
  uint64_t n_blocks;
  uint64_t* block_sizes;
  uint64_t* block_offsets;
  uint8_t* map_ptr;
  uint8_t* compressed_blocks;
  uint64_t uncombressed_block_ix;
  int16_t* uncompressed_buf;
  ZSTD_DCtx* decompressor;

  CompressedTB(std::string filename) {
    this->compressed_filename = filename;

    struct stat st;
    stat(this->compressed_filename.c_str(), &st);
    this->compressed_filesize = st.st_size;
    int fd = open(this->compressed_filename.c_str(), O_RDONLY);
    if (fd == -1) {
        printf("Could not open file %s.\n", this->compressed_filename.c_str());
        exit(EXIT_FAILURE);
    }
    this->map_ptr = (uint8_t*) mmap(0, this->compressed_filesize, PROT_READ, MAP_SHARED, fd, 0);
    if (this->map_ptr == MAP_FAILED) {
        close(fd);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }
    close(fd);

    int headersize = 3*sizeof(uint64_t);
    this->checksum = ((uint64_t*) this->map_ptr)[0];
    this->num_pos = ((uint64_t*) this->map_ptr)[1];
    this->block_size = ((uint64_t*) this->map_ptr)[2];
    this->n_blocks = ceil((double) num_pos / block_size);
    
    bool large_blocks = block_size > uint64_t(UINT16_MAX);
    this->block_sizes = (uint64_t*) malloc(n_blocks * sizeof(uint64_t));
    if (large_blocks) {
      uint32_t* block_sizes_ptr = (uint32_t*) &this->map_ptr[headersize];
      for (uint64_t i = 0; i < n_blocks; i++) this->block_sizes[i] = block_sizes_ptr[i];
      this->compressed_blocks = &this->map_ptr[headersize + this->n_blocks*sizeof(uint32_t)];
    } else {
      uint16_t* block_sizes_ptr = (uint16_t*) &this->map_ptr[headersize];
      for (uint64_t i = 0; i < n_blocks; i++) this->block_sizes[i] = block_sizes_ptr[i];
      this->compressed_blocks = &this->map_ptr[headersize + this->n_blocks*sizeof(uint16_t)];
    }
    

    this->block_offsets = (uint64_t*) malloc(this->n_blocks * sizeof(uint64_t));
    this->uncombressed_block_ix = UINT64_MAX;
    this->uncompressed_buf = (int16_t*) malloc(block_size * sizeof(int16_t));

    uint64_t offset = 0;
    for (uint64_t block_ix = 0; block_ix < this->n_blocks; block_ix++) {
      this->block_offsets[block_ix] = offset;
      offset += this->block_sizes[block_ix];
    }

    this->decompressor = ZSTD_createDCtx();

    // std::cout << "CompressedEGTB" << this->egtb_id << ": filesize=" << this->filesize << ", num_pos=" << this->num_pos << ", block_size=" << this->block_size << ", n_blocks=" << this->n_blocks << std::endl;
  }

  ~CompressedTB() {
    // std::cout << "~CompressedEGTB" << this->egtb_id << std::endl;
    int unmap = munmap(this->map_ptr, this->compressed_filesize);
    if (unmap == -1) {
        perror(("Error munmapping the file " + this->compressed_filename).c_str());
        exit(EXIT_FAILURE);
    }
    free(this->block_sizes);
    free(this->block_offsets);
    free(this->uncompressed_buf);


    ZSTD_freeDCtx(this->decompressor);
  }

  int16_t get(uint64_t ix) {
    uint64_t block_ix = ix / block_size;
    uint64_t ix_in_block = ix % block_size;
    uint64_t offset = block_offsets[block_ix];
    uint64_t size = std::min(block_size, num_pos - block_ix*block_size) * sizeof(int16_t);
    if (uncombressed_block_ix != block_ix) {
      uncombressed_block_ix = block_ix;
      size_t res = ZSTD_decompressDCtx(decompressor,
        uncompressed_buf, size,
        &compressed_blocks[offset], block_sizes[block_ix]);
      assert (!ZSTD_isError(res));
    }
    return uncompressed_buf[ix_in_block];
  }

  void decompress_to_array(int nthreads, int16_t* TB) {
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (uint64_t start_ix = 0; start_ix < num_pos; start_ix += block_size) {
      uint64_t block_ix = start_ix / block_size;
      uint64_t offset = block_offsets[block_ix];
      uint64_t size = std::min(block_size, num_pos - start_ix) * sizeof(int16_t);
      size_t res = ZSTD_decompress(
        &TB[start_ix], size,
        &compressed_blocks[offset], block_sizes[block_ix]);
      assert (!ZSTD_isError(res));
    }
    // assert (compute_checksum(TB, num_pos, nthreads) == checksum);
  }

  void decompress_to_file(std::string filename) {
    FILE* f = fopen(filename.c_str(), "wb");
    if (f == NULL) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }
    for (uint64_t start_ix = 0; start_ix < num_pos; start_ix += block_size) {
        uint64_t block_ix = start_ix / block_size;
        uint64_t offset = block_offsets[block_ix];
        uint64_t size = std::min(block_size, num_pos - start_ix) * sizeof(int16_t);
        size_t res = ZSTD_decompress(
            uncompressed_buf, size,
            &compressed_blocks[offset], block_sizes[block_ix]);
        assert (!ZSTD_isError(res));
        fwrite(uncompressed_buf, sizeof(uint8_t), size, f);
    }
    fclose(f);
  }
};

#endif