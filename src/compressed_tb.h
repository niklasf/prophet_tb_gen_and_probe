#ifndef COMPRESSED_TB_H_INCLUDED
#define COMPRESSED_TB_H_INCLUDED

#ifdef OMP
#include <omp.h>
#endif
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
/*
.bz FORMAT:

HEADER 3 x UINT64 = 24 bytes
checksum UINT64
num_pos UINT64
block_size UINT64

BLOCK_SIZES
if block_size > 65535:
  n_blocks * UINT32
else
  n_blocks * UINT16
where  n_blocks = ceil((double) num_pos / block_size);

COMPRESSED DATA
individual compressed blocks concatenated at byte level
*/

namespace Prophet {

uint64_t compute_checksum(int16_t* TB, uint64_t num_pos, [[maybe_unused]] int nthreads);
uint64_t block_compress_TB(int16_t* TB, uint64_t num_pos, int nthreads, int compression_level, uint64_t block_size, std::string compressed_filename, bool write, bool verbose);

struct DecompressCtx {
  uint64_t buf_size;
  int16_t* uncompressed_buf;
  int egtb_ix;
  uint64_t uncompressed_block_ix;
  ZSTD_DCtx* decompressor;
  ZSTD_DStream* stream;
  uint64_t probe_count;
  DecompressCtx(uint64_t buf_size_ = 2097152) {
      buf_size = buf_size_;
    uncompressed_buf = (int16_t*) malloc(buf_size * sizeof(int16_t));
#if 0
    uncompressed_block_ix = UINT64_MAX;
    egtb_ix = 1001;
    decompressor = ZSTD_createDCtx();
#else
    stream = ZSTD_createDStream();
#endif
    probe_count = 0;
  }
  ~DecompressCtx() {
    free(uncompressed_buf);
#if 0
    ZSTD_freeDCtx(decompressor);
#else
    ZSTD_freeDStream(stream);
#endif
  }
};

struct CompressedTB {
  int egtb_ix;
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
  DecompressCtx* decompress_ctx;

  CompressedTB(int egtb_ix_, std::string filename) {
    // std::cout << "CompressedEGTB " << this << ": " << filename << std::endl;
    this->egtb_ix = egtb_ix_;
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
    uint64_t offset = 0;
    for (uint64_t block_ix = 0; block_ix < this->n_blocks; block_ix++) {
      this->block_offsets[block_ix] = offset;
      offset += this->block_sizes[block_ix];
    }

    decompress_ctx = nullptr;
    // std::cout << "CompressedEGTB" << this->egtb_id << ": filesize=" << this->filesize << ", num_pos=" << this->num_pos << ", block_size=" << this->block_size << ", n_blocks=" << this->n_blocks << std::endl;
  }

  ~CompressedTB() {
    // std::cout << "~CompressedEGTB " << this << ": " << this->compressed_filename << std::endl;
    int unmap = munmap(this->map_ptr, this->compressed_filesize);
    if (unmap == -1) {
        perror(("Error munmapping the file " + this->compressed_filename).c_str());
        exit(EXIT_FAILURE);
    }
    free(this->block_sizes);
    free(this->block_offsets);
    if (decompress_ctx != nullptr) delete decompress_ctx;
  }

  // not thread-safe
  inline int16_t get_value(uint64_t ix) {

    if (decompress_ctx == nullptr)
      decompress_ctx = new DecompressCtx(block_size);

    return get_value_dctx(ix, decompress_ctx);
  }

  // thread-safe
  int16_t get_value_dctx(uint64_t ix, DecompressCtx* dctx) {
    assert (dctx->buf_size >= block_size);
    dctx->probe_count++;

    uint64_t block_ix = ix / block_size;
    uint64_t ix_in_block = ix % block_size;
    uint64_t offset = block_offsets[block_ix];
    uint64_t size = std::min(block_size, num_pos - block_ix*block_size) * sizeof(int16_t);

#if 0
    if (dctx->egtb_ix != egtb_ix || dctx->uncompressed_block_ix != block_ix) {
      dctx->egtb_ix = egtb_ix; // DecompressCtx may be shared across CompressTB
      dctx->uncompressed_block_ix = block_ix;
      size_t res = ZSTD_decompressDCtx(dctx->decompressor,
        dctx->uncompressed_buf, size,
        &compressed_blocks[offset], block_sizes[block_ix]);
      assert (!ZSTD_isError(res));
    }
#else
    ZSTD_inBuffer in_buffer;
    in_buffer.src = &compressed_blocks[offset];
    in_buffer.size = size;
    in_buffer.pos = 0;

    ZSTD_outBuffer out_buffer;
    out_buffer.dst = dctx->uncompressed_buf;
    out_buffer.size = std::min(ix_in_block + 1, dctx->buf_size) * sizeof(int16_t);
    out_buffer.pos = 0;

    ZSTD_initDStream(dctx->stream);

    while (in_buffer.pos < in_buffer.size && out_buffer.pos < out_buffer.size) {
        size_t res = ZSTD_decompressStream(dctx->stream, &out_buffer, &in_buffer);
        assert(!ZSTD_isError(res));
    }
#endif
    return dctx->uncompressed_buf[ix_in_block];
  }


  void decompress_to_array(int nthreads, int16_t* TB) {
    #pragma omp parallel num_threads(nthreads)
    {
      DecompressCtx* dctx = new DecompressCtx(block_size);

      #pragma omp for schedule(static)
      for (uint64_t start_ix = 0; start_ix < num_pos; start_ix += block_size) {
        uint64_t block_ix = start_ix / block_size;
        uint64_t offset = block_offsets[block_ix];
        uint64_t size = std::min(block_size, num_pos - start_ix) * sizeof(int16_t);
        size_t res = ZSTD_decompressDCtx(dctx->decompressor,
          &TB[start_ix], size,
          &compressed_blocks[offset], block_sizes[block_ix]);
        assert (!ZSTD_isError(res));
      }

      #pragma omp critical
      {
        delete dctx;
      }
    }

    assert (compute_checksum(TB, num_pos, nthreads) == checksum);
  }

  void decompress_to_file(std::string filename) {
    FILE* f = fopen(filename.c_str(), "wb");
    if (f == NULL) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }
    if (decompress_ctx == nullptr)
      decompress_ctx = new DecompressCtx(block_size);
    assert (decompress_ctx->buf_size >= block_size);
    for (uint64_t start_ix = 0; start_ix < num_pos; start_ix += block_size) {
        uint64_t block_ix = start_ix / block_size;
        uint64_t offset = block_offsets[block_ix];
        uint64_t size = std::min(block_size, num_pos - start_ix) * sizeof(int16_t);
        size_t res = ZSTD_decompressDCtx(decompress_ctx->decompressor,
            decompress_ctx->uncompressed_buf, size,
            &compressed_blocks[offset], block_sizes[block_ix]);
        assert (!ZSTD_isError(res));
        fwrite(decompress_ctx->uncompressed_buf, sizeof(uint8_t), size, f);
    }
    fclose(f);
  }

  bool check_integrity([[maybe_unused]] int nthreads);
};

} // namespace Prophet

#endif
