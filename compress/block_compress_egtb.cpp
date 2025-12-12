#include "repair.h"
#include <omp.h>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

#define BLOCKSIZE uint64_t(8192)

void compress_block(int16_t* TB_BLOCK, int16_t min_nonzero_val, uint64_t size_w, uint64_t freq_local[], RDS* rds) {

  uint16_t char_size = 1000 - min_nonzero_val + 1 + 1; // [minval,...,1000] + 0

  rds->txt_len = size_w;
  rds->num_pairs = 0;

  // init sequence
  SEQ* seq = rds->seq;
  for (uint64_t i = 0; i < size_w; i++) {
    if (TB_BLOCK[i] == 0) {
      seq[i].code = 0;
    } else {
      seq[i].code = abs(TB_BLOCK[i]) - min_nonzero_val + 1;
    }
    seq[i].next = DUMMY_POS;
    seq[i].prev = DUMMY_POS;
  }

  // init hashtable
  PAIR** h_first = rds->h_first;
  uint64_t h_num = rds->h_num;
  for (uint64_t i = 0; i < primes[h_num]; i++) {
    h_first[i] = NULL;
  }
    
  // init priority queue
  uint64_t p_max = rds->p_max;
  PAIR** p_que = rds->p_que;
  for (uint64_t i = 0; i < p_max; i++) {
    p_que[i] = NULL;
  }


}

void compress_egtb(const char* filename, int nthreads) {

  // load file
  FILE *f = fopen(filename, "rb");
  if (f == NULL) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  fseek(f, 0, SEEK_END);
  uint64_t file_size = ftell(f);
  rewind(f);
  
  uint64_t num_pos = file_size / sizeof(int16_t);
  int16_t *TB = (int16_t*) malloc(num_pos * sizeof(int16_t));
  
  fread(TB, sizeof(int16_t), num_pos, f);
  fclose(f);


  // determine initial alphabet size
  uint16_t min_nonzero_val = 1000;
  #pragma omp parallel for num_threads(nthreads) schedule(static) reduction(min:min_nonzero_val)
  for (uint64_t ix = 0; ix < num_pos; ix++) {
    // 0 -> draw
    // negative even numbers -> loss in ply
    // positive odd numbers -> win in ply
    if (TB[ix] != 0) min_nonzero_val = std::min(min_nonzero_val, (uint16_t) abs(TB[ix]));
  }
  uint16_t char_size = 1000 - min_nonzero_val + 1 + 1; // [minval,...,1000] + 0
  std::cout << "Read " << num_pos << " entries (" << (double) num_pos*2 / (1024*1024*1024) << " GiB) from " << filename << " -> min_nonzero_val=" << min_nonzero_val << " and char_size=" << char_size << std::endl;
  
  

  #pragma omp parallel num_threads(nthreads)
  {
    // thread-private allocations

    uint64_t freq_local[UINT16_MAX] = {0}; // frequency count for global huffman encoding
    RDS* rds = (RDS*) malloc(sizeof(RDS));
    rds->seq = (SEQ*) malloc(BLOCKSIZE * sizeof(SEQ));

    uint64_t h_num = INIT_HASH_NUM;
    rds->h_num = h_num;
    rds->h_first = (PAIR**) malloc(sizeof(PAIR*) * primes[h_num]);

    uint64_t p_max = (uint64_t) ceil(sqrt((double) BLOCKSIZE));
    rds->p_max = p_max;
    rds->p_que = (PAIR**) malloc(sizeof(PAIR*) * p_max);;

    #pragma omp for schedule(static)
    for (uint64_t start_ix = 0; start_ix < num_pos; start_ix += BLOCKSIZE) {
      uint64_t size_w = std::min(BLOCKSIZE, num_pos - start_ix);
      compress_block(&TB[start_ix], min_nonzero_val, size_w, freq_local, rds);
    }

    uint64_t freq[UINT16_MAX] = {0};
    #pragma omp critical
    {
      for (int i = 0; i < UINT16_MAX; i++) {
        freq[i] += freq_local[i];
      }
    }

    destructRDS(rds);
  }


  /*

  uint64_t* lengths = (uint64_t*) calloc(UINT16_MAX, sizeof(uint64_t));
  
  uint64_t final_size = 0;
  double avg_num_rules = 0.0;
  double avg_seq_len = 0.0;
  uint64_t num_replace_pairs = 0;
  uint64_t n_batches = ceil((double)count/batchsize);
  for (uint64_t start = 0; start < count; start += batchsize) {
    
    uint64_t size_w = std::min(batchsize, count - start);
    
    
    // init hash table
    uint64_t h_num = INIT_HASH_NUM;
    Pair** h_first = (PAIR**) malloc(sizeof(PAIR*) * primes[h_num]);
    for (uint64_t i = 0; i < primes[h_num]; i++) {
      h_first[i] = NULL;
    }
    
    // init priority queue
    uint64_t p_max = (uint64_t) ceil(sqrt((double) size_w));
    PAIR** p_que = (PAIR**) malloc(sizeof(PAIR*) * p_max);
    for (uint64_t i = 0; i < p_max; i++) {
      p_que[i] = NULL;
    }
    
    // populate RDS object
    RDS* rds = (RDS*) malloc(sizeof(RDS));
    rds->txt_len = size_w;
    rds->seq = seq;
    rds->num_pairs = 0;
    rds->h_num = h_num;
    rds->h_first = h_first;
    rds->p_max = p_max;
    rds->p_que = p_que;
    
    initRDS_by_counting_pairs(rds);
    
    DICT* dict = createDict(rds->txt_len, CHAR_SIZE);
    
    
    uint64_t num_rules = 0;
    uint64_t num_replaced = 0;
    PAIR* max_pair;
    CODE new_code = 0;
    while ((max_pair = getMaxPair(rds)) != NULL) {
      num_rules++;
      new_code = addNewPair(dict, max_pair);
      num_replaced += replacePairs(rds, max_pair, new_code);
      num_replace_pairs++;
    }
    // std::cout << "num_pairs: " << rds->num_pairs << std::endl;
    assert(rds->num_pairs == 0);
    // CODE max_code = new_code;
    // std::cout << "Total number of replacements: " << num_replaced << std::endl;
    // std::cout << "Total number of rules: " << num_rules << std::endl;
    getCompSeq(rds, dict);
    // std::cout << "Compressed sequence length: " << comp_seq.len << std::endl;
    
    avg_num_rules += (double) num_rules / n_batches;
    avg_seq_len += (double) dict->seq_len / n_batches;
    
    // PAIR* max_pair = getMaxPair(rds);
    // replacePairs(rds, max_pair, 1001);
    
    // std::cout << "Number of distinct pairs: " << rds->num_pairs << std::endl;
    // for (uint64_t i = 0; i < p_max; i++) {
    //     PAIR* p = rds->p_que[i];
    //     if (p != NULL) {
    //         std::cout << "Frequency " << i << ": ";
    //         while (p != NULL) {
    //             std::cout << "(" << p->left << "," << p->right << ", freq=" << p->freq << ") ";
    //             p = p->p_next;
    //         }
    //         std::cout << std::endl;
    //     }
    // }
    EDICT* edict = convertDict(dict, CHAR_SIZE);
    uint64_t cfg_bits = EncodeCFG(edict, CHAR_SIZE);
    
    final_size += (cfg_bits / 8 + 1);
    // final_size += 2*comp_seq.len; 
    
    // global huff
    for (uint64_t i = 0; i < dict->seq_len; i++) freq[dict->comp_seq[i]]++;
    
    DestructEDict(edict);
    destructRDS(rds);
    
    // huff
    // uint64_t* freq = (uint64_t*) calloc(max_code+1, sizeof(uint64_t));
    // uint64_t* lengths = (uint64_t*) calloc(max_code+1, sizeof(uint64_t));
    // for (uint64_t i = 0; i < comp_seq.len; i++) freq[comp_seq.seq[i]]++;
    // Node* tree = buildTree(freq, max_code + 1);
    // extractLengths(tree, 0, lengths);
    
    // uint64_t huffcoded_seq_bits = 0;
    // uint64_t nonzero_symbols = 0;
    // for (CODE code = 0; code <= max_code; code++) {
    //   huffcoded_seq_bits += freq[code] * lengths[code];
    //   nonzero_symbols += (freq[code] > 0);
    // }
    // final_size += (huffcoded_seq_bits / 8) + 1;
    // final_size += nonzero_symbols; // optimal 1 byte length canonical
    
    // free(freq);
    // free(lengths);
    // delete tree;
  }
  
  // global huff
  Node* tree = buildTree(freq, UINT16_MAX);
  extractLengths(tree, 0, lengths);
  
  uint64_t huffcoded_seq_bits = 0;
  uint64_t nonzero_symbols = 0;
  for (CODE code = 0; code < UINT16_MAX; code++) {
    huffcoded_seq_bits += freq[code] * lengths[code];
    nonzero_symbols += (freq[code] > 0);
  }
  final_size += (huffcoded_seq_bits / 8) + 1;
  final_size += nonzero_symbols * 4; // this is upper bound
  final_size += n_batches * 2; // jump table
  
  free(freq);
  free(lengths);
  delete tree;
  
  
  std::cout << "final size: " << final_size << " " << (double) final_size / (count*2.0) << std::endl;
  std::cout << "avg_num_rules: " << avg_num_rules << ", avg_seq_len: " << avg_seq_len << ", num_replace_pairs: " << num_replace_pairs << std::endl;
  free(TB);
  
  return final_size;

  */
}


int main(int argc, char *argv[]) {
  
  // uint64_t total_size = 0;
  // uint64_t total_compressed_size = 0;
  // uint64_t total_zipped_size = 0;
  // std::string path = "5men/egtbs/";
  // for (const auto & entry : fs::directory_iterator(path)) {
  //   if (entry.path().extension() != ".egtb") continue;
  //   uint64_t size = entry.file_size();
  //   total_size += size;
  //   uint64_t zipped_size = fs::file_size(fs::path("5men/zipped/" + entry.path().filename().u8string() + ".zip"));
  //   total_zipped_size += zipped_size;
  //   uint64_t compressed_size = RunRepair(batchsize, entry.path().u8string().c_str());
  //   total_compressed_size += compressed_size;
  //   std::cout << entry.path() << " size: " << size;
  //   std::cout << " zipped size: " << zipped_size << " (" << (double) zipped_size / size << ")";
  //   std::cout << " compressed size: " << compressed_size << " (" << (double) compressed_size / size << ")";
  //   std::cout << (compressed_size <= zipped_size ? " OK" : " :(") << std::endl;
  //   std::cout << std::endl;
  // }
  // std::cout << "TOTAL" << " size: " << total_size << " zipped size: " << total_zipped_size << " (" << (double) total_zipped_size / total_size << ") compressed size: " << total_compressed_size << " (" << (double) total_compressed_size / total_size << ")" << std::endl;
  // exit(1);
  
  compress_egtb("../egtbs/5men/1pawns/KRPKQ.egtb", 1);

  return 0;
}