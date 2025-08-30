
#ifndef GEN_EGTB_H_INCLUDED
#define GEN_EGTB_H_INCLUDED

#include "kkx.h"
#include "linearize.h"
#include "eg_position.h"
#include "uci.h"

class GenEGTB {
    int pieces1[6];
    int pieces2[6];
    int n_pieces;

    int16_t* TB;
    int16_t* MIRROR_TB;

    int16_t* CAPTURE_TBs[6];
    int16_t* MIRROR_CAPTURE_TBs[6];

public:
    GenEGTB(int pieces1[6], int pieces2[6]) {
        n_pieces = 0;
        for (int i = 0; i < 6; i++) {
            this->pieces1[i] = pieces1[i];
            this->pieces2[i] = pieces2[i];
            n_pieces += pieces1[i] + pieces2[i];
        }

        this->TB = (int16_t*) calloc(sizeof(int16_t), num_positions());
        this->MIRROR_TB = (int16_t*) calloc(sizeof(int16_t), num_positions());

        // TODO
        int16_t* KK_TB = (int16_t*) calloc(sizeof(int16_t), N_KKX);
        for (int i = 0; i < 6; i++) {
            CAPTURE_TBs[i] = KK_TB;
            MIRROR_CAPTURE_TBs[i] = KK_TB;
        }

    }

    int num_pieces() const;
    uint64_t num_positions() const;
    void gen();
};

int GenEGTB::num_pieces() const {
    return n_pieces;
}

uint64_t GenEGTB::num_positions() const {
    // TODO
    uint64_t n = N_KKX;
    int k = num_pieces();
    uint64_t s = 62;
    for (int i = 0; i < k; i++) {
        n *= s;
        s--;
    }
    return n;
}




inline int16_t WIN_IN(int16_t level) { return 1000 - level; }
inline int16_t LOSS_IN(int16_t level) { return -1000 + level; }
#define MAYBELOSS -1001
#define UNUSEDIX -1002

void GenEGTB::gen() {
    uint64_t NPOS = num_positions();


    EGPosition pos;

    int16_t LEVEL = 0;

    uint64_t N_UNUSED = 0;
    uint64_t N_LEVEL_POS = 0;
    for (uint64_t ix = 0; ix < NPOS; ix++) {
        pos.reset();
        pos_at_ix(pos, ix, BLACK, pieces1, pieces2);
        if (ix_from_pos(pos) != ix) {
            TB[ix] = UNUSEDIX;
            N_UNUSED++;
            continue;
        }
        if (pos.is_legal_checkmate()) {
            TB[ix] = LOSS_IN(0);
            N_LEVEL_POS++;
        }
    }
    std::cout << "Checkmate count: " << N_LEVEL_POS << std::endl;
    std::cout << N_UNUSED << " indices unused in TB" << std::endl;


    N_UNUSED = 0;
    N_LEVEL_POS = 0;
    for (uint64_t ix = 0; ix < NPOS; ix++) {
        pos.reset();
        pos_at_ix(pos, ix, WHITE, pieces1, pieces2);
        if (ix_from_pos(pos) != ix) {
            MIRROR_TB[ix] = UNUSEDIX;
            N_UNUSED++;
            continue;
        }
        if (pos.is_legal_checkmate()) {
            MIRROR_TB[ix] = LOSS_IN(0);
            N_LEVEL_POS++;
        }
    }
    std::cout << "Mirror Checkmate count: " << N_LEVEL_POS << std::endl;
    std::cout << N_UNUSED << " indices unused in MIRROR_TB" << std::endl;

}

#endif