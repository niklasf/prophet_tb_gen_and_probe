
#ifndef KKP_H_INCLUDED
#define KKP_H_INCLUDED

#include <cstdint>
#include <iostream>
#include "types.h"
#include "bitboard.h"
#include "uci.h"


#define N_KKP 3612
int8_t KKP_KTM_SQ[N_KKP];
int8_t KKP_KNTM_SQ[N_KKP];
int16_t KKP_IX_TABLE[64][64];

// 3612 positions
void _enumerate_kkp(int print_sqs) {
    if (print_sqs == 1) {
        std::cout << "int8_t KKP_KTM_SQ[N_KKP] = {" << std::endl;
    } else if (print_sqs == 2) {
        std::cout << "int8_t KKP_KNTM_SQ[N_KKP] = {" << std::endl;
    }
    int count = 0;
    for (Square ktm_sq = SQ_A1; ktm_sq <= SQ_H8; ++ktm_sq) {
        for (Square kntm_sq = SQ_A1; kntm_sq <= SQ_H8; ++kntm_sq) {
            if (ktm_sq != kntm_sq && (PseudoAttacks[KING][ktm_sq] & kntm_sq) == 0) {
                if (print_sqs == 1) {
                    printf("%2d,", int(ktm_sq));
                } else if (print_sqs == 2) {
                    printf("%2d,", int(kntm_sq));
                }
                KKP_KTM_SQ[count] = ktm_sq;
                KKP_KNTM_SQ[count] = kntm_sq;
                KKP_IX_TABLE[ktm_sq][kntm_sq] = count;
                count++;
            }
        }
        if (print_sqs > 0) { std::cout << std::endl; } 
    }
    if (print_sqs > 0) {
        std::cout << "};" << std::endl;
    } else {
        std::cout << "init_kkp_table count: " << count << std::endl;
    }
}
void enumerate_kkp() {
    std::cout << "#define N_KKP 3612" << std::endl;;
    _enumerate_kkp(1);
    _enumerate_kkp(2);
}
void init_kkp_table() {
    for (Square ktm_sq = SQ_A1; ktm_sq <= SQ_H8; ++ktm_sq) {
        for (Square kntm_sq = SQ_A1; kntm_sq <= SQ_H8; ++kntm_sq) {
            KKP_IX_TABLE[ktm_sq][kntm_sq] = -1;
        }
    }
    _enumerate_kkp(0);
}


inline int16_t get_kkp_ix(Square ktm, Square kntm) {
    int16_t kkx_ix = KKP_IX_TABLE[ktm][kntm];
    if (kkx_ix == -1) {
        printf("Tried to access KKP_IX_TABLE[%d][%d] = -1\n", ktm, kntm);
        assert(false);
    } else {
        return kkx_ix;
    }
}

#endif