
#ifndef GEN_EGTB_H_INCLUDED
#define GEN_EGTB_H_INCLUDED

#include "kkx.h"
#include "linearize.h"
#include "eg_position.h"
#include "eg_movegen.h"
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

    uint64_t N_LEVEL_POS = 0;

    int16_t* LOSS_TB;
    int16_t* WIN_TB;
    Color LOSS_COLOR;

    for (int mirrored = 0; mirrored <= 1; ++mirrored) {
        if (mirrored) {
            LOSS_TB = MIRROR_TB;
            LOSS_COLOR = WHITE;
        } else {
            LOSS_TB = TB;
            LOSS_COLOR = BLACK;
        }
        uint64_t N_UNUSED = 0;
        uint64_t N_CHECKMATE = 0;
        for (uint64_t ix = 0; ix < NPOS; ix++) {
            pos.reset();
            pos_at_ix(pos, ix, LOSS_COLOR, pieces1, pieces2);
            if (ix_from_pos(pos) != ix) {
                LOSS_TB[ix] = UNUSEDIX;
                N_UNUSED++;
                continue;
            }
            if (pos.is_legal_checkmate()) {
                LOSS_TB[ix] = LOSS_IN(0);
                N_CHECKMATE++;
                N_LEVEL_POS++;
            }
        }
        std::cout << "Checkmate count: " << N_CHECKMATE << std::endl;
        std::cout << N_UNUSED << " unused indices" << std::endl;
    }

    while (N_LEVEL_POS > 0) {
        LEVEL++;
        N_LEVEL_POS = 0;

        for (int mirrored = 0; mirrored <= 1; ++mirrored) {
            if (mirrored) {
                LOSS_TB = MIRROR_TB;
                LOSS_COLOR = WHITE;
                WIN_TB = TB;
            } else {
                LOSS_TB = TB;
                LOSS_COLOR = BLACK;
                WIN_TB = MIRROR_TB;
            }

            for (uint64_t ix = 0; ix < NPOS; ix++) {

                // for all checkmate in LEVEL-1
                if (LOSS_TB[ix] == LOSS_IN(LEVEL-1)) {
                    pos.reset();
                    pos_at_ix(pos, ix, LOSS_COLOR, pieces1, pieces2);
                    for (Move move : EGMoveList<REVERSE>(pos)) {
                        pos.do_rev_move(move); // non-capture

                        uint64_t win_ix = ix_from_pos(pos);
                        if (WIN_TB[win_ix] == 0) {
                            WIN_TB[win_ix] = WIN_IN(LEVEL);

                            if (N_LEVEL_POS == 0) { std::cout << "WIN in " <<  LEVEL << ": " << pos.fen() << ", ix: " << win_ix << std::endl; }
                            N_LEVEL_POS++;

                            for (Move move2 : EGMoveList<REVERSE>(pos)) {
                                pos.do_rev_move(move2); // non-capture

                                uint64_t maybe_loss_ix = ix_from_pos(pos);
                                if (LOSS_TB[maybe_loss_ix] == 0) {
                                    LOSS_TB[maybe_loss_ix] = MAYBELOSS;
                                }

                                pos.undo_rev_move(move2);
                            }
                        }

                        pos.undo_rev_move(move);
                    }
                }
            }
        }

        std::cout << N_LEVEL_POS << " positions at level " << LEVEL << std::endl;

        LEVEL++;
        N_LEVEL_POS = 0;
        for (int mirrored = 0; mirrored <= 1; ++mirrored) {
            if (mirrored) {
                LOSS_TB = MIRROR_TB;
                LOSS_COLOR = WHITE;
                WIN_TB = TB;
            } else {
                LOSS_TB = TB;
                LOSS_COLOR = BLACK;
                WIN_TB = MIRROR_TB;
            }

            for (uint64_t ix = 0; ix < NPOS; ix++) {
                if (LOSS_TB[ix] == MAYBELOSS) {
                    pos.reset();
                    pos_at_ix(pos, ix, LOSS_COLOR, pieces1, pieces2);

                    // check that all forward moves lead to checkmate in <= -(LEVEL-1)
                    EGMoveList moveList = EGMoveList<FORWARD>(pos);
                    int16_t max_val = LOSS_IN(0);
                    if (moveList.size() == 0) {
                        assert(false);
                    } else {
                        for (Move move : moveList) {
                            Piece capture = pos.do_move(move);
                            if (capture) {
                                max_val = std::max(max_val, (int16_t) 0);
                                // TODO
                            } else {
                                uint64_t fwd_ix = ix_from_pos(pos);
                                int16_t val = -WIN_TB[fwd_ix];
                                max_val = std::max(max_val, val);
                            }
                            pos.undo_move(move, capture);
                        }
                    }
                    if (max_val >= 0) {
                        LOSS_TB[ix] = 0;
                    } else {
                        if (N_LEVEL_POS == 0) { std::cout << "LOSS in " <<  LEVEL << ": " << pos.fen() << ", ix: " << ix << std::endl; }
                        N_LEVEL_POS++;
                        LOSS_TB[ix] = LOSS_IN(LEVEL);;
                    }
                }
            }
        }

        std::cout << N_LEVEL_POS << " positions at level " << LEVEL << std::endl;

        if (N_LEVEL_POS == 0) {
            break;
        }
    }
}

#endif