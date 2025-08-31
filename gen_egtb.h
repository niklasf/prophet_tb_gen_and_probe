
#ifndef GEN_EGTB_H_INCLUDED
#define GEN_EGTB_H_INCLUDED

#include "kkx.h"
#include "linearize.h"
#include "eg_position.h"
#include "eg_movegen.h"
#include "uci.h"
#include <iostream>
#include <fstream>
#include <string>

uint64_t compute_num_positions(const int stm_pieces[6], const int sntm_pieces[6]) {
    // TODO
    uint64_t n = N_KKX;
    
    int k = 0;
    for (int i = 0; i < 6; i++) {
        k += stm_pieces[i] + sntm_pieces[i];
    }

    uint64_t s = 62;
    for (int i = 0; i < k; i++) {
        n *= s;
        s--;
    }
    return n;
}

std::string get_egtb_identifier(int stm_pieces[6], int sntm_pieces[6]) {
    std::ostringstream os;
    for (int* pieces: {stm_pieces, sntm_pieces}) {
        os << "K";
        for (PieceType pt = QUEEN; pt >= PAWN; --pt) {
            for (int i = 0; i < pieces[pt]; i++) {
                os << PieceToChar[pt];
            }
        }
    }
    return os.str();
}

std::string get_filename(int stm_pieces[6], int sntm_pieces[6]) {
    std::ostringstream os;
    os << "egtbs/";
    os << get_egtb_identifier(stm_pieces, sntm_pieces);
    os << ".egtb";
    return os.str();
}


void store_egtb(int16_t* TB, int stm_pieces[6], int sntm_pieces[6]) {
    uint64_t NPOS = compute_num_positions(stm_pieces, sntm_pieces);
    std::string filename = get_filename(stm_pieces, sntm_pieces);
    std::ofstream outputFileStream;
    outputFileStream.open(filename, std::ios::out|std::ios::binary);
    for(uint64_t i=0; i<NPOS; i++)
        outputFileStream.write((char*) &TB[i], sizeof(int16_t));
}

void load_egtb(int16_t* TB, int stm_pieces[6], int sntm_pieces[6]) {
    uint64_t NPOS = compute_num_positions(stm_pieces, sntm_pieces);
    std::string filename = get_filename(stm_pieces, sntm_pieces);
    std::ifstream inputFileStream;
    inputFileStream.open(filename, std::ios::in|std::ios::binary);
    for(uint64_t i=0; i<NPOS; i++)
        inputFileStream.read((char*) &TB[i], sizeof(int16_t));
}


class GenEGTB {
    int wpieces[6];
    int bpieces[6];
    int n_pieces;

    int16_t* WTM_TB;
    int16_t* BTM_TB;

    int16_t* WTM_CAPTURE_TBs[6];
    int16_t* BTM_CAPTURE_TBs[6];

public:
    GenEGTB(int wpieces[6], int bpieces[6]) {
        this->n_pieces = 0;
        for (int i = 0; i < 6; i++) {
            this->wpieces[i] = wpieces[i];
            this->bpieces[i] = bpieces[i];
            this->n_pieces += wpieces[i] + bpieces[i];
        }

        this->WTM_TB = (int16_t*) calloc(sizeof(int16_t), num_positions());
        this->BTM_TB = (int16_t*) calloc(sizeof(int16_t), num_positions());

        for (PieceType pt = PAWN; pt <= QUEEN; ++pt) {
            if (wpieces[pt] > 0) {
                wpieces[pt]--;
                this->BTM_CAPTURE_TBs[pt] = (int16_t*) calloc(sizeof(int16_t), compute_num_positions(bpieces, wpieces));
                load_egtb(this->BTM_CAPTURE_TBs[pt], bpieces, wpieces);
                std::cout << "Load " << get_egtb_identifier(bpieces, wpieces) << std::endl;
                wpieces[pt]++;
            }
            if (bpieces[pt] > 0) {
                bpieces[pt]--;
                this->WTM_CAPTURE_TBs[pt] = (int16_t*) calloc(sizeof(int16_t), compute_num_positions(wpieces, bpieces));
                load_egtb(this->WTM_CAPTURE_TBs[pt], wpieces, bpieces); 
                std::cout << "Load " << get_egtb_identifier(wpieces, bpieces) << std::endl;           
                bpieces[pt]++;
            }
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
    return compute_num_positions(wpieces, bpieces);
}



inline int16_t WIN_IN(int16_t level) { return 1000 - level; }
inline int16_t LOSS_IN(int16_t level) { return -1000 + level; }
#define MAYBELOSS -1001
#define UNUSEDIX -1002

void GenEGTB::gen() {
    std::cout << "Generate " << get_egtb_identifier(wpieces, bpieces) << " and " << get_egtb_identifier(bpieces, wpieces) << std::endl;

    uint64_t NPOS = num_positions();


    EGPosition pos;

    int16_t LEVEL = 0;

    uint64_t N_LEVEL_POS = 0;

    int16_t* LOSS_TB;
    int16_t* WIN_TB;
    int16_t** CAPTURE_TBs;
    Color LOSS_COLOR;

    for (int wtm = 0; wtm <= 1; ++wtm) {
        if (wtm) {
            LOSS_TB = WTM_TB;
            LOSS_COLOR = WHITE;
        } else {
            LOSS_TB = BTM_TB;
            LOSS_COLOR = BLACK;
        }
        uint64_t N_UNUSED = 0;
        uint64_t N_CHECKMATE = 0;
        for (uint64_t ix = 0; ix < NPOS; ix++) {
            pos.reset();
            pos_at_ix(pos, ix, LOSS_COLOR, wpieces, bpieces);
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

        for (int wtm = 0; wtm <= 1; ++wtm) {
            if (wtm) {
                LOSS_TB = WTM_TB;
                LOSS_COLOR = WHITE;
                WIN_TB = BTM_TB;
            } else {
                LOSS_TB = BTM_TB;
                LOSS_COLOR = BLACK;
                WIN_TB = WTM_TB;
            }

            for (uint64_t ix = 0; ix < NPOS; ix++) {

                // for all checkmate in LEVEL-1
                if (LOSS_TB[ix] == LOSS_IN(LEVEL-1)) {
                    pos.reset();
                    pos_at_ix(pos, ix, LOSS_COLOR, wpieces, bpieces);
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
        for (int wtm = 0; wtm <= 1; ++wtm) {
            if (wtm) {
                LOSS_TB = WTM_TB;
                LOSS_COLOR = WHITE;
                WIN_TB = BTM_TB;
                CAPTURE_TBs = WTM_CAPTURE_TBs;
            } else {
                LOSS_TB = BTM_TB;
                LOSS_COLOR = BLACK;
                WIN_TB = WTM_TB;
                CAPTURE_TBs = BTM_CAPTURE_TBs;
            }

            for (uint64_t ix = 0; ix < NPOS; ix++) {
                if (LOSS_TB[ix] == MAYBELOSS) {
                    pos.reset();
                    pos_at_ix(pos, ix, LOSS_COLOR, wpieces, bpieces);

                    // check that all forward moves lead to checkmate in <= -(LEVEL-1)
                    EGMoveList moveList = EGMoveList<FORWARD>(pos);
                    int16_t max_val = LOSS_IN(0);
                    if (moveList.size() == 0) {
                        assert(false);
                    } else {
                        for (Move move : moveList) {
                            Piece capture = pos.do_move(move);
                            if (capture) {
                                // max_val = std::max(max_val, (int16_t) 0);
                                uint64_t fwd_ix = ix_from_pos(pos);
                                int16_t val = -CAPTURE_TBs[type_of(capture)][fwd_ix];
                                max_val = std::max(max_val, val);
                            } else {
                                uint64_t fwd_ix = ix_from_pos(pos);
                                int16_t val = -WIN_TB[fwd_ix];
                                max_val = std::max(max_val, val);
                            }
                            pos.undo_move(move, capture);
                        }
                    }
                    if (max_val >= 0) {
                        // assert(max_val == 0);
                        LOSS_TB[ix] = 0;
                    } else {
                        if (N_LEVEL_POS == 0) { std::cout << "LOSS in " <<  LEVEL << ": " << pos.fen() << ", ix: " << ix << std::endl; }
                        N_LEVEL_POS++;
                        // assert(max_val == LOSS_IN(LEVEL-1));
                        // LOSS_TB[ix] = LOSS_IN(LEVEL);
                        LOSS_TB[ix] = max_val + 1;
                    }

                }
            }
        }

        std::cout << N_LEVEL_POS << " positions at level " << LEVEL << std::endl;

        if (N_LEVEL_POS == 0) {
            break;
        }
    }
    for (int wtm = 0; wtm <= 1; ++wtm) {
        int16_t* _TB = wtm ? WTM_TB : BTM_TB;
        std::string egtb = wtm ? get_egtb_identifier(wpieces, bpieces) : get_egtb_identifier(bpieces, wpieces);

        uint64_t wins = 0;
        uint64_t draws = 0;
        uint64_t losses = 0;
        for (uint64_t ix = 0; ix < NPOS; ix++) {
            if (_TB[ix] == UNUSEDIX) { continue; }
            wins += (_TB[ix] > 0);
            draws += (_TB[ix] == 0);
            losses += (_TB[ix] < 0);
        }
        std::cout << wins << " wins in " << egtb << std::endl;
        std::cout << draws << " draws in " << egtb << std::endl;
        std::cout << losses << " losses in " << egtb << std::endl;
    }

    store_egtb(WTM_TB, wpieces, bpieces);
    store_egtb(BTM_TB, bpieces, wpieces);
    std::cout << "\n\n" << std::endl;
}

#endif