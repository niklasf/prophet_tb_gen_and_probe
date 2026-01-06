#ifndef PROPHET_H
#define PROPHET_H

#define NO_OMP

#include "types.h"
#include "values.h"
#include <string>
#include "compressed_tb.h"

#ifndef TYPES_H_INCLUDED
enum Piece : std::int8_t {
    NO_PIECE,
    W_PAWN = PAWN,     W_KNIGHT, W_BISHOP, W_ROOK, W_QUEEN, W_KING,
    B_PAWN = PAWN + 8, B_KNIGHT, B_BISHOP, B_ROOK, B_QUEEN, B_KING,
    PIECE_NB = 16
};
enum Square : int8_t {
    SQ_A1, SQ_B1, SQ_C1, SQ_D1, SQ_E1, SQ_F1, SQ_G1, SQ_H1,
    SQ_A2, SQ_B2, SQ_C2, SQ_D2, SQ_E2, SQ_F2, SQ_G2, SQ_H2,
    SQ_A3, SQ_B3, SQ_C3, SQ_D3, SQ_E3, SQ_F3, SQ_G3, SQ_H3,
    SQ_A4, SQ_B4, SQ_C4, SQ_D4, SQ_E4, SQ_F4, SQ_G4, SQ_H4,
    SQ_A5, SQ_B5, SQ_C5, SQ_D5, SQ_E5, SQ_F5, SQ_G5, SQ_H5,
    SQ_A6, SQ_B6, SQ_C6, SQ_D6, SQ_E6, SQ_F6, SQ_G6, SQ_H6,
    SQ_A7, SQ_B7, SQ_C7, SQ_D7, SQ_E7, SQ_F7, SQ_G7, SQ_H7,
    SQ_A8, SQ_B8, SQ_C8, SQ_D8, SQ_E8, SQ_F8, SQ_G8, SQ_H8,
    SQ_NONE,

    SQUARE_ZERO = 0,
    SQUARE_NB   = 64
};
#endif

void init_prophet_tb(std::string path);
void deinit_prophet_tb();

// returns v=0 if draw or illegal position
// returns  1000-v if win in v plies
// returns -1000+v if loss in v plies
int16_t probe(Piece pieces[6], Square squares[6]);

// for better performance reuse DecompressCtx
// for multi-threading use one DecompressCtx per thread
DecompressCtx* CreateDecompressCtx();
int16_t probe_dctx(Piece pieces[6], Square squares[6], DecompressCtx* dctx);

int dtm(int16_t v) {
    if (v == 0) {
        return 0;
    } else if (v > 0) {
        return WIN_IN(0) - v;
    } else {
        return v - LOSS_IN(0);
    }
}
#endif