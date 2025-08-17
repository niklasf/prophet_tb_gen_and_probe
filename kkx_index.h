
#ifndef KKX_INDEX_H_INCLUDED
#define KKX_INDEX_H_INCLUDED

#include <vector>

#include "kkx.h"
#include "eg_position.h"

class KKXIndex {
private:
    std::vector<PieceType> stm_pieces;
    std::vector<PieceType> sntm_pieces;
public:
    KKXIndex(std::vector<PieceType>& stm_pieces, std::vector<PieceType>& sntm_pieces) {
        this->stm_pieces = stm_pieces;
        this->sntm_pieces = sntm_pieces;
    }

    int n_pieces() const;
    uint64_t num_positions() const;

    void pos_at_ix(EGPosition &pos, uint64_t ix);
};

int KKXIndex::n_pieces() const {
    return stm_pieces.size() + sntm_pieces.size();
}

uint64_t KKXIndex::num_positions() const {
    uint64_t n = N_KKX;
    int k = n_pieces();
    uint64_t s = 62;
    for (int i = 0; i < k; i++) {
        n *= s;
        s--;
    }
    return n;
}

void KKXIndex::pos_at_ix(EGPosition &pos, uint64_t ix) {
    uint64_t kkx = ix % N_KKX;
    ix = ix / N_KKX;
    pos.set_side_to_move(BLACK);


    Square ktm_sq = Square(KKX_KTM_SQ[kkx]);
    Square kntm_sq = Square(KKX_KNTM_SQ[kkx]);
    // occupied_sqs[0] = ktm_sq;
    // occupied_sqs[1] = kntm_sq;
    pos.put_piece(B_KING, ktm_sq);
    pos.put_piece(W_KING, Square(kntm_sq));

    uint64_t s = 62;
    for (PieceType p : stm_pieces) {
        Square sq = Square(ix % s);
        ix = ix / s;
        Bitboard occupied = pos.pieces();
        while (square_bb(sq) & occupied) ++sq;
        pos.put_piece(make_piece(BLACK, p), sq);
        s--;
    }

    for (PieceType p : sntm_pieces) {
        Square sq = Square(ix % s);
        ix = ix / s;
        Bitboard occupied = pos.pieces();
        while (square_bb(sq) & occupied) ++sq;
        pos.put_piece(make_piece(WHITE, p), sq);
        s--;
    }
}

#endif