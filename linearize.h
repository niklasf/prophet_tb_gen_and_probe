
#ifndef LINEARIZE_H_INCLUDED
#define LINEARIZE_H_INCLUDED

#include "kkx.h"
#include "eg_position.h"
#include "uci.h"

void pos_at_ix(EGPosition &pos, uint64_t ix, Color stm, int wpieces[6], int bpieces[6]) {
    pos.set_side_to_move(stm);

    bool allondiag = true;
    Square occupied_sqs[6];

    if (wpieces[PAWN] + bpieces[PAWN] == 0) {
        uint64_t kkx = ix % N_KKX;
        ix = ix / N_KKX;

        Square ktm_sq = Square(KKX_KTM_SQ[kkx]);
        allondiag = allondiag && (ktm_sq & DiagBB);
        Square kntm_sq = Square(KKX_KNTM_SQ[kkx]);
        allondiag = allondiag && (kntm_sq & DiagBB);
        if (ktm_sq < kntm_sq) {
            occupied_sqs[0] = ktm_sq;
            occupied_sqs[1] = kntm_sq;
        } else {
            occupied_sqs[0] = kntm_sq;
            occupied_sqs[1] = ktm_sq;
        }
        pos.put_piece(make_piece(stm, KING), ktm_sq);
        pos.put_piece(make_piece(~stm, KING), kntm_sq);
    } else {
        std::cout << "Pawns not supported!\n";
        assert(false);
    }

    Piece pieces[4] = {NO_PIECE,NO_PIECE,NO_PIECE,NO_PIECE};
    int piece_count = 0;
    for (Color c: {~stm, stm}) {
        for (PieceType p : {QUEEN, ROOK, BISHOP, KNIGHT, PAWN}) {
            int* c_pieces = (c == WHITE) ? wpieces : bpieces;
            if (c_pieces[p] > 1) { std::cout << "Mutliiple pieces of same type not supported!\n"; assert(false); }
            if (piece_count >= 4) { std::cout << "More than 4 pieces not supported!\n"; assert(false); }
            if (c_pieces[p] == 0) { continue; }
            pieces[piece_count] = make_piece(c, p);
            piece_count++;
        }
    }

    uint64_t s = 62;
    int i = 2;

    for (Piece p : pieces) {
        if (p == NO_PIECE) { break; }
        // this assumes no two pieces of same kind
        Square orig_sq = Square(ix % s);
        ix = ix / s;
        Square sq = orig_sq;

        int k = 0;
        for (int j = 0; j < i; j++) {
            if (occupied_sqs[j] <= sq) {
                ++sq;
                k = j + 1;
            }
        }
        for (int j = i-1; j >= k; j--) {
            occupied_sqs[j+1] = occupied_sqs[j];
        }
        
        if (allondiag && !(sq & DiagBB)) {
            allondiag = false;
            if (sq & AboveDiagBB) {
                // put sq on bottom diag
                sq = Square(((sq >> 3) | (sq << 3)) & 63);
            }
        }

        i++;
        occupied_sqs[k] = sq;

        pos.put_piece(p, sq);
        s--;
    }
}

inline void maybe_update_swap(Square sq, int8_t flip, bool& allondiag, int8_t& swap) {
    // if this changes swap, we do not need to swap previous pieces since they are all on diagonal anyways
    if (allondiag) {
        if (!((sq ^ flip) & DiagBB)) {
            allondiag = false;
            swap = ((sq ^ flip) & AboveDiagBB) ? 3 : 0;
        }
    }
}

inline Square transform(const Square sq, int8_t flip, int8_t swap) {
    int8_t sq_ix = int8_t(sq) ^ flip;
    return Square(((sq_ix >> swap) | (sq_ix << swap)) & 63);
}

uint64_t ix_from_pos(EGPosition const &pos) {
    Color stm = pos.side_to_move();

    bool allondiag = true;
    Square occupied_sqs[6];


    uint64_t ix;
    uint64_t multiplier;
    int8_t flip  = 0;
    int8_t swap = 0;
    if (pos.count<PAWN>(WHITE) + pos.count<PAWN>(BLACK) == 0) {

        Square orig_ktm_sq = pos.square<KING>(stm);

        flip = ((orig_ktm_sq & RightHalfBB) ? 7 : 0) ^ ((orig_ktm_sq & TopHalfBB) ? 56 : 0);
        swap = 0;
        maybe_update_swap(orig_ktm_sq, flip, allondiag, swap);
        Square ktm_sq = transform(orig_ktm_sq, flip, swap);


        Square orig_kntm_sq = pos.square<KING>(~stm);
        maybe_update_swap(orig_kntm_sq, flip, allondiag, swap);
        Square kntm_sq = transform(orig_kntm_sq, flip, swap);

        int16_t kkx_ix = get_kkx_ix(orig_ktm_sq, orig_kntm_sq);


        if (ktm_sq < kntm_sq) {
            occupied_sqs[0] = ktm_sq;
            occupied_sqs[1] = kntm_sq;
        } else {

            occupied_sqs[0] = kntm_sq;
            occupied_sqs[1] = ktm_sq;
        }

        ix = kkx_ix;
        multiplier = N_KKX;
    
        
    } else {
        std::cout << "Pawns not supported!\n";
        assert(false);
    }


    uint64_t s = 62;
    int i = 2;


    for (Color c: {~stm, stm}) {
        for (PieceType p: {QUEEN, ROOK, BISHOP, KNIGHT}) {
            Bitboard pieceBB = pos.pieces(c, p);
            assert(!more_than_one(pieceBB));
            if (pieceBB) {
                Square orig_sq = lsb(pieceBB);
                maybe_update_swap(orig_sq, flip, allondiag, swap); // TODO: do not update swap when we have pawns
                Square sq = transform(orig_sq, flip, swap);

                Square before_sq = sq;

                int k = i;
                for (int j = i-1; j >= 0; j--) {
                    if (occupied_sqs[j] < sq) {
                        --sq;
                    } else {
                        k = j;
                        occupied_sqs[j+1] = occupied_sqs[j];
                    }
                }
                

                occupied_sqs[k] = before_sq;
                i++;


                ix += ((uint64_t) sq) * multiplier;
                multiplier *= s;

                s--;
            }
        }
    }
    return ix;
}


void transform_to(const EGPosition &pos, EGPosition &pos2) {

    if (pos.count<PAWN>(WHITE) + pos.count<PAWN>(BLACK) > 0) {
        std::cout << "Pawns not supported!\n";
        assert(false);
    }

    Color stm = pos.side_to_move();
    bool allondiag = true;
    Square orig_ktm_sq = pos.square<KING>(stm);
    int8_t flip = ((orig_ktm_sq & RightHalfBB) ? 7 : 0) ^ ((orig_ktm_sq & TopHalfBB) ? 56 : 0);
    int8_t swap = 0;
    maybe_update_swap(orig_ktm_sq, flip, allondiag, swap);
    pos2.put_piece(make_piece(stm, KING), transform(orig_ktm_sq, flip, swap));

    Square orig_kntm_sq = pos.square<KING>(~stm);
    maybe_update_swap(orig_kntm_sq, flip, allondiag, swap);
    pos2.put_piece(make_piece(~stm, KING), transform(orig_kntm_sq, flip, swap));

    for (PieceType pt: {QUEEN, ROOK, BISHOP, KNIGHT}) {
        for (Color c: {~stm, stm}) {
            Bitboard bb = pos.pieces(c, pt);
            if (bb) {
                Square sq = lsb(bb);
                maybe_update_swap(sq, flip, allondiag, swap);
                pos2.put_piece(make_piece(c,pt), transform(sq, flip, swap));
            }
        }
    }
    pos2.set_side_to_move(pos.side_to_move());
}


#endif