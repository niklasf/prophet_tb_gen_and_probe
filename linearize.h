
#ifndef LINEARIZE_H_INCLUDED
#define LINEARIZE_H_INCLUDED

#include "kkx.h"
#include "eg_position.h"
#include "uci.h"

void pos_at_ix_kkx(EGPosition &pos, uint64_t ix, Color stm, int wpieces[6], int bpieces[6]) {
    assert (wpieces[PAWN] + bpieces[PAWN] == 0);
    pos.set_side_to_move(stm);

    bool allondiag = true;
    Square occupied_sqs[6];

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

    Piece pieces[4] = {NO_PIECE,NO_PIECE,NO_PIECE,NO_PIECE};
    int piece_count = 0;
    for (Color c: {~stm, stm}) {
        for (PieceType p : {QUEEN, ROOK, BISHOP, KNIGHT}) {
            int* c_pieces = (c == WHITE) ? wpieces : bpieces;
            if (c_pieces[p] > 1) { std::cout << "Mutliple pieces of same type not supported!\n"; assert(false); }
            if (piece_count >= 4) { std::cout << "More than 6 pieces not supported!\n"; assert(false); }
            if (c_pieces[p] == 0) { continue; }
            pieces[piece_count] = make_piece(c, p);
            piece_count++;
        }
    }

    uint64_t n_available_squares = 62;
    int i = 2;

    // this also generates position where first off-diagonal piece is not on bottom half
    // but these indexes are unused anyways
    for (Piece p : pieces) {
        if (p == NO_PIECE) { break; }
        // this assumes no two pieces of same kind
        Square orig_sq = Square(ix % n_available_squares);
        ix = ix / n_available_squares;
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
        i++;
        occupied_sqs[k] = sq;

        pos.put_piece(p, sq);
        n_available_squares--;
    }
}

void pos_at_ix_kkp(EGPosition &pos, uint64_t ix, Color stm, int wpieces[6], int bpieces[6]) {
    assert (wpieces[PAWN] + bpieces[PAWN] > 0);
    pos.set_side_to_move(stm);

    Square occupied_sqs[6];

    Piece pieces[6] = {NO_PIECE,NO_PIECE,NO_PIECE,NO_PIECE,NO_PIECE,NO_PIECE};
    int piece_count = 0;
    for (Color c: {~stm, stm}) {
        PieceType p = PAWN;
        int* c_pieces = (c == WHITE) ? wpieces : bpieces;
        if (c_pieces[p] > 1) { std::cout << "Mutliple pieces of same type not supported!\n"; assert(false); }
        if (piece_count >= 6) { std::cout << "More than 6 pieces not supported!\n"; assert(false); }
        if (c_pieces[p] == 0) { continue; }
        pieces[piece_count] = make_piece(c, p);
        piece_count++;
    }
    for (Color c: {~stm, stm}) {
        pieces[piece_count] = make_piece(c, KING);
        piece_count++;
        for (PieceType p : {QUEEN, ROOK, BISHOP, KNIGHT}) {
            int* c_pieces = (c == WHITE) ? wpieces : bpieces;
            if (c_pieces[p] > 1) { std::cout << "Mutliple pieces of same type not supported!\n"; assert(false); }
            if (piece_count >= 6) { std::cout << "More than 6 pieces not supported!\n"; assert(false); }
            if (c_pieces[p] == 0) { continue; }
            pieces[piece_count] = make_piece(c, p);
            piece_count++;
        }
    }


    uint64_t n_available_pawn_squares = 48;
    uint64_t n_available_squares = 64;
    int i = 0;

    for (Piece p : pieces) {
        if (p == NO_PIECE) { break; }
        // this assumes no two pieces of same kind
        Square orig_sq;
        if (i == 0) {
            // first pawn
            uint64_t pix = ix % 24;
            orig_sq = Square(pix + (pix >> 2) * 4 + 8); // put on left side of board
            ix = ix / 24;
        } else if (type_of(p) == PAWN) {
            orig_sq = Square(ix % n_available_pawn_squares + 8);
            ix = ix / n_available_pawn_squares;
        } else {
            orig_sq = Square(ix % n_available_squares);
            ix = ix / n_available_squares;
        }
        
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
        i++;
        occupied_sqs[k] = sq;

        pos.put_piece(p, sq);
        n_available_squares--;
        n_available_pawn_squares--;
    }
}

void pos_at_ix(EGPosition &pos, uint64_t ix, Color stm, int wpieces[6], int bpieces[6]) {
    if (wpieces[PAWN] + bpieces[PAWN] > 0) {
        pos_at_ix_kkp(pos, ix, stm, wpieces, bpieces);
    } else {
        pos_at_ix_kkx(pos, ix, stm, wpieces, bpieces);
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

uint64_t ix_from_pos_kkx(EGPosition const &pos) {
    assert (pos.count<PAWN>(WHITE) + pos.count<PAWN>(BLACK) == 0);

    Color stm = pos.side_to_move();

    bool allondiag = true;
    Square occupied_sqs[6];

    Square orig_ktm_sq = pos.square<KING>(stm);

    int8_t flip = ((orig_ktm_sq & RightHalfBB) ? 7 : 0) ^ ((orig_ktm_sq & TopHalfBB) ? 56 : 0);
    int8_t swap = 0;
    maybe_update_swap(orig_ktm_sq, flip, allondiag, swap);
    Square ktm_sq = transform(orig_ktm_sq, flip, swap);


    Square orig_kntm_sq = pos.square<KING>(~stm);
    maybe_update_swap(orig_kntm_sq, flip, allondiag, swap);
    Square kntm_sq = transform(orig_kntm_sq, flip, swap);

    int16_t kkx_ix = get_kkx_ix(orig_ktm_sq, orig_kntm_sq);

    uint64_t ix = kkx_ix;
    uint64_t multiplier = N_KKX;

    if (ktm_sq < kntm_sq) {
        occupied_sqs[0] = ktm_sq;
        occupied_sqs[1] = kntm_sq;
    } else {

        occupied_sqs[0] = kntm_sq;
        occupied_sqs[1] = ktm_sq;
    }


    uint64_t n_available_squares = 62;
    int i = 2;

    for (Color c: {~stm, stm}) {
        for (PieceType p: {QUEEN, ROOK, BISHOP, KNIGHT}) {
            Bitboard pieceBB = pos.pieces(c, p);
            assert(!more_than_one(pieceBB));
            if (pieceBB) {
                Square orig_sq = lsb(pieceBB);
                maybe_update_swap(orig_sq, flip, allondiag, swap);
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
                multiplier *= n_available_squares;
                n_available_squares--;
            }
        }
    }
    return ix;
}

uint64_t ix_from_pos_kkp(EGPosition const &pos) {
    assert (pos.count<PAWN>(WHITE) + pos.count<PAWN>(BLACK) > 0);
    Color stm = pos.side_to_move();

    Square occupied_sqs[6];

    int8_t flip  = 0;

    // first pawn constrained to left side of board, no vertical symmetry
    for (Color c: {~stm, stm}) {
        if (pos.count<PAWN>(c) > 0) {
            flip = pos.square<PAWN>(c) & RightHalfBB ? 7 : 0;
            break;
        }
    }

    uint64_t ix = 0;
    uint64_t multiplier = 1;

    bool first_pawn = true;
    uint64_t n_available_pawn_squares = 48;
    uint64_t n_available_squares = 64;
    int i = 0;
    for (Color c: {~stm, stm}) {
        Bitboard pieceBB = pos.pieces(c, PAWN);
        assert(!more_than_one(pieceBB));
        if (pieceBB) {
            Square sq = transform(lsb(pieceBB), flip, 0);

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

            if (first_pawn) {
                first_pawn = false;
                uint64_t sq_ix = (uint64_t) sq - 8;
                sq_ix = sq_ix - (sq_ix >> 3) * 4; // map to 0,1,...,23
                ix += sq_ix * multiplier;
                multiplier *= 24;
            } else {
                ix += ((uint64_t) sq - 8) * multiplier;
                multiplier *= n_available_pawn_squares;
            }
            n_available_pawn_squares--;
            n_available_squares--;
        }
    }


    for (Color c: {~stm, stm}) {
        for (PieceType p: {KING, QUEEN, ROOK, BISHOP, KNIGHT}) {
            Bitboard pieceBB = pos.pieces(c, p);
            assert(!more_than_one(pieceBB));
            if (pieceBB) {
                Square sq = transform(lsb(pieceBB), flip, 0);

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
                multiplier *= n_available_squares;
                n_available_squares--;
            }
        }
    }
    return ix;
}

uint64_t ix_from_pos(EGPosition const &pos) {
    if (pos.count<PAWN>(WHITE) + pos.count<PAWN>(BLACK) > 0) {
        return ix_from_pos_kkp(pos);
    } else {
        return ix_from_pos_kkx(pos);
    }
}

void transform_to_canoncial(const EGPosition &pos, EGPosition &pos2) {
    Color stm = pos.side_to_move();

    bool allondiag = true;
    int8_t flip = 0;
    int8_t swap = 0;

    Square orig_ktm_sq = pos.square<KING>(stm);

    if (pos.count<PAWN>(WHITE) + pos.count<PAWN>(BLACK) == 0) {
        flip = ((orig_ktm_sq & RightHalfBB) ? 7 : 0) ^ ((orig_ktm_sq & TopHalfBB) ? 56 : 0);
    } else {
        for (Color c: {~stm, stm}) {
            if (pos.count<PAWN>(c) > 0) {
                flip = pos.square<PAWN>(c) & RightHalfBB ? 7 : 0;
                break;
            }
        }
        allondiag = false; // disable diagonal symmetries
    }

    maybe_update_swap(orig_ktm_sq, flip, allondiag, swap);
    pos2.put_piece(make_piece(stm, KING), transform(orig_ktm_sq, flip, swap));

    Square orig_kntm_sq = pos.square<KING>(~stm);
    maybe_update_swap(orig_kntm_sq, flip, allondiag, swap);
    pos2.put_piece(make_piece(~stm, KING), transform(orig_kntm_sq, flip, swap));

    for (Color c: {~stm, stm}) {
        for (PieceType pt: {QUEEN, ROOK, BISHOP, KNIGHT, PAWN}) {
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

void transform_to(const EGPosition &pos, EGPosition &pos2, int8_t h_flip, int8_t v_flip, int8_t swap) {

    for (Square sq = SQ_A1; sq <= SQ_H8; ++sq) {
        Piece p = pos.piece_on(sq);
        if (p)
            pos2.put_piece(p, transform(sq, h_flip ^ v_flip, swap));
    }

    pos2.set_side_to_move(pos.side_to_move());

}

#endif