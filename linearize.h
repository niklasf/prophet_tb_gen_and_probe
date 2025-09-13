
#ifndef LINEARIZE_H_INCLUDED
#define LINEARIZE_H_INCLUDED

#include "kkx.h"
#include "eg_position.h"
#include "uci.h"
#include "triangular_indexes.h"

void insert_increment_sq(Square occupied_sqs[6], int& n_occupied_sqs, Square& sq) {
    int k = 0;
    for (int j = 0; j < n_occupied_sqs; j++) {
        if (occupied_sqs[j] <= sq) {
            ++sq;
            k = j + 1;
        }
    }
    for (int j = n_occupied_sqs-1; j >= k; j--) {
        occupied_sqs[j+1] = occupied_sqs[j];
    }
    n_occupied_sqs++;
    occupied_sqs[k] = sq;
}

int insert_count_lt_squares(Square occupied_sqs[6], int& n_occupied_sqs, Square sq) {
    int k = 0;
    for (int j = 0; j < n_occupied_sqs; j++) {
        k += (occupied_sqs[j] < sq);
    }
    occupied_sqs[n_occupied_sqs] = sq;
    n_occupied_sqs++;
    return k;
}

void pos_at_ix_kkx(EGPosition &pos, uint64_t ix, Color stm, int wpieces[6], int bpieces[6]) {
    assert (wpieces[PAWN] + bpieces[PAWN] == 0);
    pos.set_side_to_move(stm);

    bool is_diag_symmetric = true;
    Square occupied_sqs[6];

    uint64_t kkx = ix % N_KKX;
    ix = ix / N_KKX;

    Square ktm_sq = Square(KKX_KTM_SQ[kkx]);
    is_diag_symmetric = is_diag_symmetric && (ktm_sq & DiagBB);
    Square kntm_sq = Square(KKX_KNTM_SQ[kkx]);
    is_diag_symmetric = is_diag_symmetric && (kntm_sq & DiagBB);

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
    int piece_counts[4] = {0,0,0,0};
    int total_piece_count = 0;
    int i = 0;
    for (Color c: {~stm, stm}) {
        int* c_pieces = (c == WHITE) ? wpieces : bpieces;
        for (PieceType p : {QUEEN, ROOK, BISHOP, KNIGHT}) {
            if (c_pieces[p] == 0) { continue; }
            pieces[i] = make_piece(c, p);
            piece_counts[i] = c_pieces[p];
            i++;
            total_piece_count += piece_counts[i];
            if (total_piece_count >= 4) { std::cout << "More than 6 pieces not supported!\n"; assert(false); }
        }
    }


    int n_occupied_sqs = 2;

    int sqs_ixs[4];

    // this also generates position where first off-diagonal piece is not on bottom half
    // but these indexes are unused anyways
    for (int l = 0; l < 4; l++) {
        Piece p = pieces[l];
        int piece_count = piece_counts[l];
        if (p == NO_PIECE) { break; }
        uint64_t s = number_of_ordered_tuples(64 - n_occupied_sqs, piece_count);

        uint64_t tril_ix = ix % s;
        ix = ix / s;

        tril_from_linear(piece_count, tril_ix, sqs_ixs);
        for (int j = 0; j < piece_count; j++) {
            Square sq = Square(sqs_ixs[j]-j);

            insert_increment_sq(occupied_sqs, n_occupied_sqs, sq);

            pos.put_piece(p, sq);
        }
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
    int n_occupied_sqs = 0;

    for (Piece p : pieces) {
        if (p == NO_PIECE) { break; }
        // this assumes no two pieces of same kind
        Square sq;
        if (n_occupied_sqs == 0) {
            // first pawn
            uint64_t pix = ix % 24;
            sq = Square(pix + (pix >> 2) * 4 + 8); // put on left side of board
            ix = ix / 24;
        } else if (type_of(p) == PAWN) {
            sq = Square(ix % n_available_pawn_squares + 8);
            ix = ix / n_available_pawn_squares;
        } else {
            sq = Square(ix % n_available_squares);
            ix = ix / n_available_squares;
        }

        insert_increment_sq(occupied_sqs, n_occupied_sqs, sq);

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

inline Square transform(const Square sq, int8_t flip, int8_t swap) {
    int8_t sq_ix = int8_t(sq) ^ flip;
    return Square(((sq_ix >> swap) | (sq_ix << swap)) & 63);
}

inline Square maybe_update_swap_and_transform(Square sq, int8_t flip, bool& is_diag_symmetric, int8_t& swap) {
    // if this changes swap, we do not need to swap previous pieces since they are all on diagonal anyways
    if (is_diag_symmetric) {
        if (!((sq ^ flip) & DiagBB)) {
            is_diag_symmetric = false;
            swap = ((sq ^ flip) & AboveDiagBB) ? 3 : 0;
        }
    }
    return transform(sq, flip, swap);
}

inline Bitboard maybe_update_swap_and_transform_bb(Bitboard piecesBB, int8_t flip, bool& is_diag_symmetric, int8_t& swap) {
    Bitboard b = 0;
    Bitboard flipped_b = 0;
    while (piecesBB) {
        Square sq = pop_lsb(piecesBB);
        b |= square_bb(transform(sq, flip, 0));
        flipped_b |= square_bb(transform(sq, flip, 3));
    }

    if (is_diag_symmetric) {

        if (flipped_b != b) {
            is_diag_symmetric = false;
            Bitboard lower = (b & BelowDiagBB) & ~(flipped_b & BelowDiagBB);
            Bitboard lower_flipped = (flipped_b & BelowDiagBB) & ~(b & BelowDiagBB);
            if (lower_flipped == 0) {
                swap = 0;
            } else if (lower == 0) {
                swap = 3;
            } else {
                swap = lsb(lower_flipped) < lsb(lower) ? 3 : 0;
            }
        }
    }
    
    if (!swap) {
        return b;
    } else {
        return flipped_b;
    }
}


uint64_t ix_from_pos_kkx(EGPosition const &pos) {
    assert (pos.count<PAWN>(WHITE) + pos.count<PAWN>(BLACK) == 0);

    Color stm = pos.side_to_move();

    bool is_diag_symmetric = true;
    Square occupied_sqs[6];

    Square orig_ktm_sq = pos.square<KING>(stm);

    // transformation to put king in left lower quadrant
    int8_t flip = ((orig_ktm_sq & RightHalfBB) ? 7 : 0) ^ ((orig_ktm_sq & TopHalfBB) ? 56 : 0);
    // transformation to put first piece that breaks diagonal symmetry below diagonal
    int8_t swap = 0;

    Square ktm_sq = maybe_update_swap_and_transform(orig_ktm_sq, flip, is_diag_symmetric, swap);

    Square orig_kntm_sq = pos.square<KING>(~stm);
    Square kntm_sq = maybe_update_swap_and_transform(orig_kntm_sq, flip, is_diag_symmetric, swap);
    
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


    int sqs_ixs[4];
    uint64_t n_available_squares = 62;
    int n_occupied_sqs = 2;

    for (Color c: {~stm, stm}) {
        for (PieceType p: {QUEEN, ROOK, BISHOP, KNIGHT}) {
            Bitboard pieceBB = pos.pieces(c, p);
            if (pieceBB) {
                Bitboard transformedBB = maybe_update_swap_and_transform_bb(pieceBB, flip, is_diag_symmetric, swap);
                int piece_count = 0;
                while (transformedBB) {
                    Square sq = pop_lsb(transformedBB);
                    int k = insert_count_lt_squares(occupied_sqs, n_occupied_sqs, sq);
                    sqs_ixs[piece_count] = (int) sq - k + piece_count;
                    piece_count++;
                }
                uint64_t tril_ix = tril_to_linear(piece_count, sqs_ixs);
                ix += tril_ix * multiplier;
                multiplier *= number_of_ordered_tuples(n_available_squares, piece_count);
                n_available_squares -= piece_count;                
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
    int n_occupied_sqs = 0;
    for (Color c: {~stm, stm}) {
        Bitboard pieceBB = pos.pieces(c, PAWN);
        assert(!more_than_one(pieceBB));
        if (pieceBB) {
            Square sq = transform(lsb(pieceBB), flip, 0);
            uint64_t k = insert_count_lt_squares(occupied_sqs, n_occupied_sqs, sq);

            if (first_pawn) {
                first_pawn = false;
                uint64_t sq_ix = (uint64_t) sq - k - 8;
                sq_ix = sq_ix - (sq_ix >> 3) * 4; // map to 0,1,...,23
                ix += sq_ix * multiplier;
                multiplier *= 24;
            } else {
                ix += ((uint64_t) sq - k - 8) * multiplier;
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

                uint64_t k = insert_count_lt_squares(occupied_sqs, n_occupied_sqs, sq);

                ix += ((uint64_t) sq - k) * multiplier;
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

    bool is_diag_symmetric = true;
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
        is_diag_symmetric = false; // disable diagonal symmetries
    }

    Square ktm_sq = maybe_update_swap_and_transform(orig_ktm_sq, flip, is_diag_symmetric, swap);
    pos2.put_piece(make_piece(stm, KING), ktm_sq);

    Square orig_kntm_sq = pos.square<KING>(~stm);
    Square ktnm_sq = maybe_update_swap_and_transform(orig_kntm_sq, flip, is_diag_symmetric, swap);
    pos2.put_piece(make_piece(~stm, KING), ktnm_sq);

    for (Color c: {~stm, stm}) {
        for (PieceType pt: {QUEEN, ROOK, BISHOP, KNIGHT, PAWN}) {
            Bitboard bb = pos.pieces(c, pt);
            if (bb) {
                Bitboard transformedBB = maybe_update_swap_and_transform_bb(bb, flip, is_diag_symmetric, swap);
                while (transformedBB) {
                    pos2.put_piece(make_piece(c,pt), pop_lsb(transformedBB));
                }
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