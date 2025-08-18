
#ifndef EG_MOVEGEN_H_INCLUDED
#define EG_MOVEGEN_H_INCLUDED

#include "eg_position.h"

enum EGGenType {
    FORWARD, // all legal moves except castling and en-passant (for now)
    FWD_EVASIONS, // if in check all legal evasions
    FWD_NON_EVASIONS,
    REVERSE // all non-capturing backward moves (no castling, no en-passant)
};


template<Direction offset>
inline Move* splat_pawn_moves(Move* moveList, Bitboard to_bb) {
    while (to_bb)
    {
        Square to   = pop_lsb(to_bb);
        *moveList++ = Move(to - offset, to);
    }
    return moveList;
}

inline Move* splat_moves(Move* moveList, Square from, Bitboard to_bb) {
    while (to_bb)
        *moveList++ = Move(from, pop_lsb(to_bb));
    return moveList;
}

template<EGGenType Type, Direction D>
Move* make_promotions(Move* moveList, Square to) {
    *moveList++ = Move::make<PROMOTION>(to - D, to, QUEEN);
    *moveList++ = Move::make<PROMOTION>(to - D, to, ROOK);
    *moveList++ = Move::make<PROMOTION>(to - D, to, BISHOP);
    *moveList++ = Move::make<PROMOTION>(to - D, to, KNIGHT);
    return moveList;
}

template<Color Us,EGGenType Type>
Move* generate_pawn_moves(const EGPosition& pos, Move* moveList, Bitboard target, Bitboard checkers) {

    constexpr Color     Them     = ~Us;
    constexpr Bitboard  TRank7BB = (Us == WHITE ? Rank7BB : Rank2BB);
    constexpr Bitboard  TRank3BB = (Us == WHITE ? Rank3BB : Rank6BB);
    constexpr Direction Up       = pawn_push(Us);
    constexpr Direction UpRight  = (Us == WHITE ? NORTH_EAST : SOUTH_WEST);
    constexpr Direction UpLeft   = (Us == WHITE ? NORTH_WEST : SOUTH_EAST);

    const Bitboard emptySquares = ~pos.pieces();
    const Bitboard enemies      = Type == FWD_EVASIONS ? checkers : pos.pieces(Them);

    Bitboard pawnsOn7    = pos.pieces(Us, PAWN) & TRank7BB;
    Bitboard pawnsNotOn7 = pos.pieces(Us, PAWN) & ~TRank7BB;

    // Single and double pawn pushes, no promotions
    {
        Bitboard b1 = shift<Up>(pawnsNotOn7) & emptySquares;
        Bitboard b2 = shift<Up>(b1 & TRank3BB) & emptySquares;

        if constexpr (Type == FWD_EVASIONS)  // Consider only blocking squares
        {
            b1 &= target;
            b2 &= target;
        }

        moveList = splat_pawn_moves<Up>(moveList, b1);
        moveList = splat_pawn_moves<Up + Up>(moveList, b2);
    }

    // Promotions and underpromotions
    if (pawnsOn7)
    {
        Bitboard b1 = shift<UpRight>(pawnsOn7) & enemies;
        Bitboard b2 = shift<UpLeft>(pawnsOn7) & enemies;
        Bitboard b3 = shift<Up>(pawnsOn7) & emptySquares;

        if constexpr (Type == FWD_EVASIONS)
            b3 &= target;

        while (b1)
            moveList = make_promotions<Type, UpRight>(moveList, pop_lsb(b1));

        while (b2)
            moveList = make_promotions<Type, UpLeft>(moveList, pop_lsb(b2));

        while (b3)
            moveList = make_promotions<Type, Up>(moveList, pop_lsb(b3));
    }
   
    {
        Bitboard b1 = shift<UpRight>(pawnsNotOn7) & enemies;
        Bitboard b2 = shift<UpLeft>(pawnsNotOn7) & enemies;

        moveList = splat_pawn_moves<UpRight>(moveList, b1);
        moveList = splat_pawn_moves<UpLeft>(moveList, b2);

        // TODO: en-passant
    }

    return moveList;
}


template<Color Us, PieceType Pt>
Move* generate_moves(const EGPosition& pos, Move* moveList, Bitboard target) {

    static_assert(Pt != KING && Pt != PAWN, "Unsupported piece type in generate_moves()");

    Bitboard bb = pos.pieces(Us, Pt);

    while (bb)
    {
        Square   from = pop_lsb(bb);
        Bitboard b    = attacks_bb<Pt>(from, pos.pieces()) & target;

        moveList = splat_moves(moveList, from, b);
    }

    return moveList;
}

template<Color Us, EGGenType Type>
Move* generate_all(const EGPosition& pos, Move* moveList, Bitboard checkers) {
    const Square ksq = pos.square<KING>(Us);
    Bitboard     target;

    // Skip generating non-king moves when in double check
    if (Type != FWD_EVASIONS || !more_than_one(checkers)) {
        target = Type == FWD_EVASIONS ? between_bb(ksq, lsb(checkers)): pos.pieces(Us);

        moveList = generate_pawn_moves<Us, Type>(pos, moveList, target, checkers);
        moveList = generate_moves<Us, KNIGHT>(pos, moveList, target);
        moveList = generate_moves<Us, BISHOP>(pos, moveList, target);
        moveList = generate_moves<Us, ROOK>(pos, moveList, target);
        moveList = generate_moves<Us, QUEEN>(pos, moveList, target);
    }

    Bitboard b = attacks_bb<KING>(ksq) & (Type == FWD_EVASIONS ? ~pos.pieces(Us) : target);

    moveList = splat_moves(moveList, ksq, b);

    return moveList;
}

template<EGGenType>
Move* generate(const EGPosition& pos, Move* moveList);

template<>
Move* generate<FORWARD>(const EGPosition& pos, Move* moveList) {

    Color    us     = pos.side_to_move();
    Bitboard pinned = pos.blockers_for_king(us) & pos.pieces(us);
    Square   ksq    = pos.square<KING>(us);
    Move*    cur    = moveList;
    Bitboard checkers = pos.checkers(us);

    if (us == WHITE)
        moveList = checkers ? generate_all<WHITE,FWD_EVASIONS>(pos, moveList, checkers) : generate_all<WHITE,FWD_NON_EVASIONS>(pos, moveList, checkers);
    else
        moveList = checkers ? generate_all<BLACK,FWD_EVASIONS>(pos, moveList, checkers) : generate_all<BLACK,FWD_NON_EVASIONS>(pos, moveList,  checkers);
    
    while (cur != moveList) {
        if (pinned & cur->from_sq() && !(line_bb(cur->from_sq(), cur->to_sq()) & ksq)) {
            *cur = *(--moveList);
        } else if (cur->from_sq() == ksq && pos.attackers_to_exist(cur->to_sq(), pos.pieces() ^ ksq, ~us)) {
            *cur = *(--moveList);
        } else {
            ++cur;
        }
    }
    return moveList;
}


template<EGGenType T>
struct EGMoveList {
    explicit EGMoveList(const EGPosition& pos) :
        last(generate<T>(pos, moveList)) {}
    const Move* begin() const { return moveList; }
    const Move* end() const { return last; }
    size_t      size() const { return last - moveList; }

private:
    Move moveList[MAX_MOVES], *last;
};

#endif