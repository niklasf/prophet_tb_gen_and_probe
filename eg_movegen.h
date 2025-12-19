
#ifndef EG_MOVEGEN_H_INCLUDED
#define EG_MOVEGEN_H_INCLUDED

#include "eg_position.h"

enum EGGenType {
    FORWARD, // all legal moves except castling and en-passant (for now)
    FWD_EVASIONS, // if in check all legal evasions
    FWD_NON_EVASIONS,
    REVERSE, // all non-capturing backward moves (no castling, no en-passant)
};


template<EGGenType>
Move* generate(const EGPosition& pos, Move* moveList, PieceType captured, PieceType promotion);

template<EGGenType T>
struct EGMoveList {
    explicit EGMoveList(const EGPosition& pos, PieceType captured = NO_PIECE_TYPE, PieceType promotion = NO_PIECE_TYPE) :
        last(generate<T>(pos, moveList, captured, promotion)) {}
    const Move* begin() const { return moveList; }
    const Move* end() const { return last; }
    size_t      size() const { return last - moveList; }

private:
    Move moveList[MAX_MOVES], *last;
};

#endif