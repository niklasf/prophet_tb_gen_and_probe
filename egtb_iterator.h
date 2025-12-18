#ifndef EGTB_ITERATOR_H_INCLUDED
#define EGTB_ITERATOR_H_INCLUDED

#include <vector>
#include <unordered_set>
#include "egtb.h"

void place_piece(Piece p, int* pieces1, int* pieces2) {
    int* pieces = color_of(p) == WHITE ? pieces1 : pieces2;
    pieces[type_of(p)]++;
}
void unplace_piece(Piece p, int* pieces1, int* pieces2) {
    int*pieces = color_of(p) == WHITE ? pieces1 : pieces2;
    pieces[type_of(p)]--;
}

std::vector<std::string> get_egtb_identifiers(int MIN_PIECE_COUNT, int MAX_PIECE_COUNT, int MIN_PAWN_COUNT = 0, int MAX_PAWN_COUNT = 4) {
    
    std::vector<int> pieces1(6);
    std::vector<int> pieces2(6);

    pieces1 = {0, 0, 0, 0, 0, 0};
    pieces2 = {0, 0, 0, 0, 0, 0};

    std::unordered_set<std::string> egtbs_set = {};
    std::vector<std::string> egtbs = {};
    Piece PIECES_ARR[] = {NO_PIECE, W_PAWN, W_KNIGHT, W_BISHOP, W_ROOK, W_QUEEN, B_PAWN, B_KNIGHT, B_BISHOP, B_ROOK, B_QUEEN};
    for (int piece_count = MIN_PIECE_COUNT; piece_count <= MAX_PIECE_COUNT; piece_count++) {
        for (int pawn_count = MIN_PAWN_COUNT; pawn_count <= std::min(piece_count,MAX_PAWN_COUNT); pawn_count++ ) {
            for (Piece p1 : PIECES_ARR) {
                for (Piece p2 : PIECES_ARR) {
                    for (Piece p3 : PIECES_ARR) {
                        for (Piece p4 : PIECES_ARR) {
                            if ((p1 != NO_PIECE) + (p2 != NO_PIECE) + (p3 != NO_PIECE) + (p4 != NO_PIECE) != piece_count) continue;
                            if ((piece_count == 0) && (p1 != NO_PIECE) + (p2 != NO_PIECE) + (p3 != NO_PIECE) + (p4 != NO_PIECE) > 0) continue;
                            if ((piece_count == 1) && (p2 != NO_PIECE) + (p3 != NO_PIECE) + (p4 != NO_PIECE) > 0) continue;
                            if ((piece_count == 2) && (p3 != NO_PIECE) + (p4 != NO_PIECE) > 0) continue;
                            if ((piece_count == 3) && (p4 != NO_PIECE) > 0) continue;

                            if ((type_of(p1) == PAWN) + (type_of(p2) == PAWN) + (type_of(p3) == PAWN) + (type_of(p4) == PAWN) != pawn_count) {
                                continue;
                            }

                            if (p1 != NO_PIECE) place_piece(p1, &pieces1[0], &pieces2[0]);
                            if (p2 != NO_PIECE) place_piece(p2, &pieces1[0], &pieces2[0]);
                            if (p3 != NO_PIECE) place_piece(p3, &pieces1[0], &pieces2[0]);
                            if (p4 != NO_PIECE) place_piece(p4, &pieces1[0], &pieces2[0]);

                            std::string id = get_egtb_identifier(&pieces1[0], &pieces2[0]);
                            if (egtbs_set.insert(id).second) {
                                egtbs.push_back(id);
                            }

                            if (p1 != NO_PIECE) unplace_piece(p1, &pieces1[0], &pieces2[0]);
                            if (p2 != NO_PIECE) unplace_piece(p2, &pieces1[0], &pieces2[0]);
                            if (p3 != NO_PIECE) unplace_piece(p3, &pieces1[0], &pieces2[0]);
                            if (p4 != NO_PIECE) unplace_piece(p4, &pieces1[0], &pieces2[0]);
                        }
                    }
                }
            }
        }
    }

    return egtbs;
}
#endif