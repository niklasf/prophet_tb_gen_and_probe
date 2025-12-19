#include <iostream>
#include <vector>
#include "bitboard.h"
#include "values.h"
#include "kkx.h"
#include "triangular_indexes.h"
#include "egtb_iterator.h"
#include "eg_movegen.h"
// #include "linearize.h"
// #include "eg_position.h"
// #include "eg_movegen.h"
// #include "gen_egtb.h"
#include <fstream>
#include <unordered_map>

std::string egtb_id_from_pos( EGPosition pos) {
    std::ostringstream os;
    for (Color c : {pos.side_to_move(), ~pos.side_to_move()}) {
        for (PieceType pt : {KING, QUEEN, ROOK, BISHOP, KNIGHT, PAWN}) {
            for (int i = 0; i < pos.count(c, pt); i++) {
                os << PieceToChar[pt];
            }
        }
    }
    return os.str();
}

std::vector<Move> get_pv(std::unordered_map<std::string, EGTB*> id_to_egtb, EGPosition pos) {
    int16_t val = id_to_egtb[egtb_id_from_pos(pos)]->query_postion(pos);
    std::vector<Move> pv;
    while (val != LOSS_IN(0)) {
        // std::cout << val << std::endl;
        val = -val;
        if (val > 0) val++;
        if (val < 0) val--;
        bool found = false;
        for (Move move : EGMoveList<FORWARD>(pos)) {
            UndoInfo u = pos.do_move(move);
            int16_t move_val = id_to_egtb[egtb_id_from_pos(pos)]->query_postion(pos);
            // std::cout << move_to_uci(move) << " " << move_val << std::endl;
            pos.undo_move(move, u);
            if (move_val == val) {
                pv.push_back(move);
                pos.do_move(move);
                // std::cout << move_to_uci(move) << std::endl;
                found = true;
                break;
            }
        }
        assert (found);
    }
    return pv;
}

int main() {
    Bitboards::init();
    init_kkx_table();
    init_kkp_table();
    init_tril();

    std::string folder = "egtbs";
    std::vector<std::string> egtb_ids = get_egtb_identifiers(0, 4);
    std::cout << "EGTB count: " << egtb_ids.size() << std::endl;

    uint64_t total_position_count = 0;
    int count = 0;
    
    std::ofstream file("mates.txt");


    std::unordered_map<std::string, EGPosition> longest_mates = {};
    std::unordered_map<std::string, EGTB*> id_to_egtb = {};

    for (std::string egtb_id : egtb_ids) {
        count++;
        EGTB* egtb = new EGTB(egtb_id);
        std::cout << count << ". " << egtb_id  << ": " << egtb->num_pos << std::endl;
        total_position_count += egtb->num_pos;

        if (!egtb->exists(folder)) {
            std::cout << "NOT FOUND." << std::endl;
            continue;
        }
        id_to_egtb[egtb_id] = egtb;
        // egtb.maybe_decompress_and_load_egtb(folder);
        egtb->init_compressed_tb(folder);

        int16_t longest_mate = WIN_IN(0) + 1;
        uint64_t longest_mate_ix = 0;
        for (uint64_t ix = 0; ix < egtb->num_pos; ix++) {
            int16_t val = egtb->get_value(ix);
            assert (IS_SET(val));
            if (0 < val && val < longest_mate) {
                longest_mate = val;
                longest_mate_ix = ix;
            }
        }
        if (longest_mate == WIN_IN(0) + 1) {
            std::cout << "no win." << std::endl;
            file << count << ". " << egtb_id  << ": " << "no win." << "\n";
        } else {
            EGPosition pos;
            pos.reset();
            egtb->pos_at_ix(pos, longest_mate_ix, WHITE);
            longest_mates[egtb_id] = pos;
            int mate_in_plies = WIN_IN(0) - egtb->get_value(longest_mate_ix);
            std::cout << longest_mate_ix << " " << pos.fen() << " mate in " << mate_in_plies << " plies (" << mate_in_plies / 2 + 1 << " moves) " << egtb_id_from_pos(pos) << std::endl;
            file << count << ". " << egtb_id  << ": " << pos.fen() << " mate in " << mate_in_plies << " plies (" << mate_in_plies / 2 + 1 << " moves)" << "\n";
        }

        // egtb.free_tb();

    }

    for (const auto& [egtb_id, pos] : longest_mates) {
        int mate_in_plies = WIN_IN(0) - id_to_egtb[egtb_id_from_pos(pos)]->query_postion(pos);
        std::cout << egtb_id << ": " << pos.fen() << " " << mate_in_plies;
        for (Move move : get_pv(id_to_egtb, pos)) {
            std::cout << " " << move_to_uci(move) ;
        }
        std::cout << std::endl;
    }

    std::cout << "total_position_count: " << total_position_count << std::endl;
}