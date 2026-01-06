#include "prophet.h"
#include "bitboard.h"
#include "kkx.h"
#include "triangular_indexes.h"
#include "egtb.h"
#include "egtb_ids.h"
#include "eg_movegen.h"

std::unordered_map<std::string, EGTB*> id_to_egtb = {};

std::string egtb_id_from_pos(const EGPosition pos) {
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

std::string get_mirror_id(const std::string& egtb_id) {
    size_t firstK = egtb_id.find('K');
    size_t secondK = egtb_id.find('K', firstK + 1);
    std::string stuff1 = egtb_id.substr(firstK + 1, secondK - firstK - 1);
    std::string stuff2 = egtb_id.substr(secondK + 1);
    return "K" + stuff2 + "K" + stuff1;
}

int16_t query_position(EGPosition pos, DecompressCtx* dctx) {
    std::string egtb_id = egtb_id_from_pos(pos);

    if (id_to_egtb[egtb_id] != nullptr) {
        return id_to_egtb[egtb_id]->query_postion_dctx(pos, dctx);
    } else {
        EGMoveList movelist = EGMoveList<FORWARD>(pos);
        if (movelist.size() == 0) {
            if (pos.stm_in_check()) {
                return LOSS_IN(0);
            } else {
                return 0;
            }
        }
        int16_t max_val = LOSS_IN(0);
        for (Move move : movelist) {
            UndoInfo u = pos.do_move(move);
            int16_t val = query_position(pos, dctx);
            max_val = std::max(max_val, (int16_t) -val);
            pos.undo_move(move, u);
        }
        if (max_val > 0) max_val--;
        if (max_val < 0) max_val++;
        return max_val;
    }
}


void init_prophet_tb(std::string path) {
    Bitboards::init();
    init_kkx_table();
    init_kkp_table();
    init_tril();
    init_egtb_id_to_ix();

    int count = 0;
    for (std::string egtb_id : get_egtb_identifiers(0, 4)) {
        EGTB* egtb = new EGTB(egtb_id);
        if (egtb->exists(path)) {
            egtb->init_compressed_tb(path);
            id_to_egtb[egtb_id] = egtb;
            count++;
        } else {
            delete egtb;
        }
    }
    std::cout << "Init ProphetTB: Loaded " << count << " tables." << std::endl;
}

void deinit_prophet_tb() {
    for (std::string egtb_id : get_egtb_identifiers(0, 4)) {
        delete id_to_egtb[egtb_id];
        id_to_egtb[egtb_id] = nullptr;
    }
}


DecompressCtx* CreateDecompressCtx() {
    return new DecompressCtx(32768);
}
int16_t probe_dctx(Piece pieces[6], Square squares[6], DecompressCtx* dctx) {
    EGPosition pos;
    pos.reset();
    for (int i = 0; i < 6; i++) {
        if (pieces[i] != NO_PIECE) {
            pos.put_piece(pieces[i], squares[i]);
        }
    }
    return query_position(pos, dctx);
}


int16_t probe(Piece pieces[6], Square squares[6]) {
    DecompressCtx* dctx = new DecompressCtx(32768);
    int16_t val = probe_dctx(pieces, squares, dctx);
    delete dctx;
    return val;
}