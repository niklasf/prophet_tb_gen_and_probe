#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <prophet.h>
#include "misc.h"
#ifdef OMP
#include <omp.h>
#endif

uint64_t probe_count = 0;
int16_t probe_position_dctx_counted(EGPosition pos, DecompressCtx* dctx) {
    #pragma omp atomic
    probe_count++;
    return probe_position_dctx(pos, dctx);
}

std::vector<Move> get_mate_line(EGPosition pos, DecompressCtx* dctx) {
    int16_t val = probe_position_dctx_counted(pos, dctx);
    std::vector<Move> pv;
    while (val != LOSS_IN(0)) {
        val = -val;
        if (val > 0) val++;
        if (val < 0) val--;
        bool found = false;
        for (Move move : EGMoveList<FORWARD>(pos)) {
            UndoInfo u = pos.do_move(move);
            int16_t move_val = probe_position_dctx_counted(pos, dctx);
            pos.undo_move(move, u);
            if (move_val == val) {
                pv.push_back(move);
                pos.do_move(move);
                found = true;
                break;
            }
        }
        assert (found);
    }
    return pv;
}

struct CSVEntry {
    std::string egtb_id;
    std::string fen;
    std::string ply_str;
};

#define LOAD_HALF 1

int main(int argc, char *argv[]) {
    assert (argc > 0);
    int nthreads = atoi(argv[1]);
    #ifndef OMP
    assert(nthreads == 1);
    #endif
    std::cout << "Probing with " << nthreads << " threads." << std::endl;


    std::string folder = "egtbs";
    init_prophet_tb(folder);
    

#if LOAD_HALF
    for (std::string egtb_id : get_egtb_identifiers(4, 4)) {
        std::string mirror_egtb_id = get_mirror_id(egtb_id);
        if (mirror_egtb_id != egtb_id) {
            if (id_to_egtb[egtb_id] == nullptr || id_to_egtb[mirror_egtb_id] == nullptr) continue;
            if (id_to_egtb[egtb_id]->CTB->compressed_filesize < id_to_egtb[mirror_egtb_id]->CTB->compressed_filesize) {
                delete id_to_egtb[mirror_egtb_id];
                id_to_egtb[mirror_egtb_id] = nullptr;
            } else {
                delete id_to_egtb[egtb_id];
                id_to_egtb[egtb_id] = nullptr;
            }
        }
    }
#endif

    
    uint64_t count = 0;
    uint64_t compressed_filesize = 0;
    for (std::string egtb_id : get_egtb_identifiers()) {
        if (id_to_egtb[egtb_id] != nullptr) {
            count++;
            compressed_filesize += id_to_egtb[egtb_id]->CTB->compressed_filesize;
        }
    }
    std::cout << "Loaded " << count << " EGTBs (" << (int) ceil((double) compressed_filesize / (1024*1024*1024)) << "GiB)" << std::endl;

    std::ifstream file("longest_mates.csv");std::string line;

    std::vector<CSVEntry> entries;
    std::getline(file, line); // header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string egtb_id, fen, ply_str;

        std::getline(ss, egtb_id, ',');
        std::getline(ss, fen, ',');
        std::getline(ss, ply_str, ',');

        entries.push_back({egtb_id, fen, ply_str});
    }

    TimePoint t0 = now();
    #pragma omp parallel num_threads(nthreads)
    {
        DecompressCtx* dctx = new DecompressCtx();

        #pragma omp for schedule(static)
        for (uint i = 0; i < entries.size(); i++) {
            CSVEntry entry = entries[i];
            // std::cout << i << ". ID: " << entry.egtb_id << " FEN: " << entry.fen << " Plies: " << entry.ply_str << std::endl;
            if (entry.ply_str != "") {
                EGPosition pos;
                pos.from_fen(entry.fen);
                get_mate_line(pos, dctx);
                // for (Move move : get_pv(pos)) {
                //     std::cout << " " << move_to_uci(move) ;
                // }
                // std::cout << std::endl;
            }
        }

        #pragma omp critical
        {
            delete dctx;
        }
    }

    TimePoint t1 = now();
    std::cout << "Finished in " << (double) (t1-t0) / 1000 << "s. Probe count = " << probe_count << " (" << (double) (t1-t0) / probe_count << "ms/probe)" << std::endl;

}
