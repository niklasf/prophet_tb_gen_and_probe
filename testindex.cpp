

void count_broken(EGTB* egtb) {

    int piece_count = 2;
    for (int i = 0; i < 6; i++) piece_count += egtb->stm_pieces[i] + egtb->sntm_pieces[i];

    uint64_t n_sntm_in_check = 0;
    uint64_t n_piece_collision = 0;
    uint64_t n_illegal_ep = 0;
    uint64_t n_unused = 0;
    uint64_t highest_used_ix = 0;

    #pragma omp parallel for num_threads(32) schedule(static,64) reduction(+:n_sntm_in_check) reduction(+:n_piece_collision) reduction(+:n_illegal_ep) reduction(+:n_unused) reduction(max:highest_used_ix)
    for (uint64_t ix = 0; ix < egtb->num_pos; ix++) {
        EGPosition pos;
        pos.reset();
        egtb->pos_at_ix(pos, ix, WHITE);
        if (popcount(pos.pieces()) != piece_count) {
            n_piece_collision++;
        } else if (pos.sntm_in_check()) {
            n_sntm_in_check++;
        } else if ((ix >= egtb->num_nonep_pos) && !pos.check_ep(pos.ep_square())) {
            n_illegal_ep++;
        } else if (egtb->ix_from_pos(pos) != ix) {
            n_unused++;
        } else {
            highest_used_ix = ix;
        }
    }

    printf("n_sntm_in_check=%ld (%.2f), n_piece_collision=%ld (%.2f), n_illegal_ep=%ld (%.2f), n_unused=%ld (%.2f), highest_used_ix=%ld\n",
        n_sntm_in_check, n_sntm_in_check * 100.0 / egtb->num_pos,
        n_piece_collision, n_piece_collision * 100.0 / egtb->num_pos,
        n_illegal_ep, n_illegal_ep * 100.0 / egtb->num_pos,
        n_unused, n_unused * 100.0 / egtb->num_pos,
        highest_used_ix
    );
}

void test_index() {

    //                0, P, N, B, R, Q
    int wpieces[6] = {0, 2, 0, 0, 0, 0};
    int bpieces[6] = {0, 1, 0, 0, 0, 0};
    bool check_ep = false;

    if (check_ep) assert(wpieces[PAWN] > 0 && bpieces[PAWN] > 0);

    EGTB wtm_table = EGTB(wpieces, bpieces);
    std::cout << wtm_table.id << ": nonep: " << wtm_table.num_nonep_pos << ", ep: " << wtm_table.num_ep_pos << std::endl;
    // count_broken(&wtm_table);

    EGTB btm_table = EGTB(bpieces, wpieces);
    std::cout << btm_table.id << ": nonep: " << btm_table.num_nonep_pos << ", ep: " << btm_table.num_ep_pos << std::endl;
    // count_broken(&btm_table);


    uint64_t count = 0;
    for (Color stm: {WHITE, BLACK}) {

        EGTB egtb = stm == WHITE ? wtm_table : btm_table;


        PieceType pts[4] = {NO_PIECE_TYPE, NO_PIECE_TYPE, NO_PIECE_TYPE, NO_PIECE_TYPE};
        Color cs[4] = {COLOR_NB, COLOR_NB, COLOR_NB, COLOR_NB};

        int i = 0;
        for (Color c: {~stm, stm}) {
            int* c_pieces = (c == WHITE) ? wpieces : bpieces;
            for (PieceType pt = QUEEN; pt >= PAWN; --pt) {
                for (int j = 0; j < c_pieces[pt]; ++j) {
                    if (i >= 4) { std::cout << "Too many pieces!\n"; assert(false); }
                    pts[i] = pt;
                    cs[i] = c;
                    i++;
                }
            }
        }

        EGPosition pos1;
        EGPosition pos2;
        EGPosition pos3;
        EGPosition pos4;


        for (Square k1_sq = SQ_A1; k1_sq <= SQ_H8; ++k1_sq) {
            for (Square k2_sq = SQ_A1; k2_sq <= SQ_H8; ++k2_sq) {
                if (k2_sq == k1_sq) { continue; }
                if ((PseudoAttacks[KING][k1_sq] & k2_sq) != 0) { continue; }
                for (Square p1_sq = (pts[0] != NO_PIECE_TYPE) ? SQ_A1 : SQ_NONE; p1_sq <= ((pts[0] != NO_PIECE_TYPE) ? SQ_H8 : SQ_NONE); ++p1_sq) {
                    if (pts[0] == PAWN && !(p1_sq & PawnSquaresBB)) { continue; }
                    if (p1_sq != SQ_NONE && (p1_sq == k1_sq || p1_sq == k2_sq)) { continue; }
                    for (Square p2_sq = (pts[1] != NO_PIECE_TYPE) ? SQ_A1 : SQ_NONE; p2_sq <= ((pts[1] != NO_PIECE_TYPE) ? SQ_H8 : SQ_NONE); ++p2_sq) {
                        if (pts[1] == PAWN && !(p2_sq & PawnSquaresBB)) { continue; }
                        if (p2_sq != SQ_NONE && (p2_sq == k1_sq || p2_sq == k2_sq || p2_sq == p1_sq)) { continue; }
                        for (Square p3_sq = (pts[2] != NO_PIECE_TYPE) ? SQ_A1 : SQ_NONE; p3_sq <= ((pts[2] != NO_PIECE_TYPE) ? SQ_H8 : SQ_NONE); ++p3_sq) {
                            if (pts[2] == PAWN && !(p3_sq & PawnSquaresBB)) { continue; }
                            if (p3_sq != SQ_NONE && (p3_sq == k1_sq || p3_sq == k2_sq || p3_sq == p1_sq || p3_sq == p2_sq)) { continue; }
                            for (Square p4_sq = (pts[3] != NO_PIECE_TYPE) ? SQ_A1 : SQ_NONE; p4_sq <= ((pts[3] != NO_PIECE_TYPE) ? SQ_H8 : SQ_NONE); ++p4_sq) {
                                if (pts[3] == PAWN && !(p4_sq & PawnSquaresBB)) { continue; }
                                if (p4_sq != SQ_NONE && (p4_sq == k1_sq || p4_sq == k2_sq || p4_sq == p1_sq || p4_sq == p2_sq || p4_sq == p3_sq)) { continue; }
                                Square ep_sq_1 = (!check_ep ? SQ_NONE : (stm == WHITE) ? SQ_A6 : SQ_A3);
                                Square ep_sq_2 = (!check_ep ? SQ_NONE : (stm == WHITE) ? SQ_H6 : SQ_H3);
                                for (Square ep_sq =ep_sq_1; ep_sq <= ep_sq_2; ++ep_sq) {

                                    pos1.reset();
                                    pos2.reset();
                                    pos3.reset();
                                    pos4.reset();

                                    pos1.set_side_to_move(stm);
                                    pos1.put_piece(make_piece(~stm, KING), k1_sq);
                                    pos1.put_piece(make_piece(stm,KING), k2_sq);
                                    if (p1_sq != SQ_NONE) pos1.put_piece(make_piece(cs[0],pts[0]), p1_sq);
                                    if (p2_sq != SQ_NONE) pos1.put_piece(make_piece(cs[1],pts[1]), p2_sq);
                                    if (p3_sq != SQ_NONE) pos1.put_piece(make_piece(cs[2],pts[2]), p3_sq);
                                    if (p4_sq != SQ_NONE) pos1.put_piece(make_piece(cs[3],pts[3]), p4_sq);
                                    if (ep_sq != SQ_NONE) {
                                        if (!pos1.check_ep(ep_sq)) { continue;}
                                        pos1.set_ep_square(ep_sq);
                                    }

                                    if (pos1.sntm_in_check()) { continue; }

                                    // std::cout << pos1;

                                    uint64_t ix = egtb.ix_from_pos(pos1);

                                    // bool verbose = true;
                                    bool verbose = false;

                                    // bool debug = true;
                                    bool debug = false;
                                    // bool debug = verbose;


                                    if (verbose) std::cout << pos1;

                                    if (verbose) std::cout << "A ix from pos\n";
                                    if (verbose) std::cout << "B pos at ix " << ix << "\n";
                                    egtb.pos_at_ix(pos2, ix, stm);
                                    if (verbose) std::cout << "C transform to canonical\n";
                                    transform_to_canoncial(pos1, pos3);
                                    if (verbose) std::cout << "D ix from canonical\n";
                                    uint64_t ix2 = egtb.ix_from_pos(pos3);
                                    if (verbose) std::cout << "E transformed pos at ix " << ix2 << "\n";
                                    egtb.pos_at_ix(pos4, ix2, stm);
                                    if (verbose) std::cout << "F\n";

                                    if (!pos2.is_equal(pos3) || ix != ix2 || !pos3.is_equal(pos4) || debug) {
                                        std::cout << pos1 << std::endl;
                                        std::cout << "vs at ix " << ix << std::endl;
                                        std::cout << pos2 << std::endl;
                                        std::cout << "vs transformed " << std::endl;
                                        std::cout << pos3 << std::endl;
                                        std::cout << "vs at ix " << ix2 << std::endl;
                                        std::cout << pos4 << std::endl;
                                        exit(1);
                                    }
                                    count++;
                                    // std::cout << "\n";
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "Checked " << count << " positions" << std::endl;
        count = 0;
    }
    
}