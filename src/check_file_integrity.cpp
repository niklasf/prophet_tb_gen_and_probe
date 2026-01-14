#include <iostream>
#include <vector>
#include <prophet.h>
#include "misc.h"
#ifdef OMP
#include <omp.h>
#endif


int main(int argc, char *argv[]) {
    assert (argc == 3);
    std::string folder = argv[1];
    int nthreads = atoi(argv[2]);
    #ifndef OMP
    assert(nthreads == 1);
    #endif

    init_prophet_tb(folder);

    TimePoint t0 = now();
    std::vector<std::string> egtb_ids = get_egtb_identifiers(0, 4);

    bool all_ok = true;
    int count = 0;
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic) reduction(&&:all_ok)
    for (uint i = 0; i < egtb_ids.size(); i++) {
        std::string egtb_id = egtb_ids[i];
        EGTB egtb = EGTB(egtb_id);

        if (!egtb.exists(folder)) {
            #pragma omp critical
            {
                count++;
                std::cout << "Finished check " << count << "/" << egtb_ids.size() << ": " << egtb_id << " does not exist." << std::endl;
            }
            continue;
        }
        egtb.init_compressed_tb(folder);
        bool ok = egtb.CTB->check_integrity(1);
        all_ok = all_ok && ok;

        #pragma omp critical
        {
            count++;
            std::cout << "Finished check " << count << "/" << egtb_ids.size() << ": " << egtb_id << (ok ? " OK" : " ERROR") << std::endl;
        }

    }
    
    TimePoint t1 = now();
    std::cout << "Finished in " << (t1-t0) / 1000 << " seconds." << std::endl;
    if (all_ok)
        std::cout << "All ok!" << std::endl;
    else
        std::cout << "There were some errors :(" << std::endl;

}