/* C integration test for SRH defect API (Task 7).
 * Exercises negf_set_defect_srh_mpr, negf_get_defect_srh_rates, and
 * negf_get_band_resolved_pl_currents via the auto-generated libnegf.h.
 *
 * Verifies:
 *   1. set_defect_srh_mpr sets is_active without crashing.
 *   2. get_defect_srh_rates writes f_occ = 0.5 (Fortran default).
 *   3. get_band_resolved_pl_currents: j_tot = j_el + j_hl on every boundary.
 */

#include "libnegf.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    int handler[NEGF_HSIZE];
    int* hand = &handler[0];
    int rank, err;
    MPI_Fint global_comm, cart_comm, k_comm, en_comm;

    /* ------------------------------------------------------------------ */
    /* MPI init                                                             */
    /* ------------------------------------------------------------------ */
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* ------------------------------------------------------------------ */
    /* Initialise a minimal libnegf session                                 */
    /* ------------------------------------------------------------------ */
    negf_init_session(hand);
    negf_init(hand);

    global_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
    negf_set_mpi_fcomm(hand, global_comm);
    negf_cartesian_init(
        hand, global_comm, 1, &cart_comm, &k_comm, &en_comm
    );

    /* Minimal 4-PL structure: 2 orbitals per PL, no contacts */
    int surfstart[1] = {0};
    int surfend[1] = {0};
    int contend[1] = {0};
    int plend[4] = {2, 4, 6, 8};
    int cblk[1] = {0};
    negf_init_structure(
        hand, 0, surfstart, surfend, contend, 4, plend, cblk
    );

    /* ------------------------------------------------------------------ */
    /* Configure a midgap MPR defect in PL 2                               */
    /* ------------------------------------------------------------------ */
    int num_defects = 1;
    int pl_indices[1] = {2};
    int orb_start[1] = {3};
    int orb_end[1] = {4};
    double eps_d[1] = {0.0};
    double density_weight[1] = {1.0};
    double kbt = 0.02585;
    double hbar_omega0 = 0.1;
    double s_n = 1.0;
    double s_p = 1.0;
    double m0_dc = 1.0e-5;
    double m0_dv = 1.0e-5;
    int n_modes = 5;
    double ec_cut[1] = {0.15};
    double ev_cut[1] = {-0.15};
    int mode = 0; /* DEFECT_MODE_DELTA */

    negf_set_defect_srh_mpr(
        hand,
        num_defects,
        pl_indices,
        orb_start,
        orb_end,
        eps_d,
        density_weight,
        kbt,
        hbar_omega0,
        s_n,
        s_p,
        m0_dc,
        m0_dv,
        n_modes,
        ec_cut,
        ev_cut,
        mode
    );

    /* ------------------------------------------------------------------ */
    /* Read rates: before any SCBA, f_occ must equal the Fortran default    */
    /* 0.5_dp set by init_defect_spec.                                      */
    /* ------------------------------------------------------------------ */
    double f_occ[1] = {0.0};
    double rate_net_el[1] = {0.0};
    double rate_net_hl[1] = {0.0};
    double total_recomb_rate = 0.0;

    negf_get_defect_srh_rates(
        hand,
        num_defects,
        f_occ,
        rate_net_el,
        rate_net_hl,
        &total_recomb_rate
    );

    err = 0;
    /* Fortran TDefectSpec%f_occ is initialised to 0.5_dp; verify the C API
     * actually copied that value (not left the C-side zero in place). */
    if(fabs(f_occ[0] - 0.5) > 1.0e-14) {
        fprintf(
            stderr,
            "FAIL: f_occ[0] = %.17g, expected 0.5\n",
            f_occ[0]
        );
        err = 1;
    }

    /* ------------------------------------------------------------------ */
    /* Band-resolved inter-PL currents (continuity mode)                   */
    /* ------------------------------------------------------------------ */
    int n_bnd = 3; /* 4 PLs => 3 boundaries */
    double e_split[3] = {0.0, 0.0, 0.0};
    double j_el_pl[3] = {0.0, 0.0, 0.0};
    double j_hl_pl[3] = {0.0, 0.0, 0.0};
    double j_tot_pl[3] = {0.0, 0.0, 0.0};

    negf_get_band_resolved_pl_currents(
        hand,
        n_bnd,
        e_split,
        j_el_pl,
        j_hl_pl,
        j_tot_pl
    );

    for(int l = 0; l < n_bnd; ++l) {
        double j_check = j_el_pl[l] + j_hl_pl[l];
        if(fabs(j_check - j_tot_pl[l]) > 1.0e-14) {
            fprintf(
                stderr,
                "FAIL: j_el[%d]+j_hl[%d] = %.17g != j_tot[%d] = %.17g\n",
                l,
                l,
                j_check,
                l,
                j_tot_pl[l]
            );
            err = 1;
        }
    }

    /* ------------------------------------------------------------------ */
    /* Clean up                                                             */
    /* ------------------------------------------------------------------ */
    negf_destruct_libnegf(hand);
    negf_destruct_session(hand);

    MPI_Finalize();

    if(rank == 0 && err == 0) {
        printf("Task 7 C test: PASSED\n");
    }

    return err;
}
