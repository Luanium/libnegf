/* C integration test for SRH defect API (Tasks 7 and 8).
 *
 * Part A (Task 7 checks):
 *   1. set_defect_srh_mpr sets is_active without crashing.
 *   2. get_defect_srh_rates writes f_occ = 0.5 (Fortran default).
 *   3. get_band_resolved_pl_currents: j_tot = j_el + j_hl on every boundary.
 *
 * Part B (Task 8 end-to-end solve):
 *   4. Read HR.dat/HI.dat, set up 2-contact system, forward bias.
 *   5. Configure midgap MPR defect, call negf_solve_landauer.
 *   6. Verify 0 < f_occ < 1, |rate_el - rate_hl| < 1e-10,
 *      total_recomb_rate > 0.
 *   7. MPI broadcast from rank 0 and verify all ranks agree < 1e-11.
 */

#include "libnegf.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    int handler[NEGF_HSIZE];
    int *hand = &handler[0];
    int rank, err;
    MPI_Fint global_comm, cart_comm, k_comm, en_comm;

    /* ------------------------------------------------------------------ */
    /* MPI init                                                             */
    /* ------------------------------------------------------------------ */
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    err = 0;

    /* ================================================================== */
    /* Part A: Task 7 API checks (minimal structure, no H/S)              */
    /* ================================================================== */
    {
        /* Initialise a minimal libnegf session */
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

        /* Configure a midgap MPR defect in PL 2 */
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

        /* Before any SCBA, f_occ must equal the Fortran default 0.5_dp */
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

        if(fabs(f_occ[0] - 0.5) > 1.0e-14) {
            fprintf(
                stderr,
                "FAIL A1: f_occ[0] = %.17g, expected 0.5\n",
                f_occ[0]
            );
            err = 1;
        }

        /* Band-resolved inter-PL currents (continuity mode) */
        int n_bnd = 3; /* 4 PLs => 3 boundaries */
        double e_split[3] = {0.0, 0.0, 0.0};
        double j_el_pl[3] = {0.0, 0.0, 0.0};
        double j_hl_pl[3] = {0.0, 0.0, 0.0};
        double j_tot_pl[3] = {0.0, 0.0, 0.0};

        negf_get_band_resolved_pl_currents(
            hand, n_bnd, e_split, j_el_pl, j_hl_pl, j_tot_pl
        );

        for(int l = 0; l < n_bnd; ++l) {
            double j_check = j_el_pl[l] + j_hl_pl[l];
            if(fabs(j_check - j_tot_pl[l]) > 1.0e-14) {
                fprintf(
                    stderr,
                    "FAIL A2: j_el[%d]+j_hl[%d]=%.17g"
                    " != j_tot[%d]=%.17g\n",
                    l,
                    l,
                    j_check,
                    l,
                    j_tot_pl[l]
                );
                err = 1;
            }
        }

        negf_destruct_libnegf(hand);
        negf_destruct_session(hand);
    }

    /* ================================================================== */
    /* Part B: Task 8 end-to-end solve with HR.dat/HI.dat                 */
    /* ================================================================== */
    {
        struct lnparams params;
        char realmat[7] = "HR.dat";
        char imagmat[7] = "HI.dat";

        /* Same structure as c_int_elph_deph: 6 PLs, 2 contacts */
        int surfstart2[2] = {61, 81};
        int surfend2[2] = {60, 80};
        int contend2[2] = {80, 100};
        int plend2[6] = {10, 20, 30, 40, 50, 60};
        int cblk2[2] = {6, 1};

        negf_init_session(hand);
        negf_init(hand);

        global_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
        negf_set_mpi_fcomm(hand, global_comm);
        negf_cartesian_init(
            hand, global_comm, 1, &cart_comm, &k_comm, &en_comm
        );

        negf_read_hs(hand, &realmat[0], &imagmat[0], 0);
        negf_set_s_id(hand, 100, 1);
        negf_init_contacts(hand, 2);
        negf_init_structure(
            hand,
            2,
            &surfstart2[0],
            &surfend2[0],
            &contend2[0],
            6,
            &plend2[0],
            &cblk2[0]
        );

        /* Forward bias: mu[0]=-0.3, mu[1]=0.3 */
        negf_get_params(hand, &params);
        /* Wide energy window covering band states of HR.dat/HI.dat and
         * defect resonances eps_d +/- l * hbar_omega0 (l = 1..b_modes). */
        params.emin = -2.0;
        params.emax = 2.0;
        params.estep = 0.05;
        params.kbt_t[0] = 0.001;
        params.kbt_t[1] = 0.001;
        params.mu[0] = -0.3;
        params.mu[1] = 0.3;
        params.verbose = 10;
        negf_set_params(hand, &params);

        /* Midgap defect in PL 3 (central), orbitals 21-30 */
        int b_pl_indices[1] = {3};
        int b_orb_start[1] = {21};
        int b_orb_end[1] = {30};
        double b_eps_d[1] = {0.0};
        double b_dw[1] = {1.0};
        double b_kbt = 0.001;
        double b_omega0 = 0.05;
        double b_sn = 2.0;
        double b_sp = 2.0;
        double b_m0dc = 1.0e-3;
        double b_m0dv = 1.0e-3;
        int b_modes = 5;
        double b_ec_cut[1] = {0.1};
        double b_ev_cut[1] = {-0.1};

        negf_set_defect_srh_mpr(
            hand,
            1,
            b_pl_indices,
            b_orb_start,
            b_orb_end,
            b_eps_d,
            b_dw,
            b_kbt,
            b_omega0,
            b_sn,
            b_sp,
            b_m0dc,
            b_m0dv,
            b_modes,
            b_ec_cut,
            b_ev_cut,
            0
        );

        negf_solve_landauer(hand);

        /* Read out rates */
        double b_f_occ[1] = {0.0};
        double b_rate_el[1] = {0.0};
        double b_rate_hl[1] = {0.0};
        double b_total = 0.0;

        negf_get_defect_srh_rates(
            hand, 1, b_f_occ, b_rate_el, b_rate_hl, &b_total
        );

        /* Verify 0 < f_occ < 1 */
        if(rank == 0) {
            if(b_f_occ[0] <= 0.0 || b_f_occ[0] >= 1.0) {
                fprintf(
                    stderr,
                    "FAIL B1: f_occ[0]=%.17g not in (0,1)\n",
                    b_f_occ[0]
                );
                err = 1;
            }
            /* |rate_el - rate_hl| < 1e-10 */
            if(fabs(b_rate_el[0] - b_rate_hl[0]) >= 1.0e-10) {
                fprintf(
                    stderr,
                    "FAIL B2: |rate_el-rate_hl|=%.3e >= 1e-10\n",
                    fabs(b_rate_el[0] - b_rate_hl[0])
                );
                err = 1;
            }
            /* total_recomb_rate > 0 */
            if(b_total <= 0.0) {
                fprintf(
                    stderr,
                    "FAIL B3: total_recomb_rate=%.17g <= 0\n",
                    b_total
                );
                err = 1;
            }
        }

        /* MPI broadcast rank-0 results and verify all ranks agree */
        double bcast_focc = b_f_occ[0];
        double bcast_total = b_total;
        MPI_Comm en_comm_c = MPI_Comm_f2c(en_comm);
        MPI_Bcast(&bcast_focc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&bcast_total, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(fabs(b_f_occ[0] - bcast_focc) > 1.0e-11) {
            fprintf(
                stderr,
                "FAIL B4 rank %d: f_occ=%.17g vs bcast=%.17g\n",
                rank,
                b_f_occ[0],
                bcast_focc
            );
            err = 1;
        }
        if(fabs(b_total - bcast_total) > 1.0e-11) {
            fprintf(
                stderr,
                "FAIL B5 rank %d: total=%.17g vs bcast=%.17g\n",
                rank,
                b_total,
                bcast_total
            );
            err = 1;
        }

        negf_destruct_libnegf(hand);
        negf_destruct_session(hand);
    }

    MPI_Finalize();

    if(rank == 0 && err == 0) {
        printf("Task 7+8 C test: PASSED\n");
    }

    return err;
}
