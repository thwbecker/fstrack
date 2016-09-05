/* advect_and_search_for_state.c */
int advect_and_search_for_state(struct mod *, int, struct stt *, struct stt *, unsigned short, unsigned short, double *);
void iter_report_progress(struct stt *, int, int, int, int, int, double, double, double, unsigned short *);
/* advect.c */
void advect(struct trc *, struct mod *, int, int);
unsigned short rk_wrapper(double *, double, double, double, unsigned short, unsigned short, double, struct der_par *, unsigned short);
void compute_state_parameters(struct stt *, struct der_par *, char *);
void deal_with_discarded_tracer(struct trc *, double *);
void texture_pathline_init(struct der_par *, double);
/* average_rphi_tracers.c */
/* average_tracers.c */
double averaging_weight(double, double, double, int, int, struct dw, unsigned short, double, double);
double interpolate_dw(double, struct dw);
void read_weights_file(struct dw *, char **);
unsigned short colocated(struct str *, struct str *, double, double);
/* bailoutcrit.c */
unsigned short bailout_crit_fulfilled(struct stt *, int, struct mod *, double *, unsigned short *, double *);
/* calc_divergence.c */
void calc_divergence(struct mod *);
/* calc_isa.c */
void calc_isa(struct trc *, int, struct mod *);
void calc_pipar(double *, double, double *, double *, double *, struct der_par *);
/* calc_lpo_from_streamline.c */
/* calc_lyapunov.c */
void calc_lyapunov(struct mod *, int);
/* calc_strain.c */
void calc_vel_and_strain_rate_units(double, double *, struct der_par *, double *, double *, unsigned short);
void calc_strainrate_and_vorticity(double *, double *, double *, unsigned short);
void calc_vel_and_strain_rate(double, double *, struct der_par *, double *, double *, unsigned short);
double calc_max_dt_from_strain_rate(double, double *, struct der_par *);
double char_strain_rate_from_vg(double *);
double char_strain_rate_from_e(double *, double *);
void calc_el_strain(double *, double *);
void calc_el_strain_(double *, double *);
void calc_cauchy_strain(double *, double *);
void calc_l2_left_stretch_(double *, double *);
void calc_l2_left_stretch(double *, double *);
double max_strain_from_def(double *, double *);
double max_strain_from_left_stretch(double *, double *);
double mstrain(double *);
/* cart_tester.c */
void cyl_strain_out(FILE *, double *, double, int);
void cart_to_shift_cyl(double *, double *, int);
void strain_eigen_out(FILE *, double *, double);
void advance_argument(int *, int, char **);
void help(char **);
/* cFfromdiscreteG.c */
/* cijklrotate.c */
/* cornerflow.c */
void calc_cornerflow_constants(int, double, double, double *, double *, double *, double *, char **);
double stream(double *, double, double, double, double);
void cylvel(double *, double *, double, double, double, double);
void cyl2cart(double *, double *);
void cart2cyl(double *, double *);
void cylvec2cartvec(double *, double *, double *);
/* cvec2ellipsoid.c */
/* datafit_util.c */
void svdfit_driver(double *, double *, double *, int, int, double, void (*)(void), int, double *, double *, double *);
int lfit_driver(double *, double *, double *, int, int, void (*)(void), int, double *, double *, double *);
void fit_poly(struct nr_datas *, int, int, double, double *, double *, double *, unsigned short, unsigned short, double *);
double evaluate_model(double, int, double *, void (*)(void), int);
void poly_fit_func(double, double *, int, int);
double poly_val(double, double *, int);
int red_dof(int, int);
int nr_lfit(double [], double [], double [], int, double [], int [], int, double **, double *, void (*)(void), int);
void nr_svdfit(double [], double [], double [], int, double [], int, double **, double **, double [], double *, double, void (*)(void), int);
/* deriv.c */
void rk_derivs(double, double *, double *, int, struct der_par *);
void fse_derivs_wrapper(double, double *, double *, int, struct der_par *, double *, double *, unsigned short, unsigned short, unsigned short *, unsigned short);
void init_der_par(struct der_par **);
void rk_check_phys_limit(double *, int, struct der_par *, unsigned short);
void check_phys_lim_tracer__(double *, double *);
void check_phys_lim_tracer_(double *, double *);
void check_phys_lim_tracer(double *, double *);
void check_physics_based_error(double *, double *, int *, double *, double *, double *);
int determine_n_error(int, struct der_par *);
void determine_error_scale(double, double *, double *, double, int, struct der_par *, double *);
/* extract_strain_field.c */
/* fazi2splitstat.c */
/* generate_vgm.c */
void generate_vgm_plate(double *, double, double, double, double);
void generate_vgm_sw(double *, double, double, double);
/* indexx.c */
void indexx(int, double *, int *);
/* init.c */
void initialize_tracers(struct mod *);
void allocate_and_clear_tracers(struct mod *, int);
void clear_tracer(struct trc *);
void init_tracer_depths(struct mod *);
unsigned short close_to_pole(double *);
void check_tracer_location(double *);
/* input_para.c */
void set_defaults(struct mod *);
void check_for_parameters(int, char **, struct mod *);
void advance_argument(int *, int, char **);
void phelp(char **, struct mod *);
/* invert3x3c.c */
void invert3x3c(double *, double *);
void invert3x3c_(double *, double *);
/* linalg.c */
void calc_eigensystem_sym(double *, double *, double *, unsigned short);
void calc_eigensystem_sym_(double *, double *, double *, int);
void calc_eigensystem(double *, double *, double *, double *, unsigned short, unsigned short *);
void eigen2d(double, double, double, double *, double *, double *);
void calc_sqrt_sym_matrix(double *, double *);
void assemble_sym_matrix(double *, double *, double *);
void calc_a_dot_at(double *, double *);
void calc_a_dot_at_3x3(double [3][3], double [3][3]);
void calc_at_dot_a(double *, double *);
void calc_cd_symm_part(double *, double *);
void calc_cd_asymm_part(double *, double *);
double det3x3(double *);
double trace3x3(double *);
double sec_inv3x3(double *);
double double_dot(double *, double *);
void calc_a_times_b_mat(double *, double *, double *);
void calc_a_times_b_3x3mat(double [3][3], double [3][3], double [3][3]);
double vector_diff(double *, double *, int);
void zero_vector(double *, int);
void scale_vector(double *, double, int);
void scale_vector3d(double *, double);
void copy_vector(double *, double *, int);
void copy_3dvector(double *, double *);
void copy_3dmatrix(double *, double *);
void add_a_to_b_vector(double *, double *, int);
void sub_a_from_b_vector(double *, double *, int);
void a_equals_b_minus_c_vector(double *, double *, double *, int);
double vec3ddotp(double *, double *);
double vec2ddotp(double *, double *);
double norm3d(double *);
void normalize3d(double *);
double vec_sum(double *, int);
void assign_c_for_ab_hs(double *, double **, int);
void symmat9to3x3(double *, double [3][3]);
void mat9to3x3(double *, double [3][3]);
void mat3x3to9(double [3][3], double *);
double max_vec(double *, int);
double max_abs_vec(double *, int);
void a_equals_b_vector(double *, double *, int);
void a_equals_b_vector3d(double *, double *);
void unity_matrix(double *, int);
void a_is_b_transpose_3by3(double *, double *);
void transpose_3by3_inplace(double *);
void a_is_b_transpose_nbyn(double *, double *, int);
double vec2d_crossp(double *, double *);
void transpose_nbyn_inplace(double *, int);
double norm(double *, int);
double vecdotp(double *, double *, int);
void zero_small_entries(double *, int);
double remove_trace(double *);
/* linalg_lapack.c */
void calc_exp_matrixt(double *, double, double *);
void test_null_space(double *, double **, int, int *, double);
void find_ABzero(double *, double **, int, int *, double);
/* main.c */
/* make_random_tensors.c */
/* make_var_tensor.c */
/* misc.c */
FILE *myopen(char *, char *, char *);
double save_sqrt(double);
void my_vecalloc(double **, int, char *);
void my_svecalloc(float **, int, char *);
void my_vecvalloc(float **, int, char *);
void my_ivecalloc(int **, int, char *);
void my_vecrealloc(double **, int, char *);
double myrand(long *);
double ran1(long *);
double ran2(long *);
double gasdev(long *);
double myrandnr(double, long *);
double mygauss_randnr(double, long *);
void get_lin_int_weights(double, double *, int, int *, int *, double *, double *);
double lin_inter(double, double *, double *, int);
int read_sym_6by6(double *, FILE *);
double restrict_ani_factor(int, double);
/* my_mem_util.c */
/* nr_util.c */
int *nr_ivector(long, long);
void nr_free_matrix(double **, long, long, long, long);
double **nr_matrix(long, long, long, long);
double *nr_vector(long, long);
void nr_error(char *);
void nr_free_vector(double *, long, long);
int nr_gaussj(double **, int, double **, int);
/* odeint.c */
int nr_odeint(double *, int, double, double, double, double, double, double, int *, int *, struct der_par *);
int nr_rkqs(double *, double *, int, double *, double, double, double, double *, double *, double *, struct der_par *);
void nr_rkck(double *, double *, int, double, double, double *, double *, struct der_par *);
int odeint_exit(double **, double **, double **, double, struct der_par *, int);
/* output.c */
void write_ascii_layer(double *, int, struct mod *, char *, char *, double);
void write_ascii_layer_v(float *, int, struct mod *, char *, char *, double);
void write_tracer_history(struct trc *, char *, int, unsigned short, struct mod *);
void write_all_tracer_history(struct mod *, char *, unsigned short);
void print_tracer_data(struct trc *, int, FILE *, unsigned short, char *, struct mod *);
void write_tracer_field(struct mod *, int, int, int, char *, int, unsigned short);
void print_tracer_stats(struct trc *, int, FILE *);
void print_llz(struct stt *, FILE *);
void print_llzt(struct stt *, FILE *);
void print_tllz(struct stt *, FILE *);
void process_ti_tens_frac(double *, double *, unsigned short, FILE *);
void rotate_mat_to_reg_cart_system(double *, double *, double *);
/* output_simple.c */
void print_vector(double *, int, FILE *);
void print_vector_stderr__(double *, int *);
void print_vector_wl(double *, int, FILE *, char *);
void print_sym_6by6(double *, FILE *);
void print_6by6_nice(double *, FILE *);
void print_6by6_hprec(double *, FILE *);
void print_cij_vera_sorted(double *, FILE *);
/* plot_kernel.c */
/* pmodel.c */
double pressure_model(double);
/* polvgm2cartvgm.c */
/* prem_util.c */
int prem_find_layer_x(double, double, double *, int, int, double *);
double prem_compute_pval(double *, double *, int, double);
double prem_compute_dpval(double *, double *, int, double);
double prem_vs_voigt(double, double, double, double, double);
void prem_get_rhodrho(double *, double *, double, struct prem_model *);
void prem_get_rho(double *, double, struct prem_model *);
void prem_get_pressure(double *, double, struct prem_model *);
void prem_get_values(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, struct prem_model *);
int prem_read_model(char *, struct prem_model *, unsigned short);
int prem_read_para_set(double *, int, int, FILE *);
/* readgrds.c */
void read_vel_grids(struct mod *, unsigned short, unsigned short);
void find_max_horizontal_divergence(struct mod *, double);
void resort_and_check(float *, float *, double *, int, int, unsigned short, double, unsigned short, unsigned short, float);
void read_depth_levels(struct mod *, int **);
void read_time_intervals(struct mod *);
/* sav2afactor.c */
/* sav2average.c */
/* sav2cijkl.c */
/* sav2decompose.c */
void print_tens_decomp(double *, unsigned short, int);
void compute_best_fit_hex(double *, double *, double *, int, double *);
/* sav2rotate.c */
/* sav2splitting.c */
/* savstrain2stress.c */
/* sens_handling.c */
void compute_phi_terms_from_Sav(double *, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, struct mod *, int);
void interpolate_sens(double, double, struct swss *, double *, double *, double *, double *, double *);
void read_sw_sens(struct mod *);
void compute_azideg_amp_from_afactors(double *, double *, double *);
/* series_analyze.c */
void calc_mean_and_stddev(float *, float *, int, double *, double *, double *, unsigned short, unsigned short, float *);
double mean(double *, int);
double wmean(double *, double *, int);
double std(double *, int);
void stat_moments(double *, int, double *, double *, double *, double *, double *, double *, unsigned short);
void fit_harmonic(double *, double *, double *, int, int, double, double *, double *, double *, unsigned short, unsigned short, unsigned short);
void harm_fit_func(double, double *, int, int);
/* spline.c */
void nr_spline(double *, double *, int, double, double, double *);
void nr_splint(double *, double *, double *, int, double, double *, unsigned short, double *);
/* splitting_util.c */
void analyze_splitting_azi_dep(int, int, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, unsigned short, int);
/* state_handling.c */
void add_state(struct trc *, double, struct der_par *);
void init_state(struct stt *, double, struct der_par *);
void clear_state(struct stt *, double);
void copy_state(struct stt *, struct stt *, unsigned short, char *, char *);
/* svd_util.c */
void nr_svdcmp(double **, int, int, double [], double **);
void nr_covsrt(double **, int, int [], int);
double nr_pythag(double, double);
void nr_svdvar(double **, int, double *, double **);
void nr_svbksb(double **, double [], double **, int, int, double [], double []);
/* test_stuff.c */
/* tmodel.c */
double temperature_model(double, int);
/* tracerl2cevec.c */
/* trig.c */
void rtp2xyz(double *, double *);
void xyz2rtp(double *, double *);
void polar_vec2cart_vec_at_xp(double *, double *, double *);
void cart_vec2polar_vec_at_xp(double *, double *, double *);
void polar_vec2cart_vec_at_xp_(double *, double *, double *);
void polar_to_cart_mat_at_r(double *, double *, double *);
void cart_to_polar_mat_at_r(double *, double *, double *);
void rotate_cart_mat(double *, double *, double, double, double);
void polar_to_cart_mat_at_r_(double *, double *, double *);
void cart_to_polar_mat_at_r_(double *, double *, double *);
void polar_to_cart_mat_at_r3x3(double [3][3], double [3][3], double *);
void cart_to_polar_mat_at_r3x3(double [3][3], double [3][3], double *);
void rotate_cart_mat_3x3(double [3][3], double [3][3], double, double, double);
void calc_rotmat_cart(double [3][3], double, double, double);
void calc_polar_basis_at_r(double [3][3], double *);
void calc_cart_basis_at_r(double [3][3], double *);
void my_sincos(double, double *, double *);
void my_sincosd(double, double *, double *);
void rotate_mat(double [3][3], double [3][3], double [3][3]);
void abase_vec2bbase_vec(double *, double *, double [3][3]);
void calc_Ax3x3(double [3][3], double *, double *);
void calc_xy_from_azi_deg(double, double *, double *);
void calc_azi_from_axay(double, double, double *);
double azi_from_etep(double, double);
double dip_from_er_unityvec(double);
void fix_deg_angle(double *);
double dist_on_sphere(double *, double *);
void lonlatz_from_xp(double *, double *);
void xp_from_lonlatz(double, double, double, double *);
void diff_from_orient_mean(double *, int, double *, double *, double *, int, double *);
double diff_orient_angle_azid(double, double);
double diff_orient_vector_azid(double *, double *, unsigned short);
