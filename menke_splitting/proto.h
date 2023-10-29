/* ah_cross_conv_1.c */
int main(int, char *[]);
/* ah_cross_conv_1_group.c */
int main(int, char *[]);
/* ah_cross_conv_1t_group.c */
int main(int, char *[]);
int single_layer_anisotropy(double, double, double, double, double [3][3][3][3], double, double, double, double, double *, double *, double *, double *, double *, double *, int);
int euler(double, double, double, double [3][3]);
int eqn_of_plane(double [3], double [3], double [3], double [4]);
int eqn_of_plane2(double [3], double [3], double [4]);
double dot(double [3], double [3]);
int cross(double [3], double [3], double [3]);
int ray_exits_layer(double [4], double [3], double [3], double [3], double [3], double *);
int normal_to_plane(double [4], double [3]);
int rotate3x3x3x3(double [3][3][3][3], double [3][3], double [3][3][3][3]);
int cmul(double, double, double, double, double *, double *);
int load3x3x3x3(char *, double [3][3][3][3], double *, double [3], double [3], double [3]);
int isevanescent(double [3][2]);
int slowness(double [3][3][3][3], double, double, double, double [3][2], double [3][2], double [3][3][2], double [3][3][2]);
int orthogonalize(double [3][2], double [3][2]);
int sortroots(double [6][2], double [3][2], double [3][2]);
int cmproot(double [2], double [2]);
int findroot(double [10][2], int, double [2], double, double *);
int evalpoly(double [10][2], int, double [2], double [2], double [2]);
int monodiv(double [10][2], int, double [2][2], double [10][2], int *);
int cdiv(double, double, double, double, double *, double *);
int findp(double [3][3][3][3], double, double [3][2], double [3][2]);
double cdot(double [3][2], double [3][2]);
int ccross(double [3][2], double [3][2], double [3][2]);
double myrand(void);
int cnormalize(double [3][2]);
int cgauss(double *, double *, int, int, double, int *, int);
int issamecvec(double [3][2], double [3][2]);
int detpoly(double [3][3][3], double [7]);
int sumdetpoly(double [3][3][3], double [7], int, int, int, int, int, int, double);
int one_axis_tensor(double, double, double, double, double, double, double [3][3][3][3]);
int lctobc(double [6][6], double [3][3][3][3]);
int ABCDEtolc(double, double, double, double, double, double [6][6]);
int euler_blerb(void);
int write3x3x3x3(FILE *, double [3][3][3][3]);
int get_velocity(double [3][3][3][3], double, double, double, double *, double *, double *);
int do_MM(double [3][3][3][3], double [3], double [3][3]);
int eigen3(double [6], double [3], double [3][3]);
int tangent(double, double, double [3]);
int tangentr(double, double, double [3]);
int viresp(double, double, double [2], double [2]);
int scatter_fs(double [3][3][3][3], double, double [3][3][3][3], double, int, double [3][2], double [3][2], double [2], double [6][3][2], double [6][3][2], double [6][2]);
int scaleeqn3x3(double [3][3][2], double [3][2]);
/* ah_cross_conv_1_unk_pol.c */
int main(int, char *[]);
/* ah_cross_conv_2.c */
int main(int, char *[]);
/* ah_cross_conv_2_group.c */
int main(int, char *[]);
/* ah_cross_conv_spectoseis.c */
void ah_cross_conv_spectoseis_(float *, float *, float *, int *, float *, float *, float *);
void ah_cross_conv_spectoseis(float *, float *, float *, int *, float *, float *, float *);
/* ah_cross_conv_spectoseis_driver.c */
int main(int, char *[]);
/* ah_est1pulse.c */
int main(int, char *[]);
/* ah_est2pulse.c */
int main(int, char *[]);
double lagdot(float *, int, float *, int, int);
/* ah_est3pulse.c */
int main(int, char *[]);
double lagdot(float *, int, float *, int, int);
/* ah_est4pulse.c */
int main(int, char *[]);
double lagdot(float *, int, float *, int, int);
/* ah_est5pulse.c */
int main(int, char *[]);
double lagdot(float *, int, float *, int, int);
/* gauss.c */
int gauss(double *, double [], int, int, double, int *, int);
/* menke_distaz.c */
int menke_distaz(float, float, float, float, float *, float *, float *);
/* read_sac.c */
int read_sac(ahhed *, float **, char *);
/* sac_cross_conv_1.c */
int main(int, char **);
/* sac_test.c */
int main(int, char **);
/* SHahio.c */
int gethead(ahhed *, FILE *);
int puthead(ahhed *, FILE *);
int size(ahhed *);
int tohead(int, FILE *);
int getdata(ahhed *, char *, FILE *);
int putdata(ahhed *, char *, FILE *);
int putrecord(ahhed *, char *, FILE *);
int getrecord(ahhed *, char *, FILE *);
int getrecord2(ahhed *, char **, FILE *);
int gogethead(int, ahhed *, FILE *);
int gogetrecord(int, ahhed *, char *, FILE *);
int logger(char *, ahhed *);
int out_is_tty(void);
int in_is_tty(void);
char *mkdatspace(ahhed *);
int get_null_head(ahhed *);
int acpy(char *, char *, unsigned);
int ah_error(char *, char *, int);
int maxamp(ahhed *, char *);
int xdr_gethead(ahhed *, XDR *);
int xdr_puthead(ahhed *, XDR *);
int xdr_tohead(int, XDR *);
int xdr_getdata(ahhed *, char *, XDR *);
int xdr_putdata(ahhed *, char *, XDR *);
int xdr_putrecord(ahhed *, char *, XDR *);
int xdr_getrecord(ahhed *, char *, XDR *);
int xdr_getrecord2(ahhed *, char **, XDR *);
int xdr_gogethead(int, ahhed *, XDR *);
int xdr_gogetrecord(int, ahhed *, char *, XDR *);
int xdr_ahhead(XDR *, ahhed *);
