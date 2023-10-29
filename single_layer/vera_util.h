/* 

   header files for Schulte-Pelkum's tensor and splitting routines

*/
#define VERA_PREC double


#define vera_sav_to_cijkl_ftrn vera_sav_to_cijkl_ftrn_
#define vera_print_cijkl_ftrn vera_print_cijkl_ftrn_
#define vera_print_cijkl_file_ftrn vera_print_cijkl_file_ftrn_
#define vera_print_splitting_ftrn vera_print_splitting_ftrn_
#define vera_layer_split_from_tensor_ftrn vera_layer_split_from_tensor_ftrn_
#define vera_read_cijkl  vera_read_cijkl_
extern void vera_sav_to_cijkl_ftrn(VERA_PREC *,double *,VERA_PREC *);
extern void vera_print_cijkl_ftrn(VERA_PREC *,int *);
extern void vera_layer_split_from_tensor_ftrn(VERA_PREC *,VERA_PREC *,VERA_PREC *,VERA_PREC *,
					      VERA_PREC *,VERA_PREC *,VERA_PREC *,VERA_PREC *,VERA_PREC *);
extern void vera_print_splitting_ftrn(VERA_PREC *,VERA_PREC *,VERA_PREC *,VERA_PREC *,VERA_PREC *,VERA_PREC *,int *);
extern void vera_print_cijkl_file_ftrn(VERA_PREC *);
extern void vera_read_cijkl(VERA_PREC *,int *);
