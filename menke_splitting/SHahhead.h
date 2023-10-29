/* SHahhead.h */
/*	structure for data file header	--	witte, 11 June 85	*/

#define CODESIZE 6
#define CHANSIZE 6
#define STYPESIZE 8
#define COMSIZE 80
#define TYPEMIN 1
#define TYPEMAX 6
#define LOGSIZE 202
#define LOGENT 10
#define NEXTRAS 21
#define NOCALPTS 30

typedef	struct	{
		float	x;
		float	y;
		} vector;

typedef	struct	{
		float	r;
		float	i;
		} complex;

typedef struct  {
		double r;
		double i;
		} d_complex;

typedef	struct	{
		float	xx;
		float	yy;
		float	xy;
		} tensor;

struct	time	{
		short		yr;	/* year		*/
		short		mo;	/* month	*/
		short		day;	/* day		*/
		short		hr;	/* hour		*/
		short		mn;	/* minute	*/
		float		sec;	/* second	*/
		};

struct	calib	{
		complex		pole;	/* pole		*/
		complex		zero;	/* zero		*/
		};

struct	station_info	{
		char		code[6];	/* station code		*/
		char		chan[6];	/* lpz,spn, etc.	*/
		char		stype[8];	/* wwssn,hglp,etc.	*/
		float		slat;		/* station latitude 	*/
		float		slon;		/*    "    longitude 	*/
		float		elev;		/*    "    elevation 	*/
		float		DS;	/* gain	*/
		float		A0;	/* normalization */
		struct	calib	cal[NOCALPTS];	/* calibration info	*/
		};

struct	event_info	{
		float		lat;		/* event latitude	*/
		float		lon;		/*   "   longitude	*/
		float		dep;		/*   "   depth		*/
		struct	time	ot;		/*   "   origin time 	*/
		char		ecomment[80];	/*	comment line	*/
		};

struct	record_info	{
		short		type;	/* data type (int,float,...) 	*/
		long		ndata;	/* number of samples		*/
		float		delta;	/* sampling interval		*/
		float		maxamp;	/* maximum amplitude of record 	*/
		struct	time	abstime;/* start time of record section */
		float		rmin;	/* minimum value of abscissa 	*/
		char		rcomment[80];	/* comment line		*/
		char		log[LOGSIZE]; /* log of data manipulations */
		};

typedef struct {
		struct	station_info	station;	/* station info */
		struct	event_info	event;		/* event info	*/
		struct	record_info	record;		/* record info	*/
		float		extra[NEXTRAS];	/* freebies */
		} ahhed;


#define	FLOAT	1
#define	COMPLEX	2
#define	VECTOR	3
#define	TENSOR	4
#define	DOUBLE	6

#define AHHEADSIZE sizeof(ahhed)

