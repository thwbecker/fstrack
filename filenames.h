//
// filenames for fstrack
//
// input file with depths of vr.i.grd and the like
#define DFILE "vdepth.dat"
// input file with initial x y tracer locations for  SPOTTED_LAT_FROM_FILE
// .dat or model->ostring will be appended 
#define TRLFILE "tracer.lonlat"
// input file with initial x y z tracer locations for  SPOTTED_3D_FROM_FILE
// .dat or model->ostring will be appended 
#define TRLZFILE "tracer.lonlatz"
// input file with initial x y z a tracer locations and attribute for  SPOTTED_3D_WITH_ATTR
#define TRLZAFILE "tracer.lonlatza"
// .dat or model->ostring will be appended 
// input file with tracer depth layers
// .dat or model->ostring will be appended 
#define TDFILE "tdepth"
// output filename for tracer L strain output, depth and .dat will be added, or as modified by model->ostring
#define TSFILE "f.s"
// output filename for tracer L ISA strain output, depth and .dat will be added, tracer prepended see above
#define TISAFILE "f.s"
// output filename for tracer F strain matrix output, depth and .dat will be added, tracer prepended see above
#define TFFILE "f.f"

// output filename for velocity stats output
// .dat or model->ostring will be appended
#define VSFILE "vel.rms"

//
// filenames for average_tracers
//
#define WSFILE "weights.dat"

// file with times 1...n at which velocities in directories 
// 1/ ... n/ are given
#define THFILE "times.dat"
/* 
   initial strain, to be added to unity matrix, [9] format
*/
#define INIT_STRAIN_FILE "initstrain.dat"
/* 

location of the fstrack type sensitivity files fsens.period.dat
which are in format z A F L (WITHOUT the starting HOME directory)


*/
#define SW_SENS_FILE "progs/src/fstrack/sw_sens/PREMb/fsens"
/* 

filename for gamma factors, similar to vr.i.grd, vt.i.grd. vp.i.grd

*/
#define GAMMA_FILE "er"

