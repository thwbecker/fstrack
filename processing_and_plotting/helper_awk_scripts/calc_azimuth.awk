#
# calculate the azimuth between two points
#
# input:
#
# lon1 lat1 lon2 lat2 [degree]
#
# output:
#
# lon1 lat1 lon2 lat2 ... azimuth[deg]
#
# $Id: calc_azimuth.awk,v 1.1 2003/10/20 00:57:45 becker Exp becker $
#
BEGIN{
  twopi=6.2831853071795864769252867665590;
  pif=57.2957795130823;
}
{
  if((substr($1,1,1)!="#")&&(NF>=4)){
    lon[1]=$1;
    lat[1]=$2;
    lon[2]=$3;
    lat[2]=$4;
    
    rlon[1]=lon[1]/pif;
    rlat[1]=lat[1]/pif;
    rlon[2]=lon[2]/pif;
    rlat[2]=lat[2]/pif;

    slat[1]=sin(rlat[1]);clat[1]=cos(rlat[1]);
    slat[2]=sin(rlat[2]);clat[2]=cos(rlat[2]);

    azi = atan2(sin(rlon[2]-rlon[1])*clat[2],
		clat[1]*slat[2]-slat[1]*clat[2]*cos(rlon[1]-rlon[2]));
    if(azi < 0)
      azi += twopi;


    for(i=1;i<=NF;i++)
      printf("%s ",$i);

    printf("%g\n", azi* pif);

	
  }
}
