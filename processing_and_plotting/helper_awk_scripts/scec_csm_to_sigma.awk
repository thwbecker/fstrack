#
# convert SCEC CSM format to lon-lat-z stress tensor
#                                                     9                             16
#LON LAT DEP 1 SHmax_trend SHmax_mag SHmin_mag SV_mag See Sen Seu Snn Snu Suu RATIO DEV ISO
#LON LAT DEP 2 1-sigma 1-sigma 1-sigma 1-sigma DOT dot-1sig ANG s1-1sig s2-1sig s3-1sig 1sig 1sig 1sig
#
BEGIN{

}
{
    if(substr($1,1,1)!="#"){
	if($4 == 1){
	    # stress line, as opposed to uncertainty line
	    lon=$1;lat=$2;dep=$3;
	    see=$9;  sen=$10;seu=$11;
	    snn=$12; snu=$13;
	    suu=$14;
	    dev = $18;
	    print(lon,lat,dep,see,sen,seu,snn,snu,suu);
	}
    }
}