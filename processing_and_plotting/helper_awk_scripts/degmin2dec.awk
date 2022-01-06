# 35°50.75'S, 18°05.79'E
BEGIN{
    FS=",";
}
{
    lat=$1;lon=$2;
    
    split(lat,lata,"°");split(lata[2],lataa,"'");
    split(lon,lona,"°");split(lona[2],lonaa,"'");

    lon_out = lona[1] + lonaa[1]/60;
    if(lonaa[2]=="W")
	lon_out = -lon_out;
    
    lat_out = lata[1] + lataa[1]/60;
    if(lataa[2]=="S")
	lat_out = -lat_out;

    #print($1,$2)
    #print(lona[1],lonaa[1],lonaa[2],lata[1],lataa[1],lataa[2])
    print(lon_out,lat_out)
    #print("")
}