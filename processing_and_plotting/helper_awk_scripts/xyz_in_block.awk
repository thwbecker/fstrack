{
x=$8;
y=$7;
z=$9;
add=$10;

if(x< 0)x=360.0+x;

if((x >= xmin)&&(x < xmax)&&(y>=ymin)&&(y<ymax)&&(z>=zmin)&&(z<zmax))
     print(x,y,z,add);




}


