{
x=$1;
if(x<=-0.03){
  y = 0.027 + x * 1.857;
}else{
  y = 0.003 + x * 0.667;
}
printf("%g %g ",x,y);
for(i=2;i<=NF;i++)
     printf("%s ",$i);
     printf("\n");
}
