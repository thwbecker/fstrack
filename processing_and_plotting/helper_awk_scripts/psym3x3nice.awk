# print a symmetric 3x3 matrix in upper right hand side format nicely
BEGIN{

}
{
  k=1;
  for(i=1;i<=3;i++)
    for(j=i;j<=3;j++){
      a[(i-1)*3+j] = $k;
      k++;
    }
  printf("%9.2e %9.2e %9.2e \n",a[1],a[2],a[3]);
  printf("          %9.2e %9.2e \n",a[5],a[6]);
  printf("                    %9.2e \n",a[9])

}
