{
n++;
x[n]=$1;
}
END{
  for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
      print(x[i],x[nx+j]);

}
