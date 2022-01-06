BEGIN{
  n=0;
  srand();

}
{
  if(substr($1,1,1)!="#"){
    n++;
    x = rand();
    s[n] = sprintf("%.7f\t%s",x,$0);
  }
}
END{
  asort(s);
  for(i=1;i<=n;i++){
    m = split(s[i], a ," ");
    for(j=2;j<=m;j++)
      printf("%s ",a[j]);
    printf("\n");
  }

}
