 BEGIN{
  xs=0.0;ys=0.0;n=0.0;
  x2s=0.0;xys=0.0;
  if(xmax == 0)xmax =  9e99;
  if(xmin == 0)xmin = -9e99;

  printf("Linear regression, (xmin,xmax):(%g,%g)\n\n",xmin,xmax);
}
{
  if(($1 != "")&&((substr($1,1,1)!="#"))&&($2 != "")&&
     ($3 != -1)&&($1 <= xmax)&&($1 >= xmin))
    {
      x=$1;
      y=$2;
      xs += x;
      ys += y;
      xys += (x*y);
      x2s +=(x*x);
      n += 1.0;
    };
}
END{
  if(n == 0)
    {printf("No matching data points found.\n\n");}
  else
    {
      printf("Used %g data pairs.\n\n",n);
      xm = xs/n;
      ym = ys/n;
      a1 = (xys - ym * xs) / (x2s - xm * xs);
      a0 = ym - a1 * xm;
      if(a0 >= 0)
	printf("y = x * %g + %g\n\n",a1,a0);
      else
	printf("y = x * %g - %g\n\n",a1,-a0);
    }
}
