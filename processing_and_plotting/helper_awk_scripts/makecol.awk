{
if(NR==1){
  center=(-$2)/(($3-$2));
  print($0);

}
}
END{
  for(y=0.;y<255.;y+=1.){
    x=y/255.;
    r=y;
    g  = 20.*(x-center);
    g *= g;
   
    g  = exp(-g);
    if(x<center)
      b=1.0-g;
    else   
      b=0.0;
    if(x>center)
      r=1.0-g;
    else
      r=0.0;
    b*=255;
    r*=255;
    g*=255;
    
    a=255;
    printf("%i %i %i %i\n",r,g,b,a);
  }
}
