BEGIN{
  if(col==0)
    col=1;
  nc=0;
}

{
    if((substr($1,1,1)!="#") && (NF >= col) && (tolower($(col)) != "nan")){
	nc++;
	if(nc==1){
	    min=$col;
	    max=$col;
	}
	if($(col) < min){
	    min=$(col);
	    pxmin=$1;
	}
	
	if($(col) > max){
	    max=$(col);
	    pxmax=$1;
	}
    }
    
}
END{
  if(!px)
    printf("%g %g\n",min,max);
  else
    printf("%g %g %g %g \n",min,max,pxmin,pxmax);
}
