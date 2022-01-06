BEGIN{
    

}
{
    if(NR==1){
	lmax=$1;
	l=0;
	m=0;
    }else{
	if(m!=0){		# A and B term
	    tmp += $1**2 + $2**2;
	    if(psimple)
		printf("%5i %5i %11g %11g\n",l,m,$1,$2);
	    else
		printf("%04i\t%04i\t%20g\tA@_{%02i}^{%02i}\n%04i\t%04i\t%20g\tB@_{%02i}^{%02i}\n",l,m,$1,l,m,l,m,$2,l,m);
	    m++;
	    if(m>l){
		if(ppower)
		    printf("# power for degree %i: %11g with normalization: %11g\n",l,tmp,tmp/(2*l+1));
		m=0;
		l++;
	    }
	}else{
	    tmp = $1**2;	# only A term for m = 0
	    if(psimple)
		printf("%5i %5i %11g %11g\n",l,m,$1,$2);
	    else
		printf("%04i\t%04i\t%20g\tA@_{%02i}^{%02i}\n",l,m,$1,l,m);
	    m++;
	    if(m>l){
		m=0;
		l++;
	    }
	}
    }
    
}
