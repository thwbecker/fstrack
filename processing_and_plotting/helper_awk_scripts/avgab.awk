BEGIN{
}
{
  if(NR!=1){
    l=$1;
    m=$2;
    ap[l,m]=$3;
    bp[l,m]=$4;
    if($5 != ""){
      second_set=1;
      at[l,m]=$5;
      bt[l,m]=$6;
    }
  }else{
    lmax=$1;
  }
}
END{
  for(l=0;l<=lmax;l++){
    sumpa=0.0;
    sumpb=0.0;
    sumta=0.0;
    sumtb=0.0;
    for(m=0;m<=l;m++){
      sumpa += ap[l,m];
      sumpb += bp[l,m];
      if(second_set){
	sumta += at[l,m];
	sumtb += bt[l,m];
      }
    }
    if(l!=0){
      sumpa/=(l+1);
      sumpb/=(l);
      if(second_set){
	sumta/=(l+1);
	sumtb/=(l);
      }
    }
    if(second_set)
      print(l,(sumpa+sumpb)/2,(sumta+sumtb)/2,
	    sumpa,sumpb,sumta,sumtb);
    else
      print(l,(sumpa+sumpb)/2,"-",
	    sumpa,sumpb,"-","-");
  }
}
