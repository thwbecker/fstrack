#
# read two gephart output files and print best ratio/vector
#
BEGIN{
  read_best = 0;
  ntest = 1;
  neigm=9;
  nvtest=11;
  neigr[1]=0;

  if(eps == "")
    eps=1e-4;
}
{
  if(match(FILENAME,"gephartfmex5")){ # read best fit value from first line
    if(NR==1){
      best_fit = $1;
      read_best = 1;
    }
  }
  if(match(FILENAME,"gepfmexcr5.out")){
    # read best fit values for different R and vectors
    if(!read_best)
      print("best fit value was not read, list mex file first") > "/dev/stderr";
    if((NF==10)&&($10==0)&&($9==0)&&($8==0)&&($7==0)){		
      # eigenvectors
      neigr[ntest]++;
      if(neigr[ntest] > neigm)
	print("too many eigen vectors",neigr[ntest],neigm) > "/dev/stderr";
      ev[(ntest-1)*neigm+neigr[ntest]] = $0;
      nv=0;
    }else{
      # best fit values for nvtest different R values
      nv++;
      bvf[(ntest-1)*nvtest+nv] = $0;
      if(nv == nvtest){		# last one
	ntest++;
	neigr[ntest]=0;
      }
    }
  }
}
END{
  # output
  ntest--;
  fbf=0;
  if(test_output){
    # verbose test output for debuggin
    
    print("read ",ntest,"sets ");
    for(i=1;i <= ntest;i++){
      r=0.0;
      for(j=1;j <= nvtest;j++){
	printf("set: %2i r: %.1f nvec: %2i \t",i,r,neigr[i]);
	nfields = split(bvf[(i-1)*nvtest+j],val);
	if(nfields != neigr[i])
	  printf("\n\nmismatch: %i %i: %i %i\n\n",i,j,nfields,neigr[i]) > "/dev/stderr";
	for(k=1;k <= neigr[i];k++){
	  printf("%6.3f ",val[k]);
	  dx = val[k] - best_fit;if(dx<0)dx=-dx;
	  if(dx < 1e-5){
	    fbf++;
	    printf("best fit vector %2i: %s\n",fbf,ev[(i-1)*neigm+k]) > "/dev/stderr";
	  }
	}
	printf("\n");
	r += 0.1;
      }
      printf("\n");
    }
  
  }else{
    # only print vector for best fitting values
    for(i=1;i <= ntest;i++){
      r=0.0;
      for(j=1;j <= nvtest;j++){
	nfields = split(bvf[(i-1)*nvtest+j],val);
	if(nfields != neigr[i])
	  printf("mismatch: %i %i: %i %in",i,j,nfields,neigr[i]) > "/dev/stderr";
	for(k=1;k <= neigr[i];k++){
	  dx = val[k] - best_fit;
	  if(dx<0)dx=-dx;

	  if(dx < eps){
	    fbf++;
	    printf("set: %2i r: %.1f vec: %2i values: %s\n",
		   i,r,k,ev[(i-1)*neigm+k]);
	  }
	}
	r += 0.1;
      }
    }
  }
}