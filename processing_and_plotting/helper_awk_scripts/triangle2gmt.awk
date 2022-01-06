#
# convert traingle output files into different formats
#
# use this script as
#
#
# gawk --assigne elecenters=1 \
#  -f triangle2gmt.awk file.node file.ele file.edge [file.code] 
#
# where not all files need to be present for all output modes
# which are set by setting different variables to non-zero (see output part)
#
BEGIN{
}

{
  if(substr(FILENAME,
	    length(FILENAME)-3,length(FILENAME))=="node"){
# read in the nodal locations
    if(FNR == 1)
      nrp=$1;
    if(FNR>1 && $1!="" && (substr($1,1,1)!="#")){
      n=$1;
      x[n]=$2;
      y[n]=$3;
    }
  }else if(substr(FILENAME,
		  length(FILENAME)-2,
		  length(FILENAME)) == "ele"){
# read in the element connectivity
    if(FNR == 1){
      nrel=$1; 
      nrcon=$2;
    }
    if(FNR>1 && $1!="" && (substr($1,1,1)!="#")){
      n=$1;
      for(i=1;i<=nrcon;i++){
	c[1+(n-1)*nrcon+i-1]=$(i+1);
      }
    }
  }else if(substr(FILENAME,
		  length(FILENAME)-3,
		  length(FILENAME))=="edge"){
# read in the edges
    if(FNR == 1){
      nred=$1;
      nrbm=$2;
    }
    if(FNR>1 && $1!="" && (substr($1,1,1)!="#")){
      n=$1;
      for(i=1;i<=2;i++)
	e[1+(n-1)*2+i-1]=$(i+1);
      bc[n]=$(i+3);
    }
  }else if(substr(FILENAME,
		  length(FILENAME)-3,
		  length(FILENAME))=="code"){
# read in boundary codes
    code[FNR]=$1;
    readcode++;
  }
}
END{
  
  if(elecenters){		# print the center of the elements
    for(i=1;i<=nrel;i++){
      x0=y0=0.0;
      for(j=1;j<=3;j++){
	pt=c[1+(i-1)*nrcon+j-1];
	x0+=x[pt];
	y0+=y[pt];
      }
      print(x0/3,y0/3,i);
    }
  }else if(inbetween){		# print what again?
    for(i=1;i<=nrel;i++){
      for(j=4;j<=6;j++){
	pt=c[1+(i-1)*nrcon+j-1];
	print(x[pt],y[pt]);
      }
    }
  }else if(edges){		# print the edges
    for(i=1;i<=nred;i++){
      x0=0.0;
      y0=0.0;
      for(j=1;j<=2;j++){
	pt=e[1+(i-1)*2+j-1];
	x0+=x[pt];
	y0+=y[pt];
      }
      print(x0/2,y0/2,i);
    }
  }else if(tpolygon){		# print the polyon corners for polygon tpolygon
    print(nrel);
    for(i=1;i<=nrel;i++){
      print(4);
      for(j=1;j<=3;j++){
	pt=c[1+(i-1)*nrcon+j-1];
	print(x[pt],y[pt]);
      }
      j=1;			# repeat first
      pt=c[1+(i-1)*nrcon+j-1];
      print(x[pt],y[pt]);
    }
  }else{
    if(readcode==nrel){
      for(i=1;i<=nrel;i++){
	if(code[i]==1){c1=0;c2=0;c3=255;}
	else if(code[i]==2){c1=255;c2=50;c3=0;}
	else if(code[i]==3){c1=255;c2=150;c3=0;}
	else if(code[i]==5){c1=255;c2=0;c3=255;}
	else {c1=200;c2=200;c3=200;}
	printf("> -G%03i/%03i/%03i\n",c1,c2,c3);
	for(j=1;j<=3;j++){
	  pt=c[1+(i-1)*nrcon+j-1];
	  printf("%g %g\n",x[pt],y[pt]);
	}
	j=1;
	pt=c[1+(i-1)*nrcon+j-1];
	printf("%g %g\n",x[pt],y[pt]);
      }
    }else{
      for(i=1;i<=nrel;i++){
	for(j=1;j<=3;j++){
	  pt=c[1+(i-1)*nrcon+j-1];
	  printf("%g %g\n",x[pt],y[pt]);
	}
	j=1;
	pt=c[1+(i-1)*nrcon+j-1];
	printf("%g %g\n",x[pt],y[pt]);
	printf(">\n");
      }
    }
  }
}
