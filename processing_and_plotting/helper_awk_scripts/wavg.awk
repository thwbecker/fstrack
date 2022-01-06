#
#
# perform weighted average
#
# weight .... $col
#
# col defines the column to use
#
#
# avg the type of average:
#
# 1: arithmetic
# 2: geometric, log-average
# 3: harmonic, 1/val average
#
# A >= G >= H
#
BEGIN{
  if(avg==0)
    avg = 1;
  if(col==0)
    col = 2;
  x = 0.0;
  n = 0;
# ofmt="%20.10e\n";
  ofmt="%g\n";
  ws = 0.;
}
{
  if((substr($0,1,1)!="#")&&(NF>=col)){
      w = $1;
    if(avg == 1){		# Arithmetic
      x += w * $(col);
    }else if(avg == 2){		# Geometric
      if($(col) <= 0)
	print("avg.awk: error, for geometrical mean: value ",n+1," is ",$(col)) > "/dev/stderr";
      x += w * log($(col)); 
    }else if(avg == 3){		# Harmonic
      if($(col) == 0)
	print("avg.awk: error, for harmonic mean: value ",n+1," is ",$(col)) > "/dev/stderr";
      x += w/$(col); 
    }else{
      print("avg.awk: error, avg code ",avg," undefined") > "/dev/stderr";
    }
    ws += w;
    n++;
  }
}
END{
  if(avg == 1){			# artihmetic
    printf(ofmt,x/ws);
  }else if(avg == 2){		# gemoetric
    printf(ofmt,exp(x/ws));
  }else if(avg == 3){
    printf(ofmt,ws/x);
  }else{
    print("avg.awk: error, avg code ",avg," undefined") > "/dev/stderr";
  }
}
