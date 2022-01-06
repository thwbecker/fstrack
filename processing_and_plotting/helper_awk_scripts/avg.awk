#
#
# perform average
#
# col defines the column to use
#
#
# avg the type of average:
#
# 1: arithmetic x = c_1 x_1 + c_2 x_2
# 2: geometric, log-average, x = x_1^c_1 x_2^c_2 , log x = c_1 log x_1 + c_2 log x_2
# 3: harmonic, 1/val average, 1/x = (c_1/x_1 + c_2/x_2)
#
# A >= G >= H
#
BEGIN{
    if(avg=="")
	avg = 1;
    if(col=="")
	col = 1;
    sum = 0.0;
    n = 0;
# ofmt="%20.10e\n";
    if(ofmt=="")
	ofmt="%g\n";
}
{
    if((substr($1,1,1)!="#")&&(NF>=col)&&(!match(tolower($col),"nan"))){
	if(avg == 1){		# Arithmetic
	    sum += $(col);
	}else if(avg == 2){		# Geometric
	    if($(col) <= 0)
		print("avg.awk: error, for geometrical mean: value ",n+1," is ",$(col)) > "/dev/stderr";
	    sum += log($(col)); 
	}else if(avg == 3){		# Harmonic
	    if($(col) == 0)
		print("avg.awk: error, for harmonic mean: value ",n+1," is ",$(col)) > "/dev/stderr";
	    sum += 1.0/$(col); 
	}else{
	    print("avg.awk: error, avg code ",avg," undefined") > "/dev/stderr";
	}
	n++;
    }
}
END{
    if(avg == 1){			# artihmetic
	printf(ofmt,sum/n);
    }else if(avg == 2){		# gemoetric
	printf(ofmt,exp(1/n*sum));
    }else if(avg == 3){           # harmonic
	printf(ofmt,n/sum);
    }else{
	print("avg.awk: error, avg code ",avg," undefined") > "/dev/stderr";
    }
}
