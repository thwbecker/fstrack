BEGIN{
#
# create histogram using the data in column col as a sorting criterion
#
# uses 10 bins if nbin is not set
#
# if dx is set, will use fixed width bins
#
# if cols is set, will take this column to sum each entry to the bin as 
# specified by col
#
# if min or max are set, will use those, else find automatically
#
# output is in format  
# left of bin, number of entries, fraction of  total number, 
# sum of y entries in bin, sum of y entries over total sum
#
# for cols not set, columns 4 and 5 will be identical to  2 and 3
#
# if step_output is set to unity, will repeat y_i entries at x_{i+1}
# for plotting purposes
#
#
# $Id: histogram.awk,v 1.1 2007/01/13 21:15:09 becker Exp becker $
#
# number of bins
#
  if(((nbin=="")||(nbin==0))&&((dx=="")||(dx==0))){
    nbin=10;
    dx=-1;
  }
# pick column with data selection quantity
  if((col=="")||(col==0))
    col=1;
# column with summing quantity, if set to -1 will use unity, ie. give
# same answer as bin_cnt
  if((cols=="")||(cols==0)){
    cols=-1;
  }
# find only the max of distribution?
  findmax=0;
# initialize
  n=0;

  if(min==""){
      min=1e50;
      fixmin=0;
  }
  if(max==""){
      max=-1e50;
      fixmax=0;
  }
# 
  sum=0.0;
}
{
# read in data
  if((substr($1,1,1)!="#")&&(NF>=col)&&(tolower($col)!="nan")){
    n++;
    x[n]=$(col);
    if(cols == -1)
      sum += 1.0;
    else{
      y[n]=$(cols);
      sum += y[n];
    }
    if(!fixmax)
	if(x[n]>max)
	    max=x[n];
    if(!fixmin)
	if(x[n]<min)
	    min=x[n];
  }
}
END{
  if(n == 0){
    print("error, no data read") > "/dev/stderr";
  }
  range=max-min;

  if(dx <= 0){# let number of bins decide spacing
      if((nbin=="")||(nbin<=0)){
	  print("error, dx is not set and nbin is",nbin) > "/dev/stderr";
      }
      nbin=int(nbin);# make sure nbin is an integer
      spacing=range/nbin;
  }else{ # have set bin width
      spacing=dx;
      nbin=int(range/spacing)+1;
  }
  if(spacing != 0.0){
# intialize bins
    for(i=1;i<=nbin;i++){
      bin_cnt[i]=0;
      bin_sum[i]=0;
    }
# count
    for(i=1;i<=n;i++){
      indexi=int((x[i]-min)/spacing)+1;# left sided boxes
      indexi=(indexi>nbin)?nbin:indexi;# account for max
      bin_cnt[indexi]++;
      if(cols == -1)
	bin_sum[indexi] += 1.0;
      else
	bin_sum[indexi] += y[i];
    }
    if(findmax){
# find bin with max entries and output center of bin	    
      hmax=0;
      for(i=1;i<=nbin;i++)
	if(bin_cnt[i]>hmax){
	  hmax=bin_cnt[i];
	  xmax=min+((i-1)+0.5)*spacing;
	}
      print(xmax);
    }else{
      if(0){# check if we actually got all entries in a bin
	nfound=0;ssum=0.0;
	for(i=1;i<=nbin;i++){
	  nfound += bin_cnt[i];
	  ssum += bin_sum[i];
	}
	if((nfound != n)||(ssum != sum)){
	  print("sorting error") > "/dev/stderr";
	  print("n",n,"nfound",nfound)> "/dev/stderr";
	  print("sum",sum,"ssum",ssum)> "/dev/stderr";
	}
      }
# output is left of bin, number of entries, fraction of 
# total number, and sum of y entries in bin
      for(i=1;i<=nbin;i++){
	if(midpoint){
	  print(min+((i-1)+.5)*spacing,
		bin_cnt[i],bin_cnt[i]/n,
		bin_sum[i],bin_sum[i]/sum);
	}else{
	  print(min+(i-1)*spacing,bin_cnt[i],bin_cnt[i]/n,bin_sum[i],bin_sum[i]/sum);
	  if(step_output){# repeat entr
	      print(min+i*spacing,bin_cnt[i],bin_cnt[i]/n,bin_sum[i],bin_sum[i]/sum);
	      print(min+i*spacing,0,0,0,0);
	      print(min+(i-1)*spacing,0,0,0);
	      print(min+(i-1)*spacing,bin_cnt[i],bin_cnt[i]/n,bin_sum[i],bin_sum[i]/sum);
	      print(">");
	  }
	}
      }
      if(step_output)# repeat last entry as zero y
	print(min+nbin*spacing,0,0,0,0);
    }
  }else{
# min and max are the same
    print("min and max of dataset are identical")> "/dev/stderr";
    print(min,max,n);
  }
}
