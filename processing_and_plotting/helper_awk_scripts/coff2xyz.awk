#
# convert geomview (C)off format to xyz for nodes only
#
BEGIN{
  readi=0;
  nread=0;

}
{
  if(match($0,"OFF"))
    readi=1;
  if(readi){
    if(($1 !="")&&(NF==3)){
      nrnd=$1;nrel=$2;
      readi=0;
    }
  }else{
    if(nrnd)
      if(nread <= nrnd){
	if(($1 !="")&&(NF==3)){
	  nread++;
	  x[nread]=$1;y[nread]=$2;z[nread]=$3;
	}
      }
  }
}
END{
# output of nodal coordinates
  for(i=1;i<=nrnd;i++){
    print(x[i],y[i],z[i]);
  }
}
