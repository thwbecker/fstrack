BEGIN{
}
{
  if($1!="" && (substr($1,1,1)!="#")){
    if($2==""){
      print("expecting lon lat data format, lat is empty") > "/dev/stderr"
	}
    weight=cos($2*0.017453293);
    for(i=3;i<=NF;i++){
      nv=i-2;
      sum[nv] += weight*$(i);
    }
    sum_weight += weight;
  }
}
END{
  if(sum_weight != 0.0){
    for(i=1;i<=nv;i++)
      printf("%20.16e ",sum[i]/sum_weight);
    printf("\n");
  }
}

