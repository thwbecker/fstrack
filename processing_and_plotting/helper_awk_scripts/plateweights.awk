BEGIN{
  # weights in % based on 0.25 degree grids
  w[1]=11.7344; p[1]="ANT";
  w[2]=9.5846; p[2]="AUS";
  w[3]=15.2668; p[3]="AFR";
  w[4]=21.1064; p[4]="PAC";
  w[5]=13.3349; p[5]="EUR";
  w[6]=11.2891; p[6]="NAM";
  w[7]=3.2324; p[7]="NAZ";
  w[8]=0.6048; p[8]="COC";
  w[9]=0.7286; p[9]="CAR";
  w[10]=0.9961; p[10]="ARA";
  w[11]=1.1252; p[11]="PHI";
  w[12]=8.5724; p[12]="SAM";
  w[13]=2.3693; p[13]="IND";
  w[14]=0.05449; p[14]="JDF";
  n=14;
  maxnf=1;
}
{
  if(substr($1,1,1)!="#"){
    for(i=1;i<=n;i++){
      if($1==p[i]){
	if(NF-1 > maxnf)
	  maxnf=NF-1;
	for(j=2;j<=NF;j++){
	  if(tolower($j) != "nan"){
	    ws[j-1] += w[i];
	    s[j-1]  += w[i]*$j;
	  }
	}
      }
    }
  }
}
END{
  for(i=1;i<=maxnf;i++)
    if(ws[i] != 0)
      printf("%g ",s[i]/ws[i]);
  printf("\n");


}
