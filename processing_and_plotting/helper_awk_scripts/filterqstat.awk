BEGIN{

  n=0;
}{
  if($1 == "Job_Name"){name=$3;n++;}
  if($1 == "queue"){queue=$3;n++;}
  if($1 == "Job_Owner"){split($3,a,"@");own=a[1];n++;}
  if($1 == "resources_used.walltime"){time=$3;n++;}
  if($1 == "job_state"){state=$3;n++;}

  if(n==5){printf("%25s %10s %10s %2s %10s\n",name,own,queue,state,time);n=0;}

}
END{

  if(n==5){printf("%25s %10s %10s %2s %10s\n",name,own,queue,state,time);n=0;}
}
