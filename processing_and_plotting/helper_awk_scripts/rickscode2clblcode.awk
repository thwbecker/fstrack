# convert rick's short plate codes to clb's long codes
BEGIN{
}
{
  if(substrc($1,1,1)!="#"){
    for(i=1;i<=NF;i++){
      if($i=="af")c="afr";
      else if($1=="an")c="ant";
      else if($1=="ar")c="arab";
      else if($1=="at")c="anat"; # anatolia, only in med brys
      else if($1=="au")c="aus";
      else if($1=="ca")c="car";
      else if($1=="co")c="coco";
      else if($1=="cr")c="cr";
      else if($1=="ea")c="";
      else if($1=="eu")c="eur";
      else if($1=="fa")c="fara";
      else if($1=="gr")c="gr";
      else if($1=="in")c="ind";
      else if($1=="iz")c="izan";
      else if($1=="jf")c="juan";
      else if($1=="ku")c="jula";
      else if($1=="lh")c="lh";
      else if($1=="na")c="nam";
      else if($1=="nz")c="naz";
      else if($1=="pa")c="pac";
      else if($1=="ph")c="phi";
      else if($1=="px")c="phen";
      else if($1=="sa")c="sam";
      else{ 
	print("error") > "/dev/stderr";
	c="";
      };
      printf("%s ",c);
    }
    printf("\n");
  }

}
