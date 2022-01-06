BEGIN{
}
{
  if($1!=""){
    n=split($0,a,",");
    for(i=1;i<=n;i++){
      while(substr(a[i],1,1)==" ")
	a[i]=substr(a[i],2);
      while(substr(a[i],length(a[i]),1)==" ")
	a[n]=substr(a[i],1,length(a[i]-1));
    }
    printf("@Article{,\n");
    printf("  author = 	 {");
    for(i=1;i<=n-5;i++)
      if(i==1)
	printf("%s, ",a[i]);
      else
	printf("%s ",a[i]);
    printf("},\n");
    printf("  title = 	 {");
    printf("%s ",a[n-4]);
    printf("},\n");
    
    printf("  journal = 	 {%s},\n",a[n-3]);
    printf("  year = 	 {%s},\n",substr(a[n],1,4));
    printf("  OPTvolume = 	 {%s},\n",a[n-2]);
    gsub("-","--",a[n-1]);
    printf("  OPTpages = 	 {%s}\n",a[n-1]);
    printf("}\n\n\n");
  }
}
