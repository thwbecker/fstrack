  H  @   k820309    °          8.1     jEäB                                                                                                                        
       decmod.f90 DECMOD                                         
                                               
       #     @                                           #JACOBI%SIZE    #JACOBI%SUM    #JACOBI%ABS    #JACOBI%MERGE    #JACOBI%SQRT    #A 	   #D 
   #V    #NROT                                         SIZE                                      SUM                                      ABS                                      MERGE                                      SQRT                               	           
           &           &                          0                          
           
           &                                                                 
           &           &                                                              #     @                                           #EIGSRT%SIZE    #D    #V                                         SIZE                                          
 (          &                                                                 
 )          &           &                       (     `                                                
#V21KEL%SQRT    #KEL    p      p        p                                                   SQRT                                     $      
     p      p      p        p      p              (     `                           $                     
#KELC%SQRT    #C    p      p      p        p      p                                                   SQRT                                     $      
     p      p      p        p      p              (     `                           $                      
#T4    p      p      p        p      p                                                  Q      
 
    p (     p      p      p      p        p      p      p      p              (     `                           Q                      
#T4    #R    p (     p      p      p      p        p      p      p      p                                                  Q      
     p (     p      p      p      p        p      p      p      p                                                  	      
     p      p      p        p      p              (     `                           Q                      
#C    p (     p      p      p      p        p      p      p      p                                                  $      
     p      p      p        p      p                     @                              
          @                              
          @                               
          @                         !     
          @                         "     	      
  p      p      p        p      p                     @                         #     	      
  p      p      p        p      p                     @ @                       $     	      
  p      p      p        p      p                     @ @                       %     	      
  p      p      p        p      p                     @ @                       &     	      
  p      p      p        p      p                     @ @                       '           
  p      p        p              #     @                          (                   #C )                                  )     $      
     p      p      p        p      p              #     @                          *                  #DEC_ISO%DOT_PRODUCT +   #DEC_ISO%SQRT ,   #C -                                  +     DOT_PRODUCT                                ,     SQRT       D @                       -     $      
     p      p      p        p      p              #     @                          .                  #DEC_SCCA%CSHIFT /   #DEC_SCCA%TRANSPOSE 0   #DEC_SCCA%MINLOC 1   #DEC_SCCA%ACOS 2   #DEC_SCCA%SIGN 3   #DEC_SCCA%ABS 4   #DEC_SCCA%DOT_PRODUCT 5   #DEC_SCCA%SQRT 6   #C 7                                  /     CSHIFT                                0     TRANSPOSE                                1     MINLOC                                2     ACOS                                3     SIGN                                4     ABS                                5     DOT_PRODUCT                                6     SQRT       D @                       7     $      
     p      p      p        p      p              #     @                         8                  #DEC_PROJ%DOT_PRODUCT 9   #DEC_PROJ%SQRT :   #X ;   #XH <   #XD =   #DEV >   #NSYM ?                                  9     DOT_PRODUCT                                :     SQRT        @                       ;           
     p      p        p                    D                         <           
     p      p        p                    D @                       =           
     p      p        p                    D                         >     
                                  ?           ¨         fn#fn    Â   4   J   NRMOD    ö   4   J   TENSMOD    *  ¯       JACOBI+NRMOD "   Ù  1      JACOBI%SIZE+NRMOD !   
  0      JACOBI%SUM+NRMOD !   :  0      JACOBI%ABS+NRMOD #   j  2      JACOBI%MERGE+NRMOD "     1      JACOBI%SQRT+NRMOD    Í  p   a   JACOBI%A+NRMOD    =  `   a   JACOBI%D+NRMOD      p   a   JACOBI%V+NRMOD "     0   a   JACOBI%NROT+NRMOD    =  [       EIGSRT+NRMOD "     1      EIGSRT%SIZE+NRMOD    É  `   a   EIGSRT%D+NRMOD    )  p   a   EIGSRT%V+NRMOD             V21KEL+TENSMOD $   +  1      V21KEL%SQRT+TENSMOD #   \     a   V21KEL%KEL+TENSMOD    à  ¦       KELC+TENSMOD "     1      KELC%SQRT+TENSMOD    ·     a   KELC%C+TENSMOD    ;         CT4+TENSMOD    Ó  ´   a   CT4%T4+TENSMOD    	  Ï       RT4+TENSMOD    V
  ´   a   RT4%T4+TENSMOD    
     a   RT4%R+TENSMOD      Ç       T4C+TENSMOD    U     a   T4C%C+TENSMOD    Ù  0       KISO    	  0       GISO    9  0       PERCA    i  0       PERCT             DI             VO             VDI             VVO             SCC      h       X21      C       DEC_CONT    Ä     a   DEC_CONT%C    H  n       DEC_ISO $   ¶  8      DEC_ISO%DOT_PRODUCT    î  1      DEC_ISO%SQRT         a   DEC_ISO%C    £  ê       DEC_SCCA       3      DEC_SCCA%CSHIFT #   À  6      DEC_SCCA%TRANSPOSE     ö  3      DEC_SCCA%MINLOC    )  1      DEC_SCCA%ACOS    Z  1      DEC_SCCA%SIGN      0      DEC_SCCA%ABS %   »  8      DEC_SCCA%DOT_PRODUCT    ó  1      DEC_SCCA%SQRT    $     a   DEC_SCCA%C    ¨         DEC_PROJ %   ;  8      DEC_PROJ%DOT_PRODUCT    s  1      DEC_PROJ%SQRT    ¤  l   a   DEC_PROJ%X      l   a   DEC_PROJ%XH    |  l   a   DEC_PROJ%XD    è  0   a   DEC_PROJ%DEV      0   a   DEC_PROJ%NSYM 