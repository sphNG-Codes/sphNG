/* grab - Grabs specified bits of a file and writes them to another
     output file */

#include "getline.c"

double time,timeff,mean_density,max_density,cen_density,total,ke,pot,internal,ang,lin,jean,utime,udens,alpha,beta,maxtot,mintot;
char line[80];

main()
 {  mintot = 1.0E+33;
    maxtot = -1.0E+33;
    while (getline(line,sizeof(line))>0) {
       if (sscanf(line, "             density    :%le time        :%le", &udens,&utime)==2)
          ;
       else if (sscanf(line, " TIME  :%le", &time)==1)
          printf("\n%16.9e ",time);
       else if (sscanf(line, " Free Fall Time  :%le",&timeff)==1)
          printf("%16.9e %16.9e ",timeff,time*utime);
       else if (sscanf(line, " total energy :%le",&total)==1) {
          printf("%e ",total);
	  if (total<mintot) mintot = total;
	  if (total>maxtot) maxtot = total;}
       else if (sscanf(line, " kinetic energy :%le",&ke)==1)
          printf("%e ",ke);
       else if (sscanf(line, " potential energy :%le",&pot)==1)
          printf("%e ",pot);
       else if (sscanf(line, " internal energy :%le",&internal)==1)
          printf("%e ",internal);
       else if (sscanf(line, " total angular momentum :%le",&ang)==1)
          printf("%e ",ang);
       else if (sscanf(line, " total linear momentum :%le",&lin)==1)
          printf("%e ",lin);
       else if (sscanf(line, " alpha :%le",&alpha)==1)
          printf("%e ",alpha);
       else if (sscanf(line, " (z)       beta parallel :%le",&beta)==1)
          printf("%e ",beta);
       else if (sscanf(line, " Jeans number :%le",&jean)==1)
          printf("%e ",jean);
       else if (sscanf(line, " density mean :%le max:%le cen:%le",
            &mean_density, &max_density, &cen_density)==3) {
          printf("%e %e ",mean_density,mean_density*udens);
          printf("%e %e ",max_density,max_density*udens);
          printf("%e %e ",cen_density,cen_density*udens);}
       else if (sscanf(line, " density mean :%le max:%le",
            &mean_density, &max_density)==2) {
          printf("%e %e ",mean_density,mean_density*udens);
          printf("%e %e ",max_density,max_density*udens);};
     };
    printf("\n %e \n",mintot-maxtot);
 };
