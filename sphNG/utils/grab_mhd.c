/* grab - Grabs specified bits of a file and writes them to another
     output file */

//#include "getline.c"

#include <stdio.h>

double
time,timeff,mean_density,max_density,cen_density,total,ke,pot,
internal,mag,ang,jean,utime,udens,uang,umagfd,alpha,beta,maxtot,
mintot,lin,divBmax,divBmean,Jmax,Jmean,divBhmax,max_beta,min_beta,
mean_beta,fluxtot,crosshel,min_b,mean_b,max_b;
//char line[120];
char *line = NULL;
size_t len = 0;

main()
 {  mintot = 1.0E+33;
    maxtot = -1.0E+33;
//    while (getline(line,sizeof(line))>0) {
    while (getline(&line, &len, stdin) > 0) {
       if (sscanf(line, "             density    :%le time        :%le", &udens,&utime)==2)
          ;
       else if (sscanf(line, " ang. mom. :%le mag field :%le", &uang, &umagfd)==2)
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
       else if (sscanf(line, " magnetic energy :%le",&mag)==1)
          printf("%e ",mag);	  
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
       else if (sscanf(line, "  plasma beta  min :%le max :%le mean :%le",
            &min_beta, &max_beta, &mean_beta)==3) {
          printf("%e ",min_beta);
          printf("%e ",max_beta);
          printf("%e ",mean_beta);}
       else if (sscanf(line, "  div B max : %le mean: %le", &divBmax, &divBmean)==2) {
          printf("%e ",divBmax);
          printf("%e ",divBmean);}          
       else if (sscanf(line, "  curl B max : %le mean: %le", &Jmax, &Jmean)==2) {
          printf("%e ",Jmax);
          printf("%e ",Jmean);}
       else if (sscanf(line, "  divB*h/B max : %le", &divBhmax)==1)
          printf("%e ",divBhmax);
       else if (sscanf(line, " total magnetic flux  (int div B dV)   :%le",&fluxtot)==1)
          printf("%e ",fluxtot);
       else if (sscanf(line, " total cross helicity (int v.B dV) :%le",&crosshel)==1)
          printf("%e ",crosshel);  	  
       else if (sscanf(line, " density mean :%le max:%le cen:%le",
            &mean_density, &max_density, &cen_density)==3) {
          printf("%e %e ",mean_density,mean_density*udens);
          printf("%e %e ",max_density,max_density*udens);
          printf("%e %e ",cen_density,cen_density*udens);}
       else if (sscanf(line, " mag field   min :%le mean :%le max: %le", &min_b, &mean_b, &max_b)==3) {
          printf("%e %e ",min_b,min_b*umagfd);
          printf("%e %e ",mean_b,mean_b*umagfd);
          printf("%e %e ",max_b,max_b*umagfd);};
     };
    printf("\n %e \n",mintot-maxtot);
 };
