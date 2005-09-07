/* getline - Reads a line from std input into, s with max length, limit */
 
#include <stdio.h>
 
int getline(line,limit)
char line[];
int limit;
{  int c,i;
 
   for (i=0; i<limit-1 && (c=getchar())!=EOF && c!='\n'; ++i)
      line[i]=c;
 
   if (c=='\n') {
      line[i]=c;
      ++i;
   }
   line[i]='\0';
   return i;
};
  
