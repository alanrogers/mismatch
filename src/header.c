#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "memstack.h"
#include "version.h"
#include "mismatch.h"
#include "header.h"
void header(char *program, char *description, FILE *fp)
{
  char buff[80], word[20];
  int i;
  
  centerline(program,fp);
  centerline(description,fp);
  centerline("by Alan R. Rogers",fp);
  sprintf(buff, "Version %s", VERSION);
  centerline(buff,fp);
  centerline(DATE, fp);
  for(i=0; program[i] != '\0'; i++)
    word[i] = tolower(program[i]);
  word[i] = '\0';
  sprintf(buff,"Type `%s -- ' for help", word);
  centerline(buff, fp);
  fflush(fp);
}

/** center a character string on output line */
void centerline(char *s, FILE *fp)
{
  int pad;

  pad = 80 - strlen(s);
  if(pad < 0)
    pad = 0;
  else
    pad /= 2;
  putc(START_COMMENT, fp);   /* begin line with comment character */
  pad--;
  while(--pad > 0)
    putc(' ', fp);
  fputs(s, fp);
  putc('\n', fp);
}
