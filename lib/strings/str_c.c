#include <stdlib.h>
#include <string.h>

#define MAXLEN 128

double atod_ (char *str, int len) {
  char buffer[MAXLEN+1];
  int buflen;
  if (len >= buflen)
    buflen = MAXLEN;
  else
    buflen = len;
  strncpy (buffer, str, buflen);
  buffer[buflen] = 0;
  return atof (buffer);
}

int atoi_ (char *str, int len) {
  char buffer[MAXLEN+1];
  int buflen;
  if (len >= buflen)
    buflen = MAXLEN;
  else
    buflen = len;
  strncpy (buffer, str, buflen);
  buffer[buflen] = 0;
  return atoi (buffer);
}

long int atol_ (char *str, int len) {
  char buffer[MAXLEN+1];
  int buflen;
  if (len >= buflen)
    buflen = MAXLEN;
  else
    buflen = len;
  strncpy (buffer, str, buflen);
  buffer[buflen] = 0;
  return atol (buffer);
}

int strtol_ (char *str, int len) {
  char buffer[MAXLEN+1];
  int buflen;
  if (len >= buflen)
    buflen = MAXLEN;
  else
    buflen = len;
  strncpy (buffer, str, buflen);
  buffer[buflen] = 0;
  return strtol(buffer,NULL,0);
}
