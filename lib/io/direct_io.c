#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <fcntl.h>
#include <string.h>
#include <zlib.h>

static char *cstring(char *string, size_t len){
  char *result = malloc(len+1);
  strncpy(result,string,len);
  result[len]=0;
  return result;
}

size_t fd_open_() __attribute__ ((alias("fd_open__")));
size_t fd_open__(char *filename, int *mode, int filename_len){
  int fd;
  char *cfilename=cstring(filename,filename_len);
  fd=open(cfilename,*mode);
  free(cfilename);
  return (size_t)fd;
}

size_t fd_read_() __attribute__ ((alias("fd_read__")));
size_t fd_read__(int *fd, void *buffer, size_t *len){
  return read(*fd, buffer, *len);
}

size_t fd_write_() __attribute__ ((alias("fd_write__")));
size_t fd_write__(int *fd, void *buffer, size_t *len){
  return write(*fd, buffer, *len);
}

size_t fd_close_() __attribute__ ((alias("fd_close__")));
size_t fd_close__(int *fd){
  return (size_t) close(*fd);
}

size_t stream_open_() __attribute__ ((alias("stream_open__")));
size_t stream_open__(char *filename, char *mode, int filename_len, int mode_len){
  FILE *pFile;
  char *cfilename=cstring(filename,filename_len);
  char *cmode=cstring(mode,mode_len);
  pFile=fopen(cfilename,cmode);
  //fprintf(stderr,"stream_open: fopen(cfilename=\"%s\",cmode\"%s\")=%p;\n",cfilename,cmode,pFile);
  free(cmode);
  free(cfilename);
  return (size_t)pFile;
}

size_t stream_read_() __attribute__ ((alias("stream_read__")));
size_t stream_read__(FILE **ppFile, void *buffer, size_t *len){
  return fread(buffer, 1, *len, *ppFile);
}

size_t stream_write_() __attribute__ ((alias("stream_write__")));
size_t stream_write__(FILE **ppFile, void *buffer, size_t *len){
  return fwrite(buffer, 1, *len, *ppFile);
}

int stream_gets_() __attribute__ ((alias("stream_gets__")));
int stream_gets__(FILE **ppFile, char *buffer, int buffer_len){
  if (fgets(buffer,buffer_len,*ppFile)) return 1;
  else return 0;
}

size_t stream_ftell_() __attribute__ ((alias("stream_ftell__")));
size_t stream_ftell__(FILE **ppFile){
  return (size_t)ftell(*ppFile);
}

size_t stream_close_() __attribute__ ((alias("stream_close__")));
size_t stream_close__(FILE **ppFile){
  return (size_t) fclose(*ppFile);
}

size_t gz_reopen_stream_() __attribute__ ((alias("gz_reopen_stream__")));
size_t gz_reopen_stream__(FILE **ppFile, char *mode, int mode_len){
  gzFile File;
  char *cmode=cstring(mode,mode_len);
  //fprintf(stderr,"gz_reopen_stream(*ppFile=%p,cmode=\"%s\");\n",*ppFile,cmode);
  File=gzdopen(dup(fileno(*ppFile)),cmode);
  fclose(*ppFile);
  free(cmode);
  return (size_t)File;
}

size_t gz_open_() __attribute__ ((alias("gz_open__")));
size_t gz_open__(char *filename, char *mode, int filename_len, int mode_len){
  gzFile File;
  char *cfilename=cstring(filename,filename_len);
  char *cmode=cstring(mode,mode_len);
  File=gzopen(cfilename,cmode);
  //fprintf(stderr,"gzopen(cfilename=\"%s\",cmode\"%s\")=%p;\n",cfilename,cmode,File);
  free(cmode);
  free(cfilename);
  return (size_t)File;
}

size_t gz_read_() __attribute__ ((alias("gz_read__")));
size_t gz_read__(gzFile *pFile, void *buffer, size_t *len){
  size_t result= gzread(*pFile, buffer, *len);
  //fprintf(stderr,"gzread(gzFile=%p, buffer=%p, len=%d)=%d;\n", *pFile, buffer, *len, result);
  return result;
}

size_t gz_write_() __attribute__ ((alias("gz_write__")));
size_t gz_write__(gzFile *pFile, void *buffer, size_t *len){
  size_t result;
  result = gzwrite(*pFile, buffer, *len);
  //fprintf(stderr,"gzwrite(gzFile=%p, buffer=%p, len=%d)=%d;\n", *pFile, buffer, *len, result);
  return result;
}

int gz_gets_() __attribute__ ((alias("gz_gets__")));
int gz_gets__(gzFile *pFile, char *buffer, int buffer_len){
  if (gzgets(*pFile,buffer,buffer_len)) return 1;
  else return 0;
}

size_t gz_ftell_() __attribute__ ((alias("gz_ftell__")));
size_t gz_ftell__(gzFile *pFile){
  z_off_t result = gztell(*pFile);
  //fprintf(stderr,"gztell(%p)=%d, %d;\n",pFile,pos,result);
  return (size_t)result;
}

size_t gz_close_() __attribute__ ((alias("gz_close__")));
size_t gz_close__(gzFile *pFile){
  return (size_t) gzclose(*pFile);
}

