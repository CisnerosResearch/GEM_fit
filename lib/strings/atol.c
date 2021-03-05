#include <stdio.h>
int main(int argc, char *argv[]){
    fprintf(stderr,"I=%d\n",strtol(argv[1],NULL,0));
    return 0;
}
