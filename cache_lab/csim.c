#include "cachelab.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

void print_usage() {
  printf("Usage: ./csim-ref [-hv] -s <s> -E <E> -b <b> -t <tracefile>\n");
  printf("  -h              Print this help message.\n");
  printf("  -v              Optional verbose flag.\n");
  printf("  -s <s>          Number of set index bits.\n");
  printf("  -E <E>          Number of lines per set.\n");
  printf("  -b <b>          Number of block offset bits.\n");
  printf("  -t <tracefile>  Name of the valgrind trace to replay.\n");
}

typedef struct {
  unsigned int valid_bit;
  unsigned int tag; 
  clock_t time_stamp;
}cache_line;

void cache_line_insert(cache_line *cacheLine, unsigned int tag){
  cacheLine->valid_bit = 1;
  cacheLine->tag = tag;
  cacheLine->time_stamp = clock();
}

int compare_times(const void* a, const void* b) {
  const cache_line *line1 = (const cache_line *)a;
  const cache_line *line2 = (const cache_line *)b;
  if (line1->time_stamp < line2->time_stamp) return -1;
  if (line1->time_stamp > line2->time_stamp) return 1;
  return 0;
}

int invoke(cache_line *cache, unsigned int v, unsigned int E, unsigned int currentTag, unsigned int currentSet){
  int i;
  //printf("Tag=%x Set=%d",currentTag,currentSet);
  for(i=0;i<E;i++){
    if(cache[currentSet*E+i].tag == currentTag&&cache[currentSet*E+i].valid_bit==1 ){
      cache_line_insert(&cache[currentSet*E+i],currentTag);
      if(v==1) printf(" hit\n");
      else if(v==3) printf(" hit");
      return 0;
    }
  }
  if(i==E){
    int j;
    for(j=0; j<E;j++){
      if(cache[currentSet*E+j].valid_bit == 0){
        cache_line_insert(&cache[currentSet*E+j],currentTag);
        if(v==1) printf(" miss\n");
        else if(v==3) printf(" miss");
        return 1;
      }
    }
    if(j==E){
      qsort(&cache[currentSet*E], E, sizeof(cache_line), compare_times); 
      cache_line_insert(&cache[currentSet*E],currentTag);
      if(v==1) printf(" miss eviction\n");
      else if(v==3) printf(" miss eviction");
      return 2;
    }
  }
  return -1;
}

int main(int argc, char* argv[]) {
  unsigned int opt, s = -1, E = -1, b = -1,v=0, byt;
  unsigned int addr;
  char *trace = NULL, ch;
  FILE * file_stream;

  int hitNum = 0, missNum = 0, evicNum = 0;

  // Process command-line arguments
  while ((opt = getopt(argc, argv, "hvs:E:b:t:")) != -1) {
    switch (opt) {
      case 's':
        s = atoi(optarg);
        break;
      case 'E':
        E = atoi(optarg);
        break;
      case 'b':
        b = atoi(optarg);
        break;
      case 't':
        trace = optarg;
        break;
      case 'v':
        v = 1;
        break;
      case 'h':
        print_usage();
        return 0;  // Exit after printing help
      default:
        print_usage();
        return 1;  // Exit with error
    }
  }
  int row = pow(2,s);
  cache_line *cache = (cache_line *)calloc(row*E, sizeof(cache_line));

  file_stream = fopen(trace, "r");

  while(fscanf(file_stream," %c %x,%d",&ch,&addr,&byt)==3){
    unsigned int currentTag = addr>>(s+b), currentSet = (addr - (currentTag<<(s+b)))>>b;
    if(ch!='I'&&v==1) printf("%c %x,%d", ch, addr, byt);
    switch(ch){
      case 'I':
        break;
      case 'M':
        switch(invoke((cache_line*)cache, v+2, E,currentTag, currentSet)){
          case 0:
            hitNum++;
            break;
          case 2:
            evicNum++;
          case 1:
            missNum++;
            break;

          default:
            break;
        }
      case 'L':
      case 'S':
        switch(invoke((cache_line*)cache, v, E, currentTag, currentSet)){
          case 0:
            hitNum++;
            break;
          case 2:
            evicNum++;
          case 1:
            missNum++;
            break;
          default:
            break;
        }
    }
  }

  printSummary(hitNum, missNum, evicNum);
  free(cache);
  fclose(file_stream);
  return 0;

}
