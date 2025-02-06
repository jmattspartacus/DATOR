/* 
Tim Gray - tgray30@utk.edu - 01/2025
takes a *.ldf file containing presumed ORRUBA data. Converts into GEB format (type-19 events).

compile with 
gcc -o LDFConvert LDFConvert.c -O3 -std=c99

*/

#define _POSIX_C_SOURCE 1 //this is for the fileno() function

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <zlib.h>

//these are for the (LDF) channels in which the MyRIAD timestamp is stored
#define LDF_TS_LSB 1000  //lowest significant bit
#define LDF_TS_MSB 1001  //middle significant bit
#define LDF_TS_HSB 1002  //highest significant bit

#define MAX_INT64 18446744073709551615ULL

#define DATA 0x41544144
#define SCALER 0x4C414353

struct ldfdata {
    char type[4];
    unsigned int num_words;
    unsigned int data[8192];
};

struct GEBHeader {
    int type;
    int length; /* length of payload following the header, in bytes */
    unsigned long long int timestamp;
};

struct cFile { //handler for compression
  FILE *ptr;
  gzFile gzf; //compressed file  
};

bool compressed = false; //yes, this is global

int ReadLDF(FILE *fptr, struct ldfdata *data, int *startWord, int *stopWord) {
  int readsuccess = 1;

  unsigned int typeID = *((unsigned int*)(&data->type[0])); //crucify me for crimes of memory abuse
  
  if ( typeID == DATA ){ //lit. 44 41 54 41 packed as 4144 4154 in hexdump
    int good = 0;
    int wordPtr = *startWord;
    for (int word = wordPtr; word < 8192; ++word) {
      if (good == 0 && data->data[word] != 0xFFFFFFFF) { good = 1; startWord[0] = word; continue; }
      if (good == 1 && data->data[word] == 0xFFFFFFFF) { stopWord[0] = word; break; }

    }
    if (good == 1) { return 1; }
  }
  else if ( typeID == SCALER ) {
    if (startWord[0] == 0) { 
      startWord[0] = 0;
      stopWord[0] = 1023;
      data->num_words=1024;
      return 1;
    }
    /*
    if (startWord[0] == 0) {
      printf("\n");
      startWord[0] = 0;
      int charCt = 0;
      int spaceCt = 0;
      for (int word = startWord[0]; word < 8192; ++word) {
        char string[4];
        memcpy(&string[0], &data->data[word], 1);

        printf("%#x ", data->data[word]);
        for (int i=0; i<4; ++i) {
          charCt += 1;            
          if (charCt%80 == 0) {printf("\n");}
          //if (string[i] == 0 || string[i] == 127) { continue; }
          if (string[i] == 32) { spaceCt += 1; }
          else { spaceCt = 0; }
          
          printf(" %u:%c ", string[i], string[i]);
            
          //if (string[i]==10) { printf("\n"); } 
          //else {printf("%u  %c\n", string[i], string[i]);}
        }

        if (spaceCt == 20) { stopWord[0] = word; break; };
        
        if (data->data[word] == 0) { stopWord[0] = word; break; }
      }
      printf("\n");
      return 1;
    }
    */
  }
  
  while ( true ) {
    readsuccess = fread ( data, sizeof ( *data ), 1, fptr );
    if ( readsuccess == 0 ) { break; }
    typeID = *((unsigned int*)(&data->type[0])); //crucify me for crimes of memory abuse
    if ( typeID == DATA || typeID == SCALER ){ //lit. 44 41 54 41 packed as 4144 4154 in hexdump
        break;
    }
  }
  if (readsuccess == 0) { return 0; }
  startWord[0] = 0;
  stopWord[0] = 0;
  return ReadLDF(fptr, data, startWord, stopWord);
}

int GetLDFTime(struct ldfdata *data, int *startWord, int *stopWord, unsigned long long int *ldf_time) {
  *ldf_time = 0;
  for ( int word = *startWord; word <= *stopWord; ++word) {
    unsigned short int channel = (data->data[word] & 0x7FFF);
    unsigned short int value = (data->data[word] >> 16 );
  
    if (data->data[word] == 0xFFFFFFFF) {
      data->num_words = *stopWord - *startWord + 1;
      break;
    }
  }
  
  return 1;      
}


void WriteLDF(struct cFile *outfile, struct ldfdata *data, int *startWord, int *stopWord, unsigned long long int time_LDF) {
  struct GEBHeader header;
  header.timestamp = time_LDF;
  unsigned int typeID = *((unsigned int*)(&data->type[0])); //crucify me for crimes of memory abuse
  if (typeID == DATA) {
    header.type = 19;
  }
  else if (typeID == SCALER) {
    header.type = 22;
  }
  header.length = data->num_words*4;

  if (!compressed) {
    fwrite(&header, sizeof(header), 1, outfile->ptr);
    fwrite(&data->data[*startWord], header.length, 1, outfile->ptr);
  }
  else {
    gzwrite(outfile->gzf, &header, sizeof(header));
    gzwrite(outfile->gzf, &data->data[*startWord], header.length);
  }      

  //printf("writing LDF, header.length = %i, time = %5.1f\n", header.length, (double)header.timestamp/1e8/60.0);
  return;
}

void WriteGEB(struct cFile *outfile, struct cFile *file, struct GEBHeader *header) {
  unsigned int sub[8192];
  int retval;
  if (!compressed) {
    retval = fread(sub, header->length, 1, file->ptr);
    fwrite(header, sizeof(*header), 1, outfile->ptr);
    fwrite(&sub[0], header->length, 1, outfile->ptr);
  }
  else {
    retval = gzread(file->gzf, sub, header->length);
    gzwrite(outfile->gzf, header, sizeof(*header));
    gzwrite(outfile->gzf, &sub[0], header->length);
  }
  //printf("writing GEB, header.length = %i, time = %5.1f\n", header->length, (double)header->timestamp/1e8/60.0);
  return;
}

int main(int argc, const char **argv) {
  if (argc < 2) {
    printf("Useage: LDFConvert ORRUBA.ldf\n");
    exit(1);
  }

  struct cFile gebfile;
  struct cFile outfile;
 
  FILE *ldffile = fopen(argv[1], "r");
  if (!compressed) { outfile.ptr = fopen("GlobalConverted.dat", "w"); }
  else { outfile.ptr = fopen("GlobalConverted.dat.gz", "w"); outfile.gzf = gzdopen(fileno(outfile.ptr), "w"); }

  int readLDF = 1;
  struct ldfdata dataLDF;
  int startWord = 0;
  int stopWord = 0;
  struct GEBHeader header;
  unsigned long long int time_LDF = 0;

  unsigned long long int readCtrLDF = 0;
  unsigned long long int ctrLDF = 0;
  unsigned long long discards = 0;

  unsigned long long int type19 = 0;
  unsigned long long int type22 = 0;
  
  while (true) {
    int retval = 1;
    if ((ctrLDF) % 1000 == 0) {
      printf("%llu LDF events written\r", ctrLDF);
    }
    
    retval = 0;
    unsigned int typeID;
    while (retval == 0) {
      retval = ReadLDF(ldffile, &dataLDF, &startWord, &stopWord);
      typeID = *((unsigned int*)(&dataLDF.type[0])); //crucify me for crimes of memory abuse
      if (retval == 0) {
        readLDF = 0;
        time_LDF = MAX_INT64; //max 64 bit integer to ensure we never try to write GEB again
        retval = 1;
      }
      else {
        readCtrLDF += 1;
        readLDF = 0;

        retval = GetLDFTime(&dataLDF, &startWord, &stopWord, &time_LDF);
        time_LDF = 1000*(ctrLDF+1);  //no global timestamp for LDF format, so make the GEB timestamp 10 us * event No.
        
        if (retval == 0) {
          printf("\nSpurious MyRIAD, writing immediately\n");
          if (typeID == DATA) {
            type19 += 1;
          }
          else if (typeID == SCALER) {
            type22 += 1;
          }                      
          WriteLDF(&outfile, &dataLDF, &startWord, &stopWord, time_LDF);
          readLDF = 1;
          startWord = stopWord + 1;
          ctrLDF += 1;
        }
          
        //printf("%i  %llu\n", retval, time_LDF);
      }
    }

    if (time_LDF == MAX_INT64) {
      break;
    }

    WriteLDF(&outfile, &dataLDF, &startWord, &stopWord, time_LDF);
    if (typeID == DATA) {
      type19 += 1;
    }
    else if (typeID == SCALER) {
      type22 += 1;
    }            

    readLDF = 1;
    startWord = stopWord + 1;
    ctrLDF += 1;
    
  }

  fclose(ldffile);
  fclose(outfile.ptr);
  printf("\nLDF events: %llu read, %llu written", readCtrLDF, ctrLDF);
  printf("\nType 19: %llu, Type 22: %llu\n", type19, type22);
  printf("Done\n");
  
}
