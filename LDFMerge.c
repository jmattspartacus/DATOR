/* 
Tim Gray - tgray30@utk.edu - 04/2024
takes a Global.dat (or HFC.dat) and a *.ldf file containing presumed ORRUBA data.
Both of these should be time ordered. Merges the two together, ignoring any type-19
events in the Global.dat, and making new type-19 events from the *.ldf file

compile with 
gcc -o LDFMerge LDFMerge.c -O3 -std=c99

Updates:
05/2024 - Correct handling of "Spurious MyRIAD" events where each 16-bit word of the 48-bit timestamp is the same: 
          writing these out immediately
06/2024 - Support for compressed files using zlib
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
  unsigned short int lsb = 0;
  unsigned short int msb = 0;
  unsigned short int hsb = 0;
  for ( int word = *startWord; word <= *stopWord; ++word) {
    unsigned short int channel = (data->data[word] & 0x7FFF);
    unsigned short int value = (data->data[word] >> 16 );
  
    if (channel == LDF_TS_LSB) {
      lsb = value;
    }
    else if (channel == LDF_TS_MSB) {
      msb = value;
    }
    else if (channel == LDF_TS_HSB) {
      hsb = value;
    }

    if (data->data[word] == 0xFFFFFFFF) {
      data->num_words = *stopWord - *startWord + 1;
      break;
    }
  }

  if (lsb == msb && msb == hsb) {
    //printf("%i %i %i \n", hsb, msb, lsb);
    return 0;
  }
  else {
    *ldf_time += ((unsigned long long)lsb);
    *ldf_time += ((unsigned long long)msb<<16);
    *ldf_time += ((unsigned long long)hsb<<32);
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

int ReadGEB(struct cFile *file, struct GEBHeader *header, unsigned long long *discards) {
  int retval;
  unsigned int sub[8192];
  while (true) {
    if (!compressed) {
      retval = fread(header, sizeof(*header), 1, file->ptr);
    }
    else {
      retval = gzread(file->gzf, header, sizeof(*header));
    }
      
    if (retval == 0) { break; }

    if (header->type == 19) {
      if (!compressed) {
        retval = fread(sub, header->length, 1, file->ptr);
      }
      else {
        retval = gzread(file->gzf, sub, header->length);
      }
      if (retval == 0) { break; }
      //printf("type 19; length = %i\n", header->length);
      discards[0] += 1;
      continue;
    } //skip type-19 events
    else {
      break;
    }
  }
  return retval;
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
    printf("Useage: LDFMerge Global.dat[.gz] ORRUBA.ldf [outputPath.dat[.gz]]\n");
    exit(1);
  }

  char *compressed_ext = ".dat.gz";
  char *raw_ext = ".dat";
  bool compressed = false;

  if (!strcmp(&argv[1][strlen(argv[1])-7], compressed_ext)) {
    compressed = true;
  }
  else if (!strcmp(&argv[1][strlen(argv[1])-4], raw_ext)) {
    compressed = false;
  }
  else {
    printf("Unrecognized file extension for %s\n", argv[1]);
    exit(1);
  }

  struct cFile gebfile;
  struct cFile outfile;
 
  gebfile.ptr = fopen(argv[1], "r");
  if (compressed) { gebfile.gzf = gzdopen(fileno(gebfile.ptr), "r"); }
  FILE *ldffile = fopen(argv[2], "r");
  // this next block allows selecting an output path
  char *outf;
  int compressed_out = false;
  if(argc > 3){
    outf=argv[3];
    if (!strcmp(&outf[strlen(outf)-7], compressed_ext)) {
      compressed_out = true;
    }
    else if (!strcmp(&outf[strlen(outf)-4], raw_ext)) {
      compressed_out = false;
    }
    else {
      printf("Unrecognized file extension for %s\n", outf);
      exit(1);
    }
  } else {
    if (compressed){
      outf = "GlobalMerge.dat.gz";
    } else {
      outf = "GlobalMerge.dat";
    }
  }
  if (!compressed_out) { outfile.ptr = fopen(outf, "w"); }
  else { outfile.ptr = fopen(outf, "w"); outfile.gzf = gzdopen(fileno(outfile.ptr), "w"); }

  int readGEB = 1;
  int readLDF = 1;
  struct ldfdata dataLDF;
  int startWord = 0;
  int stopWord = 0;
  struct GEBHeader header;
  unsigned long long int time_LDF = 0;
  unsigned long long int last_time_LDF = 0;
  unsigned long long int time_GEB = 0;

  unsigned long long int readCtrLDF = 0;
  unsigned long long int readCtrGEB = 0;
  unsigned long long int ctrLDF = 0;
  unsigned long long int ctrGEB = 0;
  unsigned long long discards = 0;

  unsigned long long int type19 = 0;
  unsigned long long int type19_badmyriad_a = 0;
  unsigned long long int type19_badmyriad_b = 0;
  unsigned long long int type22 = 0;
  
  while (true) {
    int retval = 1;
    if ((ctrLDF + ctrGEB) % 1000 == 0) {
      printf("%llu LDF, %llu GEB events written\r", ctrLDF, ctrGEB);
    }
    if (readGEB) {
      retval = ReadGEB(&gebfile, &header, &discards);
      if (retval == 0) { // end of GEB file
        readGEB = 0;
        time_GEB = MAX_INT64; //max 64 bit integer to ensure we never try to write GEB again
      }
      else {
        readCtrGEB += 1;
        readGEB = 0;
        time_GEB = header.timestamp;
      }
    }
    if (readLDF) {
      retval = 0;
      while (retval == 0) {
        retval = ReadLDF(ldffile, &dataLDF, &startWord, &stopWord);
        if (retval == 0) {
          readLDF = 0;
          time_LDF = MAX_INT64; //max 64 bit integer to ensure we never try to write GEB again
          retval = 1;
        }
        else {
          readCtrLDF += 1;
          readLDF = 0;

          unsigned int typeID = *((unsigned int*)(&dataLDF.type[0])); //crucify me for crimes of memory abuse

          if (typeID == SCALER) { time_LDF = time_LDF < time_GEB ? time_LDF : time_GEB; }
          else {
            retval = GetLDFTime(&dataLDF, &startWord, &stopWord, &time_LDF);
            long long int delta = (long long int)time_LDF - (long long int)last_time_LDF;
            bool immediate_write = false;
            if (retval == 0) {
              //printf("\nBad MyRIAD (a), writing immediately\n");
              type19_badmyriad_a += 1;
              immediate_write = true;
            }
            else if (delta < 0) {
              printf("\nWarning! LDF file appears to be improperly time ordered!\n");
              immediate_write = true;
            }
            else if (delta > 1e9 && type19>0) {
              //printf("\nBad MyRIAD (b), writing immediately\n");
              type19_badmyriad_b += 1;
              immediate_write = true;
            }

            if (immediate_write) {
              if (typeID == DATA) {
                type19 += 1;
              }
               if (typeID == SCALER) {
                type22 += 1;
              }                      
              WriteLDF(&outfile, &dataLDF, &startWord, &stopWord, time_LDF);
              readLDF = 1;
              startWord = stopWord + 1;
              ctrLDF += 1;
            }
            else {
              last_time_LDF = time_LDF;
            }
          
            //printf("%i  %llu\n", retval, time_LDF);
          }
        }
      }
    }

    if (time_LDF == MAX_INT64 && time_GEB == MAX_INT64) {
      break;
    }

               
    if (time_LDF <= time_GEB) {
      unsigned int typeID = *((unsigned int*)(&dataLDF.type[0])); //crucify me for crimes of memory abuse
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
      //printf("writing LDF: %llu   %llu\n", time_LDF, time_GEB);
      continue;
    }
    else {
      WriteGEB(&outfile, &gebfile, &header);
      readGEB = 1;
      ctrGEB += 1;
      //printf("writing GEB: %llu   %llu\n", time_LDF, time_GEB);
      continue;
    }

    
  }

  gzclose(gebfile.gzf);
  gzclose(outfile.gzf);
  fclose(gebfile.ptr);
  fclose(ldffile);
  fclose(outfile.ptr);
  printf("\nLDF events: %llu read, %llu written\n", readCtrLDF, ctrLDF);
  printf("            %llu type-19, %llu type-22\n", type19, type22);  
  printf("  Bad MyRIAD: %llu (a), %llu (b)\n", type19_badmyriad_a, type19_badmyriad_b);  

  printf("GEB events: %llu read, %llu written, %llu type-19 discarded\n", readCtrGEB, ctrGEB, discards);
  printf("Done\n");
  
}
