#ifndef _RAIPHY_H_
#define _RAIPHY_H_

typedef struct {
  int kmerSize;
  long long ** count;
} kmer_freq_t;
  
kmer_freq_t* kiAllocKmerFreq(int kmerSize);
kmer_freq_t* kiFreeKmerFreq(kmer_freq_t* freq);
void kiClearKmerFreq(kmer_freq_t* freq);
void kiScanSeqsForFreq(alignment_t* seqs, kmer_freq_t* freq);
void kiScanFastaForFreq(char* fileName, kmer_freq_t* freq);
void kiFreqToRaiVectorOriginal(kmer_freq_t* freq, double* vector);

double* kiFreqToRai(kmer_freq_t* freq);


typedef struct {
  int nameLen;
  int kmerSize;
  int nClass;
  int nDim;
  int nSaved;
  char** names;
  double** vectors;
} rai_db_t;

rai_db_t* kiAllocRaiDb(int kmerSize);
rai_db_t* kiFreeRaiDb(rai_db_t* db);
void      kiRaiDbGrow(rai_db_t* db);
void      kiRaiDbAdd(rai_db_t* db, char* name, double* vector);
void      kiPrintRaiDb(FILE* fp, rai_db_t* db);
rai_db_t* loadDatabase(char* fileName);

int       classifySequenceOriginal(char* seq, rai_db_t* db, double* margin);
int       classifySequence(char* seq, rai_db_t* db, double* margin);

#endif /* _RAIPHY_H_ */


