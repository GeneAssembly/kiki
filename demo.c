
#include <locale.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ki.h"
#include "debug.h"

void testArgs() {
  int new_kmer_id;
  char new_kmer[100];
  int array[100], new_array[100];
  int i, n;
  for (i = 0; i < 10; i++) array[i] = i;
    
  kiUserTestArgs(10, array, "ABCDEFGHIJKLMN", 77, &n, new_array, new_kmer, &new_kmer_id);

  printf("new_array = ");
  for (i = 0; i < n; ++i) {
    printf("%d ", new_array[i]);
  }
  printf("\n");
  printf("new_kmer = '%s', new_kmer_id = %d\n", new_kmer, new_kmer_id);
}

void testArgs2() {
  int cmd, i;
  for (i = 100; i < 103; ++i) {
    cmd = i;
    kiUserCall(cmd);
  }
  testArgs();
}

void testLoadSeq() {
  seq_id_t id;

  kiUserLoadSeq("ATTGCAACCGACCCGCTCG", &id);
  printf("seq_id 1 = (%d, %d)\n", id.id, id.cpu);

  kiUserLoadSeq("TGCAACCGACCCGCTCGCA", &id);
  printf("seq_id 2 = (%d, %d)\n", id.id, id.cpu);
}

void testLoadFasta() {
  long long nSeqs = 0;
  kiUserSetHashMode(KI_HASH_ENDS);
  kiUserLoadFasta("short.fa", &nSeqs);
  kiUserLoadFasta("medium.fa", &nSeqs);
  /* kiUserLoadFasta("SMG-5.fa", &nSeqs); */
  printf("Read %lld sequences from fasta.\n", nSeqs);
}

void testMemInfo() {
  long used_sum;
  kiUserMemInfo(&used_sum);

  setlocale(LC_NUMERIC, "en_US.utf-8");
  printf("Total memory used = %'ld KB\n", used_sum >> 10);
}

void testSetK() {
  int k = 31;
  kiUserSetK(k);
}

void testGetKmer() {
  int nSeqs;
  char kmer1[] = "TGGTTTTCTTTCGGG";
  char kmer2[] = "TGCAACCGACCCGCT";
  char kmer3[] = "ATTGCAACCGACCCG";

  kiUserGetKmer(kmer1, &nSeqs);
  printf("%d sequences match '%s'.\n", nSeqs, kmer1);

  kiUserGetKmer(kmer2, &nSeqs);
  printf("%d sequences match '%s'.\n", nSeqs, kmer2);

  kiUserGetKmer(kmer3, &nSeqs);
  printf("%d sequences match '%s'.\n", nSeqs, kmer3);
}

void userTests() {
  testLoadSeq();
  testLoadFasta();
  testGetKmer();
  /* testSetK(); */
  /* testArgs(); */
  /* testArgs2(); */
  /* testMemInfo(); */
}

void dummyAssemble() {
  long long nSeqs;
  /* char* fileName = "short.fa"; */
  /* char* fileName = "short11.fa"; */
  /* char* fileName = "medium.fa"; */
  char* fileName = "medium11.fa";
  /* char* fileName = "SMG-5.fa"; */
  kiUserLoadFasta(fileName, &nSeqs);
  printf("Loaded %lld sequences from '%s'.\n", nSeqs, fileName);

  kiUserDummyAssemble(fileName);
}

void exactAssemble(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "No input fasta files specified.\n");
    exit(1);
  }
  long long nSeqs = 0, nTotalSeqs = 0;
  char* fileName = NULL;
  int i;
  for (i = 1; i < argc; ++i) {
    fileName = argv[i];
    kiUserLoadFasta(fileName, &nSeqs);
    printf("Loaded %lld sequences from '%s'.\n", nSeqs, fileName);
    nTotalSeqs += nSeqs;
  }
  if (i > 2) printf("Total sequences loaded = %lld\n", nTotalSeqs);

  char* dot = strrchr(fileName, '.');
  if (dot != NULL) *dot = '\0';
  
  kiUserDummyAssemble(fileName);
}


void hybridAssemble(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "No input fasta files specified.\n");
    exit(1);
  }
  long long nSeqs = 0, nTotalSeqs = 0;
  char* fileName = NULL;
  int i;
  for (i = 1; i < argc; ++i) {
    fileName = argv[i];
    kiUserLoadFasta(fileName, &nSeqs);
    printf("Loaded %lld sequences from '%s'.\n", nSeqs, fileName);
    nTotalSeqs += nSeqs;
  }
  if (i > 2) printf("Total sequences loaded = %lld\n", nTotalSeqs);

  char* dot = strrchr(fileName, '.');
  if (dot != NULL) *dot = '\0';
  
  /* kiUserHybridAssemble(); */
}

void userMain(int argc, char *argv[]) {
  /* userTests(); */
  /* dummyAssemble(); */
  exactAssemble(argc, argv);
}

int main(int argc, char *argv[])
{
  kiInit(&argc, &argv);
  kiStart();

  if (kiIsUserDomain())
    userMain(argc, argv);
  
  kiStop();
  kiFinalize();
  return 0;
}
