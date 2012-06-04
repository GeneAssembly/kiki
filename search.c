
#include "extern.h"

double global_time = 0.;

void checkTime(char* info) {
  if (global_time < 0.1) {
    global_time = MPI_Wtime();
    return;
  }
  double time = MPI_Wtime();
  double elapsed = time - global_time;
  global_time = time;
  kipm("Time: %s (%.2f s elapsed)\n", info, elapsed);
}


void testHash(int argc, char *argv[]) {
  char kmer1[] = "TGGTTTTCTTTCGGG";
  char kmer2[] = "TGCAACCGACCCGCT";
  char kmer3[] = "ABCAAAAAAAAAAAA";
  /* char kmer3[] = "ATBAAAAAAAAAAAA"; */
  /* char kmer3[] = "ATTGCAACCGACCCG"; */
  char* strings[3];
  strings[0] = kmer1;
  strings[1] = kmer2;
  strings[2] = kmer3;

  hashtable_t* hash = kiMakeHashtable(strings, 3, 16, false, true);
  /* int i; */
  /* for (i = 0; i < 3; ++i) { */
  /*   kiHashtableAddCopy(hash, strings[i], strlen(strings[i]), i+3, 0); */
  /* } */
  hash = kiFreeHashtable(hash);
  /* pm("nChains = %d\n", hash->nChains); */
}

void testPairwiseAlignment(int argc, char *argv[]) {
  /* char s1[] = "MSVMYKKILYPTDFSETAEIA"; */
  /* char s2[] = "SVMYSETAEIA"; */

  char s1[] = "MSVMYKKILYPTDFSETAEIALKHVKAFKTLKAEEVILLHVIDEREIKKRDIFSLLLGVAGLNKSVEEFENELKNKLTEEAKNKMENIKKELEDVGFKVKDIIVVGIPHEEIVKIAEDEGVDIIIMGSHGKTNLKEILLGSVTENVIKKSNKPVLVVKRKNS";
  /* char s2[] = "SVMYKKILYPTDFSETAEIALKHVKAFKTLKAEEVILLHVIDEREIKSVEEFENELKNKLTEEAKNKMENIKKELEDVGFKVKDIIVVGIPHEEIVKIAEDEGVDIIIMGSHGKTNLKEILLGSVTENVIKKSNKPVLVVKRKNS"; */

  /* char s1[] = "MSVMYKKILYPTDFSETAEIALKHVKAFKTLKAEEVILLHVIDEREIKKRDIFSLLLGVAGLNKSVEEFENELKNKLTEEAKNKMENIKKELEDVGFKVKDIIVVGIPHEEIVKIAEDEGVDIIIMGSHGKTNLKEILLGSVTENVIKKSNKPVLVVK"; */
  /* char s2[] = "MIFMFRKVLFPTDFSEGAYRAVEVFEKRNKMEVGEVILLHVIDEGTLEELMDGYSFFYDNAEIELKDIKEKLKEEASRKLQEKAEEVKRAFRAKNVRTIIRFGIPWDEIVKVAEEENVSLIILPSRGKLSLSHEFLGSTVMRVLRKTKKPVLIIK"; */
  
  /* char s1[] = "MSVMYKKILYPTDFSETAEIALKHVKAFKTLKAEEVILLHVIDEREIKKRDIFSLLLGVAGLNKSVEEFENELKNKLTEEAKNKMENIKKELEDVGFKVKDIIVVGIPHEEIVKIAEDEGVDIIIMGSHGKTNLKEILLGSVTENVIKKSNKPVLVVK"; */
  /* char s2[] = "MYSKILLPTDGSKQANKAAEHAIWIARESGAEIIALTVMETSSLVGLPADDLIIRLREMLEEEASRSLEAVKKLVEESGADIKLTVRTDEGSPAEAILRTVEKEGVDLVVMGTSGKHGLDRFLLGSVAEKVVRSAGCPVLVV"; */

  char s2[] = "MYRKILLATDGSECSMQAAGYAIETAAQNRAELLALTVTETYPLDNLPVEELTRKVTELFRKESEEALQKVEDLAVSLDTPVKVRKMMVDGSPAETILKVADEENVDLIVVGASGKHALERFLLGSVSEKIVRHARVPVLVVHSK"; /* 197 */

  /* char s2[] = "MFEKIMVPTDGSEYAARAEDMAIELAGRLGSVVIAVHVIDEKLIYPFDVLEDEGKEILASVQRKGREAGVQVDEVLVFGSPAHDMKKITEKTGADLVVIASHGRSGLEKLLMGSVAETTLKTVDVPVLLVK"; /\* 195 *\/ */
  

  /* char s1[] = "PHEEIVKIAEDEGVDIIIMGSHGKTNLKEILLGSVTENVIKKSNKPVLVVK"; */
  /* char s2[] = "PPDDIVDFADEVDAIAIVIGIRKRSPTGKLIFGSVARDVILKANKPVICIK"; */

  /* char s1[] = "FRKVLFPTDFSEGAY"; */
  /* char s2[] = "FRKILIPPDFSQALH"; */

  float score1 = kiPairwiseAlignmentLocal(s1, s2);
  float score2 = kiPairwiseAlignmentGlobal(s1, s2);
  float score3 = kiPairwiseAlignmentSimpleGap(s1, s2);
  /* float score = kiPairwiseAlignmentSimpleGap(s1, s2); */
  kipm("Alignment scores: local = %.2f,  global = %.2f,  simple gapped = %.2f\n", score1, score2, score3);
}


void proteinSearch(int argc, char *argv[]) {
  char usage[] = "search -p profile -i database -s cutoff";

  /* Read sequences from multiple files */
  int nSeqs = 0, nTotalSeqs = 0;
  char* fileName = NULL;
  char* profileName = NULL;

  int k = 8;
  float cutoff = 50.;
  
  int i;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-k") == 0) {
      if (i < argc) k = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-s") == 0) { 
      if (i < argc) cutoff = (float)atof(argv[++i]);
    } else if (strcmp(argv[i], "-i") == 0) { 
      while (i+1 < argc && argv[i+1][0] != '-') {
        fileName = argv[++i];
        if (kiIsDomainRoot()) {
          kiUserReadFasta(fileName, &nSeqs);
          /* kiUserLoadFasta(fileName, &nSeqs); */
          /* printf("Read %d sequences from '%s'.\n", nSeqs, fileName); */
          /* nTotalSeqs += nSeqs; */
        }
      }
    } else if (strcmp(argv[i], "-I") == 0) { 
      fileName = argv[++i];
      char faFileName[KI_NAME_BUF_SIZE];
      if (kiIsDomainRoot()) {
        FILE* fp = fopen(fileName, "r");
        if (fp != NULL) {
          while (fscanf(fp, "%s", faFileName) != EOF) {
            kipmsg(3, "faFileName = %s\n", faFileName);
            kiUserReadFasta(faFileName, &nSeqs);
          }
        }
        fclose(fp);
      }
    } else if (strcmp(argv[i], "-p") == 0 && i+1 < argc) {
      profileName = argv[++i];
    } 
  }
  if (nTotalSeqs > nSeqs && kiIsDomainRoot())
    printf("Total sequences read = %d\n", nTotalSeqs);

  if (!fileName || !profileName) {
    fprintf(stderr, "Usage: %s\n", usage);
    pmSetLevel(0);
    return;
  }

  kiUserSetK(k);
  
  alignment_t* hits = kiAllocAlignment();
  kiUserSearchProfile(profileName, cutoff, hits);
  
  kiFreeAlignment(hits);
}

void userMain(int argc, char *argv[]) {
  /* testHash(argc, argv); */
  /* testPairwiseAlignment(argc, argv); */
  proteinSearch(argc, argv);
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
