/* hgFakeAgpForNcbi - Create fake AGP file by looking at N's using NCBI's specs*/
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "dnaseq.h"
#include "fa.h"
#include <stdbool.h>


int minContigGap = 25;
int minScaffoldGap = 50000;
int unknownSize = 100;
int minRangeToCallUnknown = 2000;
int maxRangeToCallUnknown = 50000;
bool wgs = false;

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "hgFakeAgpForNcbi - Create fake AGP (v1.1) file by looking at N's using NCBI's specs\n"
      "usage:\n"
      "   hgFakeAgpForNcbi input.fa output.agp\n"
      "options:\n"
      "   -minContigGap=N Minimum size for a gap between contigs.  Default %d\n"
      "   -WGS This is a wgs project, so mark contigs as 'W'. Default, marks as draft 'D'\n"
      "   -minScaffoldGap=N Min size for a gap between scaffolds. Default %d\n"
      "   -minRangeToCallUnknown=N Min size for a gap between scaffolds. Default %d\n"
      "   -maxRangeToCallUnknown=N Min size for a gap between scaffolds. Default %d\n"
      "   -unknownSize=N Min size for a gap between scaffolds. Default %d\n"
      , minContigGap, minScaffoldGap, minRangeToCallUnknown, maxRangeToCallUnknown, unknownSize
  );
}

static struct optionSpec options[] = {
    {"minContigGap", OPTION_INT},
    {"minScaffoldGap", OPTION_INT},
    {"minRangeToCallUnknown", OPTION_INT},
    {"maxRangeToCallUnknown", OPTION_INT},
    {"unknownSize", OPTION_INT},
    {"WGS", OPTION_BOOLEAN},
    {NULL, 0},
};

void agpContigLine(FILE *f, char *name, int seqStart, int seqEnd, int partIx, int contigIx)
/* Output AGP line for contig. */
{
  fprintf(f, "%s\t", name);
  fprintf(f, "%d\t", seqStart + 1);
  fprintf(f, "%d\t", seqEnd);
  fprintf(f, "%d\t", partIx);
  if(wgs)
    fprintf(f, "%c\t", 'W');
  else
    fprintf(f, "%c\t", 'D');
  fprintf(f, "%s_%d\t", name, contigIx);
  fprintf(f, "1\t");
  fprintf(f, "%d\t", seqEnd - seqStart);
  fprintf(f, "+\n");
}

char *strNotChar(char *s, char c)
/* Return pointer to first character in s that is *not* c, or NULL
 * if none such. */
{
  char x;
  for (;;)
  {
    x = *s;
    if (x == 0)
      return NULL;
    if (x != c)
      return s;
    ++s;
  }
}

int agpGapLine(FILE *f, char *name, int seqStart, int seqEnd, int gapSize, int lineIx)
/* Write out agp line for gap. */
{
  int offset = 0;
  bool unknown = false;
  if(gapSize >= minRangeToCallUnknown && gapSize < maxRangeToCallUnknown){
    offset = unknownSize - gapSize;
    gapSize = unknownSize;
    seqEnd = seqStart + unknownSize;
  }

  fprintf(f, "%s\t", name);
  fprintf(f, "%d\t", seqStart + 1);
  fprintf(f, "%d\t", seqEnd);
  fprintf(f, "%d\t", lineIx);
  if(unknown)
    fprintf(f, "%c\t", 'U');
  else
    fprintf(f, "%c\t", 'N');
  fprintf(f, "%d\t", gapSize);
  if (gapSize >= minScaffoldGap)
    fprintf(f, "contig\tno\n");
  else
    fprintf(f, "fragment\tyes\n");
  return(offset);
}

void fakeAgpForNcbiFromSeq(struct dnaSeq *seq, FILE *f)
/* Look through sequence and produce agp file. */
{
  //reset
  int partIx = 0;
  int contigIx = 0;
  char *gapStart, *gapEnd, *contigStart;
  char *dna, *seqStart, *seqEnd;
  int gapSize;
  int offset = 0;
  dna = seqStart = contigStart = seq->dna;
  seqEnd = seqStart + seq->size;
  for (;;)
  {
    gapStart = strchr(dna, 'n');
    if (gapStart == NULL)
    {
      agpContigLine(f, seq->name, contigStart - seqStart + offset, seqEnd - seqStart + offset,
          ++partIx, ++contigIx);
      break;
    }
    gapEnd = strNotChar(gapStart, 'n');
    if (gapEnd == NULL)
      gapEnd = seqEnd;
    gapSize = gapEnd - gapStart;
    if (gapSize >= minContigGap || gapEnd == seqEnd)
    {
      if (gapStart != contigStart)
        agpContigLine(f, seq->name, contigStart - seqStart + offset, gapStart - seqStart + offset,
            ++partIx, ++contigIx);
      offset += agpGapLine(f, seq->name, gapStart - seqStart + offset, gapEnd - seqStart + offset, gapSize, ++partIx);
      if (gapEnd == seqEnd)
        break;
      contigStart = gapEnd;
    }
    dna = gapEnd;
  }
}

void hgFakeAgpForNcbi(char *faIn, char *agpOut)
/* hgFakeAgpForNcbi - Create fake AGP file by looking at N's. */
{
  struct lineFile *lf = lineFileOpen(faIn, TRUE);
  FILE *f = mustOpen(agpOut, "w");
  struct dnaSeq seq;

  while (faSpeedReadNext(lf, &seq.dna, &seq.size, &seq.name))
  {
    fakeAgpForNcbiFromSeq(&seq, f);
  }

  carefulClose(&f);
  lineFileClose(&lf);
}

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 3)
    usage();
  minContigGap = optionInt("minContigGap", minContigGap);
  minScaffoldGap = optionInt("minScaffoldGap", minScaffoldGap);
  minRangeToCallUnknown = optionInt("minRangeToCallUnknown", minRangeToCallUnknown);
  maxRangeToCallUnknown = optionInt("maxRangeToCallUnknown", maxRangeToCallUnknown);
  unknownSize = optionInt("unknownSize", unknownSize);
  if(optionExists("WGS"))
    wgs = true;
  hgFakeAgpForNcbi(argv[1], argv[2]);
  return 0;
}
