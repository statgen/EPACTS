/**
 * 3. Output in VAT format
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "Argument.h"
#include "IO.h"
#include "StringUtil.h"
#include "LogFile.h"

#include "Type.h"
#include "Codon.h"
#include "GenomeSequence.h"
#include "Gene.h"
#include "GeneFormat.h"
#include "SequenceUtil.h"
#include "FreqTable.h"
#include "StringTemplate.h"
#include "GeneAnnotation.h"
#include "BedReader.h"
#include "GenomeScore.h"
#include "ModelParser.h"
#include "TabixReader.h"
#include "GitVersion.h"

void banner(FILE* fp) {
  const char* string =
      "..............................................         \n"
      " ...      Anno(tation)                       ...       \n"
      "  ...      Xiaowei Zhan, Goncalo Abecasis     ...      \n"
      "   ...      Speical Thanks:                    ...     \n"
      "    ...      Hyun Ming Kang, Yanming Li         ...    \n"
      "     ...      zhanxw@umich.edu                    ...  \n"
      "      ...      Sep 2011                            ... \n"
      "       ................................................\n"
      "                                                       \n"
      ;
  fputs(string, fp);
#ifdef DEBUG
  const char* debug =
      "-------------------------------------------------------\n"
      "|                                                     |\n"
      "|                DEBUG  MODE                          |\n"
      "|                (slow mode)                          |\n"
      "|                                                     |\n"
      "|   Try:                                              |\n"
      "|      make clean; make release                       |\n"
      "|   then run:                                         |\n"
      "|      ./executable/anno                              |\n"
      "|   to use release/fast version                       |\n"
      "|                                                     |\n"
      "-------------------------------------------------------\n"
      ;
  fputs(debug, fp);
#endif
};

extern const char* gitVersion;

OutputAnnotationString AnnotationString; // global variable


typedef enum {VCF = 0 , PLAIN, PLINK, EPACTS} InputFileFormat;

// hold input file
// and return (chrom, pos, ref, alt) iteratively
class AnnotationInputFile{
 public:
  AnnotationInputFile(const char* inputFileName, const char* inputFormatStr) {
    // check inputFormat
    std::string inputFormat = toLower(inputFormatStr);
    if (!inputFormat.empty() &&
        inputFormat != "vcf" &&
        inputFormat != "plain" &&
        inputFormat != "plink" &&
        inputFormat != "epacts") {
      fprintf(stderr, "Unsupported input format [ %s ], we support VCF, plain, plink and EPACTS formats.\n", inputFormatStr);
      LOG << "Unsupported input format [ " << inputFormatStr << " ], we support VCF, plain, plink and EPACTS formats.\n";
      abort();
    };

    if (inputFormat == "vcf" || inputFormat.empty()) {
      // ga.annotateVCF(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = VCF;
    } else if (inputFormat == "plain") {
      // ga.annotatePlain(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = PLAIN;
    } else if (inputFormat == "plink") {
      // ga.annotatePlink(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = PLINK;
    } else if (inputFormat == "epacts") {
      // ga.annotateEpacts(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = EPACTS;
    } else{
      fprintf(stderr, "Cannot recognize input file format: %s \n", inputFileName);
      abort();
    };

    // open input files
    this->lr = new LineReader(inputFileName);
  };
  ~AnnotationInputFile() {
    this->close();
  }
  void close() {
    if (this->lr) {
      delete lr;
      this->lr = NULL;
    }
  };


  int openReferenceGenome(const char* referenceGenomeFileName) {
    return this->gs.open(referenceGenomeFileName);
  }
  // // extract one header line
  // bool extractHeader(std::string* l) {
  //   this->line.clear();
  //   bool ret = this->lr->readByLine(this->line);
  //   *l = this->line;
  //   return ret;
  // };
  // // extract one header line
  // bool extractHeader(std::string* l) {
  //   this->line.clear();
  //   bool ret = this->lr->readByLine(this->line);
  //   *l = this->line;
  //   return ret;
  // };

  // check ref and alt alleles (may switch ref and alt)
  //@return true if (1) ref match reference; (2) after switch ref and alt, match reference genome.
  bool forceReferenceStrand(const std::string& chrom,
                            const int& pos,
                            std::string* ref,
                            std::string* alt) {
    // determine ref base from reference
    bool refMatchRef = true;
    for (size_t i = 0; i < ref->size(); i++ ) {
      if ((*ref)[i] != gs[chrom][pos - 1 + i]) {
        refMatchRef = false;
        break;
      }
    }
    if (!refMatchRef) {
      bool altMatchRef = true;
      for (size_t i = 0; i < alt->size(); i++ ) {
        if ( (*alt)[i] != gs[chrom][pos - 1 + i]) {
          altMatchRef = false;
          break;
        }
      }
      if (!altMatchRef) {
        fprintf(stderr, "Ref [%s] and alt [%s] does not match reference: %s:%d\n", ref->c_str(), alt->c_str(),
                chrom.c_str(), pos);
        return false;
      } else {
        std::swap(*ref, *alt);
      }
    }
    return true;
  }
  void setCheckReference(bool b) {
    this->checkReference = b;
  };
  // if reach end or experience something wrong, will @return false
  bool extract(std::string* chrom,
               int* pos,
               std::string* ref,
               std::string* alt,
	       std::string* markerId
	       ) {
    bool ret;
    *markerId = ".";
    do {
      ret = this->lr->readLine(&this->line);
      if (ret == false)
        return ret;
    } while (this->line.empty());

    // for any line beginning with '#', store headers
    while (line[0] == '#') {
      this->header.push_back(line);
      do {
        ret = this->lr->readLine(&this->line);
        if (ret == false)
          return ret;
      } while (this->line.empty());
    }


    stringTokenize(line, "\t ", &fd);

    switch (this->format){
      case VCF:
        if (fd.size() < 5) return false;
        *chrom = fd[0];
        *pos = toInt(fd[1]);
	*markerId = fd[2];
        *ref = fd[3];
        *alt = fd[4];
        break;
      case PLAIN:
        if (fd.size() < 4) return false;
        *chrom = fd[0];
        *pos = toInt(fd[1]);
        *ref = fd[2];
        *alt = fd[3];
        break;
      case PLINK:
        if (fd.size() < 10) return false;
        *chrom = fd[0];
        *pos = toInt(fd[2]);
        *ref = fd[3];
        *alt = fd[6];

        if (!forceReferenceStrand(*chrom, *pos, ref, alt))
          return false;
        break;
      case EPACTS:
        {
          // e.g.
          // 20 139681 139681 20:139681_G/A   266     1       1       0.0018797       NA      NA
          int beg = 0;
          int sep = fd[3].find(':', beg);
          *chrom = fd[3].substr(beg, sep - beg);

          beg = sep + 1;
          sep = fd[3].find('_', beg);
          *pos = toInt(fd[3].substr(beg, sep - beg));

          beg = sep + 1;
          sep = fd[3].find('/', beg);
          *ref = toUpper(fd[3].substr(beg, sep - beg));

          beg = sep + 1;
          sep = fd[3].find_first_of(" _\t", beg);
          *alt = toUpper(fd[3].substr(beg, sep - beg));

          epactsPrefixLength = sep;
          if ( chrom->empty() || *pos <= 0 || ref->empty() || alt->empty()) {
            fprintf(stderr, "Skip line: %s ...." , fd[3].c_str());
            LOG << "Skip: " << fd[3];
            return false;
          }
        }
        break;
      default:
        fprintf(stderr, "Unknown format, quitting!\n");
        abort();
        break;
    }// end switch

    // verify reference
    if (this->checkReference) {
      std::string refFromGenome = this->gs.getBase(*chrom, *pos, *pos + ref->size());
      if ((*ref) != refFromGenome) {
        fprintf(stdout, "ERROR: Reference allele does not match genome reference [ %s:%d %s]\n", chrom->c_str(), *pos, ref->c_str());
        LOG << "ERRROR: Reference allele [" << ref <<   "]  does not match reference genome [" << refFromGenome << "] at " << *chrom << ":" << *pos << "\n";
      };
    }
    return true;
  };
  InputFileFormat getFormat() const {
    return this->format;
  };
  std::string getEpactsPrefix(const std::string& s) const{
    return s.substr(0, this->epactsPrefixLength);
  };
  const std::vector< std::string>& getHeader() const {
    return this->header;
  };
  const std::vector< std::string>& getFields() const {
    return this->fd;
  };
 private:
  bool checkReference;
  InputFileFormat format;
  LineReader* lr;
  std::vector< std::string> fd;
  std::string line;
  std::vector< std::string> header;
  GenomeSequence gs; // check if ref alleles matched reference genome
  size_t epactsPrefixLength;
}; // end class AnnotationInputFile

class AnnotationOutputFile{
 public:
  AnnotationOutputFile(const char* out):headerOutputted(false), totalVariants(0),outputFileName(out) {
    this->fout = fopen(out, "wt");
    if (!this->fout) {
      fprintf(stderr, "Cannot open otuput file %s for write.\n", out);
    };
  };
  ~AnnotationOutputFile() {
    this->close();
  }
  void close() {
    if (this->fout) {
      fprintf(stdout, "DONE: %d varaints are annotated.\n", totalVariants);
      fprintf(stdout, "DONE: Generated annotation output in [ %s ].\n", outputFileName.c_str());
      fclose(this->fout);
      this->fout = NULL;
    }
  };
  void linkToInput(const AnnotationInputFile& in) {
    this->aif =  &in;
  };
  void writeHeader() {
    writeHeader(this->aif->getHeader());
  };
  void writeHeader(const std::vector< std::string> & h) {
    for (size_t i = 0; i < h.size(); ++i) {
      fputs(h[i].c_str(), this->fout);
      fputc('\n', this->fout);
    }
  }
  void writeResult(const OrderedMap<std::string, std::string>& res) {
    //fprintf(stderr, "Need to link output to input!\n");
    assert(aif);

    if (!this->headerOutputted) {
      this->writeHeader();
      this->headerOutputted = true;
    }

    const std::vector< std::string> & field = aif->getFields();
    std::string key;
    std::string val;
    switch (aif->getFormat()) {
      case VCF:
        {
          for (size_t i = 0; i < field.size(); ++i) {
            if (i) fputc('\t', fout);
            fputs(field[i].c_str(), fout) ;
            if (i == 7) { // 7: the 8th column in 0-based index
              if (!field[i].empty())
                fputc(VCF_INFO_SEPARATOR[0], fout);
              for (size_t j = 0; j < res.size(); ++j) {
                if (j)
                  fputc(VCF_INFO_SEPARATOR[0], fout);
                res.at(j, &key, &val);
                fputs(key.c_str(), this->fout);
                fputc('=', this->fout);
                fputs(val.c_str(), this->fout);
              }
            }
          };
        }
        break;
      case PLAIN:
      case PLINK:
        {
          for (size_t i = 0; i < field.size(); ++i) {
            if (i) fputc('\t', fout);
            fputs(field[i].c_str(), fout);
          };
          for (size_t j = 0; j < res.size(); ++j) {
            fputc('\t', fout);
            res.at(j, &key, &val);
            // fputs(key.c_str(), this->fout);
            // fputc('=', this->fout);
            fputs(val.c_str(), this->fout);
          }
        }
        break;
      case EPACTS:
        {
          //fputs(field[0].substr(0, sep).c_str(), fout);
          //fputs(aif->getEpactsPrefix(field[3]).c_str(), fout);
          //fputc('_', fout);
          //std::string s;
          //res.value("ANNO", &s);
          //fputs(s.c_str(), fout);

          for (unsigned int i = 0; i < field.size(); i++ ) {
            if ( i > 0 ) fputs("\t", fout);
	    if ( i == 3 ) {
	      fputs(aif->getEpactsPrefix(field[i]).c_str(), fout);
	      fputc('_', fout);
	      std::string s;
	      res.value("ANNO", &s);
	      fputs(s.c_str(), fout);
	    }
	    else 
	      fputs(field[i].c_str(), fout) ;
          }
        }
        break;
      default:
        break;
    };
    fputc('\n', this->fout);
    ++this->totalVariants;
  };
 private:
  bool headerOutputted;
  const AnnotationInputFile* aif;
  FILE* fout;
  int totalVariants;
  std::string outputFileName;

}; // class AnnotationOutputFile

// run annotations (gene based, bed file, genomeScore, tabix database)
// and store results
class AnnotationController{
 public:
  AnnotationController(AnnotationInputFile& in): aif(in){
  };
  virtual ~AnnotationController() {
    for (size_t i = 0; i < bedReader.size() ; ++i) {
      delete bedReader[i];
    }
    for (size_t i = 0; i < genomeScore.size(); ++i ) {
      delete genomeScore[i];
    };
    for (size_t i = 0; i < tabixReader.size(); ++i ) {
      delete tabixReader[i];
    };

    
  };
  void openBedFile(const char* tag, const char* fn) {
    // check duplication
    for (size_t i = 0; i < this->bedTag.size(); ++i ) {
      if (this->bedTag[i] == tag) {
        fprintf(stderr, "ERROR: Duplicated tag [ %s ] \n", tag);
        return;
      }
    }

    // add bedFile
    BedReader* p = new BedReader;
    int ret = p->open(fn);
    if (ret < 0) {
      fprintf(stderr, "Cannot open BED file: [ %s ]\n", fn);
      delete p;
      return ;
    } else {
      fprintf(stderr, "DONE: Load %d regions from BED file\n", ret);
    };

    this->bedTag.push_back(tag);
    this->bedReader.push_back(p);
  };
  void openGenomeScoreFile(const char* tag, const char* fn){
    // check duplication
    for (size_t i = 0; i < this->genomeScoreTag.size(); ++i) {
      if (this->genomeScoreTag[i] == tag ) {
        fprintf(stderr, "ERROR: Duplicated tag [ %s ] \n", tag);
        return;
      }
    }

    // add genome score file
    GenomeScore* p = new GenomeScore(fn);
    this->genomeScoreTag.push_back(tag);
    this->genomeScore.push_back(p);
  };

  void addTabixReader(TabixReader* t) {
    this->tabixReader.push_back(t);
  };

  void annotate(std::string& chrom,
                int& pos,
                std::string& ref,
                std::string& alt,
		std::string& markerId
		) {
    this->geneAnnotation.annotate(chrom, pos, ref, alt, markerId);

    this->result.clear();
    this->result["ANNO"] = this->geneAnnotation.getTopPriorityAnnotation();
    this->result["ANNOFULL"] = this->geneAnnotation.getFullAnnotation();

    std::vector<std::string> bedString;
    for (size_t br = 0; br < this->bedReader.size(); ++ br) {
      if (this->bedReader[br]->find(chrom.c_str(), pos, &bedString)){
        if (!bedString.empty())  {
          this->result[this->bedTag[br]] = stringJoin(bedString, ',');
        } else {
          this->result[this->bedTag[br]] = "";
        }
      }
    }
    for (size_t gs = 0; gs < this->genomeScore.size(); ++ gs) {
      this->result[this->genomeScoreTag[gs]] = toStr(this->genomeScore[gs]->baseScore(chrom.c_str(), pos));
    }
    for (size_t tb = 0; tb < this->tabixReader.size(); ++ tb) {
      TabixReader& tabix = *this->tabixReader[tb];
      tabix.addAnnotation(chrom, pos, ref, alt);
      size_t s = tabix.getTag().size();
      for (size_t i = 0; i < s; ++i) {
        this->result[tabix.getTag()[i]] = tabix.getAnnotation()[i];
      }
    }
  };

  const OrderedMap<std::string, std::string>& getResult() const {
    return this->result;
  }
 public:
  AnnotationInputFile& aif;
  // various annotation types
  GeneAnnotation geneAnnotation;
 private:
  // various annotation types
  std::vector<BedReader*> bedReader;
  std::vector<std::string> bedTag;
  std::vector<GenomeScore*> genomeScore;
  std::vector<std::string> genomeScoreTag;
  std::vector<TabixReader*> tabixReader;
  OrderedMap<std::string, std::string> result; // store all types of annotation results
};


int main(int argc, char *argv[])
{
  banner(stdout);

  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Required Parameters")
      ADD_STRING_PARAMETER(pl, inputFile, "-i", "Specify input VCF file")
      ADD_STRING_PARAMETER(pl, outputFile, "-o", "Specify output VCF file")
      ADD_PARAMETER_GROUP(pl, "Gene Annotation Parameters")
      ADD_STRING_PARAMETER(pl, geneFile, "-g", "Specify gene file")
      ADD_STRING_PARAMETER(pl, referenceFile, "-r", "Specify reference genome position")      
      ADD_STRING_PARAMETER(pl, inputFormat, "--inputFormat", "Specify format (default: vcf). \"-f plain \" will use first 4 columns as chrom, pos, ref, alt")
      ADD_BOOL_PARAMETER(pl, checkReference, "--checkReference", "Check whether reference alleles matches genome reference")
      ADD_STRING_PARAMETER(pl, geneFileFormat, "-f", "Specify gene file format (default: refFlat, other options: knownGene, refGene)")
      ADD_STRING_PARAMETER(pl, priorityFile, "-p", "Specify priority of annotations")
      ADD_STRING_PARAMETER(pl, codonFile, "-c", "Specify codon file (default: codon.txt)")
      ADD_INT_PARAMETER(pl, upstreamRange, "-u", "Specify upstream range (default: 50)")
      ADD_INT_PARAMETER(pl, downstreamRange, "-d", "Specify downstream range (default: 50)")
      ADD_INT_PARAMETER(pl, spliceIntoExon, "--se", "Specify splice into extron range (default: 3)")
      ADD_INT_PARAMETER(pl, spliceIntoIntron, "--si", "Specify splice into intron range (default: 8)")
      ADD_STRING_PARAMETER(pl, outputFormat, "--outputFormat", "Specify predefined annotation words (default or epact)")
      ADD_PARAMETER_GROUP(pl, "Other Annotation Tools")
      ADD_STRING_PARAMETER(pl, genomeScore, "--genomeScore", "Specify the folder of genome score (e.g. GERP=dirGerp/,SIFT=dirSift/)")
      ADD_STRING_PARAMETER(pl, bedFile, "--bed", "Specify the bed file and tag (e.g. ONTARGET1=a1.bed,ONTARGET2=a2.bed)")
      ADD_STRING_PARAMETER(pl, tabixFile, "--tabix", "Specify the tabix file and tag (e.g. abc.txt.gz(chrom=1,pos=7,ref=3,alt=4,SIFT=7,PolyPhen=10)")      
      END_PARAMETER_LIST(pl)
      ;

  pl.Read(argc, argv);
  if (FLAG_geneFileFormat.empty()) {
    FLAG_geneFileFormat = "refFlat";
  }

  if (FLAG_inputFile.empty()) {
    pl.Help();
    fprintf(stderr, "Please specify input file\n");
    exit(1);
  }
  if (FLAG_outputFile.empty()) {
    pl.Help();
    fprintf(stderr, "Please specify output file\n");
    exit(1);
  }

  GeneAnnotationParam param;
  param.upstreamRange = FLAG_upstreamRange ? FLAG_upstreamRange : 50;
  param.downstreamRange = FLAG_downstreamRange ? FLAG_downstreamRange : 50;
  param.spliceIntoExon = FLAG_spliceIntoExon ? FLAG_spliceIntoExon : 3;
  param.spliceIntoIntron = FLAG_spliceIntoIntron ? FLAG_spliceIntoIntron : 8;

  std::string logFileName = FLAG_outputFile + ".log";
  LOG_START(logFileName.c_str());
  LOG_START_TIME;
  LOG_PARAMETER(pl);
  LOG << "Version: " << gitVersion << "\n";
  
  pl.Status();
  if (FLAG_REMAIN_ARG.size() > 0){
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    abort();
  }
  
  if (FLAG_geneFile.empty()) {
    pl.Help();
    fprintf(stderr, "Please specify gene file\n");
    exit(1);
  }
  
  if (FLAG_priorityFile.empty()) {
    fprintf(stderr, "Use default priority file: /net/fantasia/home/zhanxw/anno/priority.txt\n");
    FLAG_priorityFile = "/net/fantasia/home/zhanxw/anno/priority.txt";
  };
  if (FLAG_codonFile.empty()) {
    fprintf(stderr, "Use default codon file: /net/fantasia/home/zhanxw/anno/codon.txt\n");
    FLAG_codonFile = "/net/fantasia/home/zhanxw/anno/codon.txt";
  }
  if (FLAG_referenceFile.empty()) {
    fprintf(stderr, "Use default codon file: /net/fantasia/home/zhanxw/anno/codon.txt\n");
    FLAG_referenceFile = "/data/local/ref/karma.ref/human.g1k.v37.fa";
  }

  if (!FLAG_outputFormat.empty()) {
    AnnotationString.setFormat(FLAG_outputFormat.c_str());
  }
  else {
    AnnotationString.setFormat("default");
  };


  AnnotationInputFile aif(FLAG_inputFile.c_str(), FLAG_inputFormat.c_str());
  aif.openReferenceGenome(FLAG_referenceFile.c_str());
  aif.setCheckReference(FLAG_checkReference);

  AnnotationController controller(aif);


  controller.geneAnnotation.setAnnotationParameter(param);
  controller.geneAnnotation.openReferenceGenome(FLAG_referenceFile.c_str());
  controller.geneAnnotation.openCodonFile(FLAG_codonFile.c_str());
  controller.geneAnnotation.openPriorityFile(FLAG_priorityFile.c_str());
  // controller.geneAnnotation.setFormat(FLAG_geneFileFormat);
  controller.geneAnnotation.openGeneFile(FLAG_geneFile.c_str(), FLAG_geneFileFormat.c_str());

  if (!FLAG_bedFile.empty()) {
    fprintf(stderr, "Use bed file: %s\n", FLAG_bedFile.c_str() );
    std::vector< std::string> fd;
    std::vector< std::string> bed;
    stringTokenize(FLAG_bedFile, ",", &fd);
    for (size_t i = 0; i < fd.size(); ++i ){
      stringTokenize(fd[i], "=", &bed);
      if (bed.size() == 2) {
        controller.openBedFile(bed[0].c_str(), bed[1].c_str());
      } else {
        fprintf(stderr, "ERROR: Cannot recognized format [ %s ].\n", fd[i].c_str());
        exit(1);
      };
    };
  };
  if (!FLAG_genomeScore.empty()){
    // fprintf(stderr, "Use binary GERP score: %s\n", FLAG_genomeScore.c_str());
    // ga.addGenomeScore("GERP", FLAG_genomeScore.c_str());
    fprintf(stderr, "Use binary score file: %s\n", FLAG_genomeScore.c_str() );
    std::vector< std::string> fd;
    std::vector< std::string> gs;
    stringTokenize(FLAG_genomeScore, ",", &fd);
    for (size_t i = 0; i < fd.size(); ++i ){
      stringTokenize(fd[i], "=", &gs);
      if (gs.size() == 2) {
        controller.openGenomeScoreFile(gs[0].c_str(), gs[1].c_str());
      } else {
        fprintf(stderr, "ERROR: Cannot recognized format [ %s ].\n", fd[i].c_str());
        exit(1);
      };
    };
  }
  // parse something like:
  // abc.txt.gz(chrom=1,pos=7,ref=3,alt=4,SIFT=7,PolyPhen=10)
  if(!FLAG_tabixFile.empty()){
    fprintf(stderr, "Use tabix file: %s\n", FLAG_tabixFile.c_str() );    
    ModelParser mp;
    mp.parse(FLAG_tabixFile);
    std::string fn = mp.getName();
    int chrom, pos, ref, alt;
    mp.assign("chrom", &chrom).assign("pos", &pos).assign("ref", &ref).assign("alt", &alt);
    fprintf(stderr, "Column %d, %d, %d and %d in tabix file will be matched to chromosome, position, reference allele, alternative allele respectively.\n", chrom, pos, ref, alt);
    TabixReader* tabix = new TabixReader(fn.c_str(), chrom, pos, ref, alt);
    
    for (size_t i = 0; i < mp.getParam().size(); ++i) {
      if ( toLower(mp.getParam()[i]) == "chrom" ||
           toLower(mp.getParam()[i]) == "pos" ||
           toLower(mp.getParam()[i]) == "ref" ||
           toLower(mp.getParam()[i]) == "alt") {
        continue;
      }
      int intValue;
      if (str2int(mp.getValue(i), &intValue)) {
        tabix->addTag(mp.getParam()[i], intValue);
        fprintf(stderr, "Tag %s will be from column %d in tabix file\n", mp.getParam()[i].c_str(), intValue);
      } else {
        tabix->addTag(mp.getParam()[i], mp.getValue(i));
        fprintf(stderr, "Tag %s will be from column %s (from header) in tabix file\n", mp.getParam()[i].c_str(), mp.getValue(i).c_str());        
      }
    }
    controller.addTabixReader(tabix);
  }
  // if (inputFormat == "vcf" || FLAG_inputFormat.size() == 0) {
  //   ga.annotateVCF(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
  // } else if (inputFormat == "plain") {
  //   ga.annotatePlain(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
  // } else if (inputFormat == "plink") {
  //   ga.annotatePlink(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
  // } else if (inputFormat == "epacts") {
  //   ga.annotateEpacts(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
  // } else{
  //   fprintf(stderr, "Cannot recognize input file format: %s \n", FLAG_inputFile.c_str());
  //   abort();
  // };

  std::string chrom;
  int pos;
  std::string ref;
  std::string alt;
  std::string markerId;
  AnnotationOutputFile aof(FLAG_outputFile.c_str());
  aof.linkToInput(aif);
  while (aif.extract(&chrom, &pos, &ref, &alt, &markerId)) {
    controller.annotate(chrom, pos, ref, alt, markerId);
    aof.writeResult(controller.getResult());
  };
  // aof.writeResult(controller.getResult()); // TODO: will add this to handle when input only have comment lines
  
  // output stats
  controller.geneAnnotation.outputAnnotationStats(FLAG_outputFile.c_str());

  aof.close();
  aif.close();
  LOG << "Annotate " << FLAG_inputFile << " to " << FLAG_outputFile << " succeed!\n";

  LOG_END_TIME;
  LOG_END ;
  printf("Annotation succeed!\n");
  return 0;
}
