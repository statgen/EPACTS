#ifndef _TYPE_H_
#define _TYPE_H_


typedef enum {
  SNP = 0,
  INS,
  DEL,
  MIXED,
  SV,
  NO_VARIATION, // monomorphic site
  UNKNOWN = 99
} VARIATION_TYPE;

typedef enum {
  STRUCTURE_VARIATION = 0,
  STOP_GAIN,
  STOP_LOSS,
  START_GAIN,
  START_LOSS,
  FRAME_SHIFT,            /* Indel length is not divisible by 3 */
  CODON_GAIN,             /* Insertion length is divisible by 3 */
  CODON_LOSS,             /* Deletion length is divisible by 3 */
  CODON_REGION,           /* Just say the variant is in the Coding Region, used in Structrual Varition*/
  INSERTION,
  DELETION,
  NONSYNONYMOUS,
  SYNONYMOUS,
  ESSENTIAL_SPLICE_SITE,
  NORMAL_SPLICE_SITE,
  UTR5,
  UTR3,
  EXON,
  INTRON,
  UPSTREAM,
  DOWNSTREAM,
  SNV,                    /*SNV contains the following 6 types, it appears when there is no reference.*/
  NONCODING,
  INTERGENIC,
  MONOMORPHIC
} AnnotationType;


#endif /* _TYPE_H_ */
