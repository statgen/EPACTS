/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __CODONHELPER_H__
#define __CODONHELPER_H__

#include <map>
#include <vector>
#include <string>
#include <stdint.h>

//#include <cstdio>

class codonHelper
{
 public:
  std::map<std::string,char> cAMap; // 3-letters to AA char
  std::vector<char> nAMap;  // 0-63 to AA char
  std::map<std::string,std::vector<uint8_t> > cDMap; // 3 letters to degeneracy
  std::vector<uint8_t> nDMap;  // codon (0-63) * npos (0-2) to degeneracy
  std::vector<std::string> n2c;
  std::map<std::string,uint8_t> c2n;

  char b2n[5]; 
  int n2b[256];

  char getVarAA(const char* codon, int pos, char varnt) {
    std::string c(codon);
    c[pos] = varnt;
    return cAMap[c];
  }

  uint8_t getDegeneracy(char aa, char aaA, char aaC, char aaG, char aaT) {
    uint8_t ret = 0;
    if ( aa == aaA ) ++ret;
    if ( aa == aaC ) ++ret;
    if ( aa == aaG ) ++ret;
    if ( aa == aaT ) ++ret;
    return ret;
  }

  codonHelper() {
    // build string codon table
    cAMap["TTT"] = 'F';
    cAMap["TTC"] = 'F';
    cAMap["TCT"] = 'S';
    cAMap["TCC"] = 'S';
    cAMap["TAT"] = 'Y';
    cAMap["TAC"] = 'Y';
    cAMap["TGT"] = 'C';
    cAMap["TGC"] = 'C';
    cAMap["TTA"] = 'L';
    cAMap["TCA"] = 'S';
    cAMap["TAA"] = '*';
    cAMap["TGA"] = '*';
    cAMap["TTG"] = 'L';
    cAMap["TCG"] = 'S';
    cAMap["TAG"] = '*';
    cAMap["TGG"] = 'W';
    cAMap["CTT"] = 'L';
    cAMap["CTC"] = 'L';
    cAMap["CCT"] = 'P';
    cAMap["CCC"] = 'P';
    cAMap["CAT"] = 'H';
    cAMap["CAC"] = 'H';
    cAMap["CGT"] = 'R';
    cAMap["CGC"] = 'R';
    cAMap["CTA"] = 'L';
    cAMap["CTG"] = 'L';
    cAMap["CCA"] = 'P';
    cAMap["CCG"] = 'P';
    cAMap["CAA"] = 'Q';
    cAMap["CAG"] = 'Q';
    cAMap["CGA"] = 'R';
    cAMap["CGG"] = 'R';
    cAMap["ATT"] = 'I';
    cAMap["ATC"] = 'I';
    cAMap["ACT"] = 'T';
    cAMap["ACC"] = 'T';
    cAMap["AAT"] = 'N';
    cAMap["AAC"] = 'N';
    cAMap["AGT"] = 'S';
    cAMap["AGC"] = 'S';
    cAMap["ATA"] = 'I';
    cAMap["ACA"] = 'T';
    cAMap["AAA"] = 'K';
    cAMap["AGA"] = 'R';
    cAMap["ATG"] = 'M';
    cAMap["ACG"] = 'T';
    cAMap["AAG"] = 'K';
    cAMap["AGG"] = 'R';
    cAMap["GTT"] = 'V';
    cAMap["GTC"] = 'V';
    cAMap["GCT"] = 'A';
    cAMap["GCC"] = 'A';
    cAMap["GAT"] = 'D';
    cAMap["GAC"] = 'D';
    cAMap["GGT"] = 'G';
    cAMap["GGC"] = 'G';
    cAMap["GTA"] = 'V';
    cAMap["GTG"] = 'V';
    cAMap["GCA"] = 'A';
    cAMap["GCG"] = 'A';
    cAMap["GAA"] = 'E';
    cAMap["GAG"] = 'E';
    cAMap["GGA"] = 'G';
    cAMap["GGG"] = 'G';    

    b2n[0] = 'A'; b2n[1] = 'C'; b2n[2] = 'G'; b2n[3] = 'T'; b2n[4] = 'N';
    for(int i=0; i < 256; ++i) { 
      n2b[i] = 4;
    }
    for(int i=0; i < 4; ++i) { 
      n2b[(int)b2n[i]] = i; 
    }
    char codon[4] = {'N','N','N','\0'};

    nAMap.resize(64,'\0');
    nDMap.resize(64*3,'\0');
    n2c.resize(64);
    for(int i=0; i < 4; ++i) {
      codon[0] = b2n[i];
      for(int j=0; j < 4; ++j) {
	codon[1] = b2n[j];
	for(int k=0; k < 4; ++k) {
	  codon[2] = b2n[k];
	  int nCodon = i*16+j*4+k;
	  nAMap[nCodon] = cAMap[codon]; // build integer codon table
	  n2c[nCodon] = codon;          // convert integer codons to string
	  c2n[codon] = nCodon;          // convert string codons to integer
	}
      }
    }

    for(int i=0; i < 4; ++i) {
      codon[0] = b2n[i];
      for(int j=0; j < 4; ++j) {
	codon[1] = b2n[j];
	for(int k=0; k < 4; ++k) {
	  codon[2] = b2n[k];
	  int nCodon = i*16+j*4+k;
	  cDMap[codon].resize(3);       // build degeneracy table
	  // check the first position;
	  cDMap[codon][0] = nDMap[nCodon*3+0] = getDegeneracy(nAMap[nCodon],nAMap[0*16+j*4+k],nAMap[1*16+j*4+k],nAMap[2*16+j*4+k],nAMap[3*16+j*4+k]);
	  cDMap[codon][1] = nDMap[nCodon*3+1] = getDegeneracy(nAMap[nCodon],nAMap[i*16+0*4+k],nAMap[i*16+1*4+k],nAMap[i*16+2*4+k],nAMap[i*16+3*4+k]);
	  cDMap[codon][2] = nDMap[nCodon*3+2] = getDegeneracy(nAMap[nCodon],nAMap[i*16+j*4+0],nAMap[i*16+j*4+1],nAMap[i*16+j*4+2],nAMap[i*16+j*4+3]);
	}
      }
    }
  }
};

extern codonHelper codonHelp;

#endif
