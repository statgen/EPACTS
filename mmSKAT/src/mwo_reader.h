/*************************************************************
 *
 * Moving window project 
 * File: mwo_reader.h	
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   Deal with *.mwa file
 *
 **************************************************************/
#ifndef _MWO_READER_H        
#define _MWO_READER_H 

#include <fstream>  
#include <iostream> 

#include <algorithm>
#include <vector>

#include "DArray.h"
#include "bed_reader.h"

using namespace std;


//===============================================================
//===============================================================

class snp
{
private:
public:	
	char m_name[SNP_ID_SIZE];
	CDArray<char> m_char;

	snp(){memset(m_name,'\0',sizeof(m_name));}
	~snp(){}


};
//===============================================================
//===============================================================

class snpset   // info for column i  the columns i = 0 ... get_noc()-1
{
private:
public:	
	CDArray<snp> m_snp;      //row j  //different types of categorical. index in m_char will maped into m_integer 
	snpset(){}
	~snpset(){}

};
//===============================================================
//===============================================================

class MwoFileReader
{
public:

	MwoFileReader(char* filename, int* myerror, char* info = NULL);
	~MwoFileReader();
	void prepare_out_array_print_snpset_to_file(snpset* ss, int n, int* Z,int size, int Is_MakeFile, int* myerror, char * SNPID=NULL);
	void upload_offsets_table();
	void get_set(int set_num, int* Z,int size, int* myerror,int Is_MakeFile = 1 , char * SNPID=NULL );

	//======================================================
	//This function returns Total number of SNP Sets in current ".mwa"
	//RETURN VALUE: Total Number of SNP Sets in the MWA file. 
	//======================================================
	int get_num_of_sets(){return this->m_total_num_of_sets;}

	//======================================================
	//This function returns Total number of SNP Sets in current ".mwa"
	//Total_Num_SNPSets : Total Number of SNP Sets in the MWA file. 
	//======================================================
	void get_TotalNumberofSets(int *Total_Num_SNPSets){ *Total_Num_SNPSets = this->m_total_num_of_sets; }
	
	//======================================================
	//This function returns Total number of SNP that in work now
	//Total_Num_SNP : Total Number of SNPs in the MWA file. 
	//======================================================
	void get_TotalNumberofSNPs(int *Total_Num_SNP){ *Total_Num_SNP = this->m_total_num_of_snps; }

	//======================================================
	//This function returns Total number of IND
	//Total_Num_SNP : Total Number of individuals in the MWA file. 
	//======================================================
	void get_TotalNumberofInd(int *Total_Num_IND) { *Total_Num_IND = this->m_num_of_individuals;}
	
		
	//======================================================
	//Return a number of SNPs in the given SNP Set. 
	//SetID : SetID
	//Num_SNP : number of SNPs in the given SNP Set. 
	//======================================================
	void get_NumberofSnps(int SetID,int *Num_SNP,int* myerror);
	
	


private:
	char* m_filename; 
	std::ifstream m_file;
	std::ifstream m_infoin;
	CDArray<snpset> m_snpsets;
	
	int m_win_size;
	int m_num_of_bytes_per_line;
	int m_ovlp_size ;
	int m_num_of_different_snps; // num of snps without repetition - like in BIM file
	int m_total_num_of_snps;	 // num of snps with repetition like in BUM file + overlaping
	int m_total_num_of_sets;
	int m_num_of_individuals;
	long* m_offsetarr;
	int* m_set_size;

	void decode_byte(int* bits_val,char* buff, int* ind_count);
	void Tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters = " ");
};

#endif //_MWO_READER_H


