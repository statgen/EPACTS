/*************************************************************
 *
 * Moving window project 
 * File: bed_reader.h	
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   Deal with *.bed file
 *
 **************************************************************/
#ifndef _BED_READER_H        
#define _BED_READER_H 

#include <fstream>  
#include <iostream> 
#include "setid_bim_index.h"

using namespace std;

//===============================================================
//===============================================================

#define MY_CHAR_BIT 8


class BedFileReader
{
public:
	BedFileReader(char* f, char* m, char* fm, char* o, Hasht* ht, int* myerror, char* info = NULL);
	BedFileReader(char* filename, char* bim_file, 
							 char* fam_file, char* out_file, 
							 int win_size, int ovrlp_size, int encode_output, int* myerror, char* info = NULL);
	~BedFileReader();
	int m_set_counter; 
	int m_num_of_snps_insetid;

private:

	char* m_filename; 
	char* m_file_temp_name;
	char* m_filename_mwo;
	char* m_info_file;
	char* m_filename_bim;
	char* m_filename_fam;

	std::ifstream m_file;
	std::ifstream m_bim;
	std::ifstream m_fam;
	std::ofstream m_file_temp;
	std::ofstream m_info;
	std::ifstream m_infoi;

	std::fstream m_file_mwo;


	char* m_info_rewritten;
	std::ofstream m_info_rewr;
	long m_begin4rw; 

	int m_encode_output;

	int m_approx_line_lenght;//number of snps
	SNP_info* m_snp_sets;
	int m_line_counter; // 0 .... NumOfIndividuals
	int m_setnumber;
	int m_size_of_esi;
	int m_win_size;
	int m_ovrlp_size;

	char str[1000]; 
	char str2[1000];
	char str3[1000];

	void init(char* bim_file, char* fam_file, int* myerror);
	void read_data_and_update_temp_file(int* myerror);
	void read_data_and_create_mwo_used_hashtable(Hasht* ht, int* myerror);
	void upload_snpid_from_bim(int* myerror);
	void decode_byte(int* bits_val,int* individuals_counter, int* temp_snp_info0, int* temp_snp_info1,int snp_set_ind);
	void encode(int* temp_snp_info,char* encoded_snp_info);
};

#endif //_BED_READER_H

