// movingwindow.cpp : main project file.
//SBM = SETID BEDFILE MERGE

#include "bed_reader.h"
#include "mwo_reader.h"
#include "setid_bim_index.h"
#include "interface_to_R.h"

int main(int argc, char *argv[] )
{
	//=========  WORK WITH R INTERFACE   ===========

	int Total_Num_SNPSets;
	int Total_Num_SNP;
	int Total_Num_IND;
	int my_error;



	/*
	Generate_MWA_SetID_File((char*)"../data/Example1.bed",(char*)"../data/Example1.bim",(char*)"../data/Example1.fam",
		(char*)"../data/Example1.SetID", (char*)"../data/Example1_0124_SetID.mwa", (char*)"../data/Example1_0124_INFO.txt",&my_error);

	*/
	/*
	Generate_MWA_SetID_File((char*)"../data/reiner.bed",(char*)"../data/reiner.bim",(char*)"../data/reiner.fam",
		(char*)"../data/reiner.SetID", (char*)"../data/reiner_0124_SetID.mwa", (char*)"../data/reiner_0124_INFO.txt",&my_error);

	if(my_error != 0)
		return 0;
	Kill_MWA_SetID_File();
    */


	/* Open_MWA((char*)"../data/reiner_0124_SetID.mwa", (char*)"../data/reiner_0124_INFO.txt",&my_error); */
	Open_MWA((char*)"/home/seunggeun/Project/Kernel_Machine/Test_From_Others/Yasin/test.SSD"
	, (char*)"/home/seunggeun/Project/Kernel_Machine/Test_From_Others/Yasin/test.SSD.info",&my_error);

	if(my_error != 0)
		return 0;
	Get_TotalNumberofSets(&Total_Num_SNPSets);
	Get_TotalNumberofSNPs(&Total_Num_SNP);
	Get_TotalNumberofInd( &Total_Num_IND);

	int Num_SNP7;
	int* TheGenotypesArray7;
	int size;
	int i;
	for (i = 1; i < Total_Num_SNPSets+1; ++ i)
	{
		Get_NumberofSnps(i,&Num_SNP7, &my_error);
		if(my_error != 0) continue;
		TheGenotypesArray7 = new int [Total_Num_IND*Num_SNP7];
		size = Total_Num_IND*Num_SNP7;
		Get_Genotypes(i, TheGenotypesArray7,size, 1,&my_error); // int Is_MakeFile = 1 -> print set to file, 0 -> don't
		Get_Genotypes(i, TheGenotypesArray7,size, 1,&my_error);
		//delete [] TheGenotypesArray7;
		if(my_error != 0) continue;
	}
/*
		Get_NumberofSnps(1,&Num_SNP7, &my_error);
		//if(my_error != 0) continue;
		TheGenotypesArray7 = new int [Total_Num_IND*Num_SNP7];
		size = Total_Num_IND*Num_SNP7;
		Get_Genotypes(1, TheGenotypesArray7,size, 1,&my_error); // int Is_MakeFile = 1 -> print set to file, 0 -> don't
		delete [] TheGenotypesArray7;
		//if(my_error != 0) continue;



	for (i = 1; i < 11; ++ i)
	{
		Get_NumberofSnps(i,&Num_SNP7, &my_error);
		if(my_error != 0) continue;
		TheGenotypesArray7 = new int [Total_Num_IND*Num_SNP7];
		size = Total_Num_IND*Num_SNP7;
		Get_Genotypes(i, TheGenotypesArray7,size, 1,&my_error); // int Is_MakeFile = 1 -> print set to file, 0 -> don't
		delete [] TheGenotypesArray7;
		if(my_error != 0) continue;


	}
*/
	Close_MWA();
	return 1;
	//============================================================
	Generate_MWA_MovingWindow( (char*)"../data/CEU.bed",(char*)"../data/CEU.bim",(char*)"../data/CEU.fam",
							(char*)"../data/CEU_MW.mwa", 20, 5, (char*)"../data/MW_INFO.txt",&my_error );
	Kill_MWA_MovingWindow();
	Open_MWA((char*)"../data/CEU_MW.mwa",(char*)"../data/MW_INFO.txt",&my_error);
	Get_TotalNumberofSets(&Total_Num_SNPSets);
	Get_TotalNumberofSNPs(&Total_Num_SNP);
	Get_TotalNumberofInd( &Total_Num_IND);

	int Num_SNP25;
	Get_NumberofSnps(25,&Num_SNP25,&my_error);
	int* TheGenotypesArray25 = new int [Total_Num_IND*Num_SNP25];

	Get_Genotypes(25, TheGenotypesArray25,Total_Num_IND*Num_SNP25, 1,&my_error);
	Close_MWA();


	//============================================================
	//Generate_MWA_SetID_File((char*)"../data/CEU.bed",(char*)"../data/CEU.bim",(char*)"../data/CEU.fam"
	//	, (char*)"../data/SetID.TXT", (char*)"../data/CEU_SetID.mwa", (char*)"../data/SetID_INFO.txt",&my_error);
	//if(my_error != 0) return 0;
	//Kill_MWA_SetID_File();
	/*
	Open_MWA((char*)"../data/CEU_SetID.mwa", (char*)"../data/SetID_INFO.txt",&my_error);
	if(my_error != 0) return 0;
	Get_TotalNumberofSets(&Total_Num_SNPSets);
	Get_TotalNumberofSNPs(&Total_Num_SNP);
	Get_TotalNumberofInd( &Total_Num_IND);
	int Num_SNP8;
	int* TheGenotypesArray8;
	for (int i = 1; i < 9; ++i)
	{
		Get_NumberofSnps(i,&Num_SNP8, &my_error);
		if(my_error != 0) continue;

		TheGenotypesArray8 = new int[Total_Num_IND*Num_SNP8];
		int size = 5; //Total_Num_IND*Num_SNP8;
		Get_Genotypes(i, TheGenotypesArray8,size, 1,&my_error); // int Is_MakeFile = 1 -> print set to file, 0 -> don't
		delete [] TheGenotypesArray8;
		if(my_error != 0) continue;

	}
	Close_MWA();
	*/

   //============================================

	//Generate_MWA_MovingWindow( (char*)"../data/CEU.bed",(char*)"../data/CEU.bim",(char*)"../data/CEU.fam",
	//						(char*)"../data/CEU_MW.mwa", 20, 5, (char*)"../data/MW_INFO.txt",&my_error );
	Open_MWA((char*)"../data/CEU_MW.mwa",(char*)"../data/MW_INFO.txt",&my_error);
	Get_TotalNumberofSets(&Total_Num_SNPSets);
	Get_TotalNumberofSNPs(&Total_Num_SNP);
	Get_TotalNumberofInd( &Total_Num_IND);
	int Num_SNP9;
	Get_NumberofSnps(4,&Num_SNP9,&my_error);
	int* TheGenotypesArray9 = new int [Total_Num_IND*Num_SNP9];

	Get_Genotypes(4, TheGenotypesArray9,Total_Num_IND*Num_SNP9, 1,&my_error);
	Close_MWA();
	//Kill_MWA_MovingWindow();

	/*
		//int* TheGenotypesArray = NULL;
	//Get_Genotypes(-33, TheGenotypesArray);
	//Get_Genotypes(33, TheGenotypesArray);
	//Get_Genotypes(51, TheGenotypesArray);
	//Get_Genotypes(30, TheGenotypesArray);

	int Num_SNP1, Num_SNP2, Num_SNP3,Num_SNP4,Num_SNP5;
	Get_NumberofSnps(33,&Num_SNP1);
	Get_NumberofSnps(-14,&Num_SNP2);
	Get_NumberofSnps(5000,&Num_SNP3);
	Get_NumberofSnps(2,&Num_SNP4);
	Get_NumberofSnps(1,&Num_SNP5);
*/
	//=========  END OF WORK WITH R INTERFACE   ===========

	//======================================================================================================

	//SBM = SETID BEDFILE MERGE
/*	Hasht hash_table("../data/SetID.TXT", "../data/CEU.bim");
	BedFileReader bed_pedigree( "../data/CEU.bed","../data/CEU.bim","../data/CEU.fam","../data/CEU.mwa",&hash_table);
	MwoFileReader mwa("../data/CEU.mwa", NULL, &bed_pedigree); //it will take also "CEU.bed.INFO.txt"
	mwa.get_set(1); //some integer - enter the value base on CEU.bed.INFO.txt
	mwa.get_set(2); //some integer - enter the value base on CEU.bed.INFO.txt
	mwa.get_set(3); //some integer - enter the value base on CEU.bed.INFO.txt
*/
//	int a;
//	mwa.get_TotalNumberofSets(&a);


//

}

//=============================================================
//=============================================================
//=============================================================
//=============================================================
//=============================================================

/*
void print_usage()
{
		std::cout << "Full usage: " <<  std::endl;
		std::cout << "\t> mwa -p bedfile -m bimfile"<<  std::endl;
		std::cout << "\t      -s setfile -o outputfile" <<  std::endl;
		std::cout << "\t> mwa -p bedfile -s setfile -set set#" <<  std::endl;
		std::cout << "Minimal usage: " <<  std::endl;
		std::cout << "\t> mwa -p bedfile -s setfile" <<  std::endl <<  std::endl;

		std::cout << "Defaults:" << std::endl;

		std::cout << "\tThe default value for: set# = 1;" << std::endl;
		std::cout << "\tThe default value for: mapfile = bedfile[:-4]+.bim ;" << std::endl;
		std::cout << "\tThe default value for: outputfile = bedfile[:-4]+.mwa ;" << std::endl << std::endl;


		std::cout << "Notes:" << std::endl;
		std::cout << "\t\"-set\" - option should be one before last word,\n" <<
				     "\t\"set#\" - an integer should be last word" << std::endl;
		std::cout << "\tIf you use BED file, BIM and FAM should be in same directory -\n" <<
					 "\tthey will taken automatically from there." << std::endl;
		std::cout << "\tMWA file - the output will be created in same directory." << std::endl;
}

int main(int argc, char *argv[] )
{
	//int win_size = -999;
	//int ovrlp_size = -999;
	char* out_file = NULL;
	char* ped_file = NULL;
	char* map_file = NULL;
	char* fam_file = NULL;
	if (argc == 2 && (strcmp(argv[1], "-h")== 0 || strcmp(argv[1], "-help")== 0))
	{
		print_usage();
		exit(0);
	}

	if((argc < 3) ||
	   (argc == 3 && strcmp(argv[1], "-p") != 0) ||
	   (argc == 5 && strcmp(argv[1], "-p") != 0 && strcmp(argv[3], "-p") != 0) ||
	   (argc == 7 && strcmp(argv[1], "-p") != 0 && strcmp(argv[3], "-p") != 0 && strcmp(argv[5], "-p") != 0) ||
	   (argc == 9 && strcmp(argv[1], "-p") != 0 && strcmp(argv[3], "-p") != 0 && strcmp(argv[5], "-p") != 0 && strcmp(argv[7], "-p") != 0) ||
	   (argc == 11 && strcmp(argv[1], "-p") != 0 && strcmp(argv[3], "-p") != 0 && strcmp(argv[5], "-p") != 0 && strcmp(argv[7], "-p") != 0 && strcmp(argv[9], "-p") != 0))

	{
		print_usage();
		exit(0);
	}

	for (int i = 0; i < 6; ++ i )
	{
		char* paramflag;
		try
		{
			paramflag = argv[i * 2 + 1];
			if (strcmp(paramflag, "-p") == 0 && ped_file == NULL)
				ped_file = argv[i * 2 + 2];
			else if( strcmp(paramflag,"-m")== 0 && map_file == NULL)
				map_file = argv[i * 2 + 2];
			//else if( strcmp(paramflag,"-s")== 0 && win_size == -999)
			//	win_size = atoi(argv[i * 2 + 2]);
			//else if( strcmp(paramflag,"-l")== 0 && ovrlp_size == -999)
			//	ovrlp_size = atoi(argv[i * 2 + 2]);
			else if( strcmp(paramflag,"-o")== 0 && out_file == NULL)
				out_file = argv[i * 2 + 2];
			else if( strcmp(paramflag,"-set")== 0 )
				continue;
			else
			{
				print_usage();
				exit(0);
			}
		}
		catch(...)
		{
			if (ped_file == NULL)
			{
				std::cout << "Check your command line parameters." << std::endl;
				print_usage();
				exit(0);
			}
			if (out_file == NULL)
			{
				int l = strlen(ped_file);
				out_file = new char [l+1];
				out_file[0] = '\0';
				strncat( out_file ,ped_file,l-3 );
				strcat(out_file ,"mwa");
				out_file[l] = '\0';
			}
			if(map_file == NULL)
			{
				int l = strlen(ped_file);
				map_file = new char [l+1];
				map_file[0] = '\0';
				strncat( map_file ,ped_file,l-3 );
				if (ped_file[l-3] == 'p' && ped_file[l-2] == 'e' && ped_file[l-1] == 'd')
					strcat(map_file ,"map");
				else if (ped_file[l-3] == 'b' && ped_file[l-2] == 'e' && ped_file[l-1] == 'd')
					strcat(map_file ,"bim");
				map_file[l] = '\0';
			}
			//if(win_size == -999)
			//	win_size = 500;
			//if(ovrlp_size == -999)
			//	ovrlp_size = win_size * 0.1; //10% from 500

			break;
		}
	}

	int l = strlen(ped_file);
	if (ped_file[l-3] == 'b' && ped_file[l-2] == 'e' && ped_file[l-1] == 'd')
	{
		fam_file = new char [l+1];
		fam_file[0] = '\0';
		strncat( fam_file ,ped_file,l-3 );
		strcat( fam_file ,"fam");
		fam_file[l] = '\0';

	}

	int encode_output = 1;


	//check if output already exists
	//check if info already exists
	char str2[1000];
	memset(str2,'\0',sizeof(str2));
	strcpy (str2,ped_file);
	strcat (str2,".INFO.txt");

	std::ifstream ifile(out_file);
	std::ifstream infofile(str2);
	if (!ifile || !infofile)
	{
		//if one of above not exists run the following
		//==========================================
		//if ".bed" read the follow, and use "*.fam" - extension to "*.bed"  ,"*.bim" - likes ".map"
		//==========================================
		BedFileReader bed_pedigree(ped_file, map_file, fam_file, out_file, win_size, ovrlp_size, encode_output );
	}


	//TODO check if INFO exists, if not - don't run the foillowing
	MwoFileReader mwo(out_file);
	if (strcmp(argv[argc-2], "-set") == 0)
	{
		if (atoi(argv[argc-1]) > mwo.get_num_of_sets() || atoi(argv[argc-1]) < 0 )
		{
			std::cout << "Your set number is incorrect. " << std::endl;
			std::cout << "Total number of sets is: " << mwo.get_num_of_sets() << std::endl;
			std::cout << "Set number should be >= then 1 and <= " << mwo.get_num_of_sets() << std::endl;
			return 0;
		}
		std::cout << "The total number of sets is: " << mwo.get_num_of_sets() << std::endl;
		mwo.get_set(atoi(argv[argc-1])); // here it's also writing to the file
	}
	else
		mwo.get_set(1);
	std::cout << "The set number " << argv[argc-1] << " was printed for you." << std::endl;

    return 0;
}

*/

