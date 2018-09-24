/*************************************************************
 *
 * NPTest Project
 * File: NPsort.cpp	
 * Date: January 1, 2011
 * Author: Larissa Miropolsky
 *
 * Description:
 *   sort char**, int*, double* array    
 *
 **************************************************************/


#include <iostream>
#include <cstring> 

#include "NPsort.h"


//===============================================================================
//================= SORT FUNCTION ===============================================
//===============================================================================
//Input:
//a -		void array
//b -		indexes of sorted `a` 
//size -	size of `a` and `b` 
//type -	type of `a` - D_INT, D_DOUBLE, D_CHARSTAR
//offset -	now not in use 
//			move offset from begining of every item of `a`  - by default = 0, then starting read real data:
//			read 4 bytes for type "D_INT", read 8 bytes for type "D_DOUBLE", read until get '\0' for type "D_CHARSTAR" 
//f -		flag f==0 for ascending order, f==1 for descending order
//pt2Func - function that compare two inputs and return 1,0,-1 (if first bigger then second, equal, less respectively)
//          the default functions for this parameter: compare_strings, compare_ints, compare_doubles
//===============================================================================
void sort_data::sort(const void* a, int* b, int size, DATA2SORT type, int offset, 
					 int f, int (*pt2Func)(const void*, const void *))
{
	if (type == D_CHARSTAR)
	{
		char** arr = ((char **)a); 
		sort(arr, b, size, f, pt2Func, offset);
	}
	else if (type == D_DOUBLE)
	{
		const double* arr = (const double *)a;
		sort(arr, b, size, f, pt2Func, offset );
	}
	else if (type == D_INT)
	{
		const int* arr = (const int *)a;
		sort(arr, b, size, f, pt2Func, offset );
	}
}
//===============================================================================
//===============================================================================

int sort_data::compare_strings(const void* x, const void* y)
{
	int result = 0;
	const char* xx = (char *)x;
	const char* yy = (char *)y;
	result = strcmp(xx, yy);
	return result;
}

int sort_data::compare_ints(const void* x, const void* y)
{
	int result = 0;
	const int* xx = (int *)x;
	const int* yy = (int *)y;
	if (*xx < *yy)
		result = -1;
	else if (*xx > *yy)
		result = 1;
	return result;
}

int sort_data::compare_doubles(const void* x, const void* y)
{
	int result = 0;
	const double* xx = (double *)x;
	const double* yy = (double *)y;
	if (*xx < *yy)
		result = -1;
	else if (*xx > *yy)
		result = 1;
	return result;
}

//===============================================================================
//===============================================================================
//INTEGER array input
//===============================================================================
//===============================================================================

void sort_data::sort(const int* a, int* b, int size, int f, int (*pt2Func)(const void*, const void *), int offset)
{
	for (int i = 0; i < size; ++ i)
		b[i] = i;

	if (pt2Func == NULL)
	{
		if (f == 0)//ascending
			sort_data::quicksort_asc(b, 0, size-1, a, &sort_data::compare_ints);
		else if (f == 1)//descending order 
			sort_data::quicksort_desc(b, 0, size-1, a, &sort_data::compare_ints);
	}
	else
	{
		if (f == 0)//ascending
			sort_data::quicksort_asc(b, 0, size-1, a, pt2Func);
		else if (f == 1)//descending order 
			sort_data::quicksort_desc(b, 0, size-1, a, pt2Func);
	}
}
//===============================================================================
//===============================================================================
//DOUBLE array input
//===============================================================================
//===============================================================================
void sort_data::sort(const double* a, int* b, int size, int f, int (*pt2Func)(const void*, const void *), int offset)
{
	for (int i = 0; i < size; ++ i)
		b[i] = i;
	if (pt2Func == NULL)
	{
		if (f == 0)//ascending
			sort_data::quicksort_asc(b, 0, size-1, a, &sort_data::compare_doubles);
		else if (f == 1)//descending order 
			sort_data::quicksort_desc(b, 0, size-1, a, &sort_data::compare_doubles);
	}
	else
	{
		if (f == 0)//ascending
			sort_data::quicksort_asc(b, 0, size-1, a, pt2Func);
		else if (f == 1)//descending order 
			sort_data::quicksort_desc(b, 0, size-1, a, pt2Func);
	}
}
//==============================================================================
//===============================================================================
//CHAR* array input
//===============================================================================
//===============================================================================
void sort_data::sort(char** a, int* b, int size, int f, int (*pt2Func)(const void*, const void *),int offset)
{
	for (int i = 0; i < size; ++ i)
		b[i] = i;
	if (pt2Func == NULL)
	{	
		if (f == 0)//ascending
			sort_data::quicksort_asc(b, 0, size-1, a, &sort_data::compare_strings);
		else if (f == 1)//descending order 
			sort_data::quicksort_desc(b, 0, size-1, a, &sort_data::compare_strings);	
	}
	else
	{	
		if (f == 0)//ascending
			sort_data::quicksort_asc(b, 0, size-1, a, pt2Func);
		else if (f == 1)//descending order 
			sort_data::quicksort_desc(b, 0, size-1, a, pt2Func);	
	}
}



//===============================================================================
//===============================================================================
//=================== Quick sort implementation =================================
//===============================================================================
//===============================================================================

//===============================================================================
//===============================================================================
//ascending order for int array
//===============================================================================
//===============================================================================

void sort_data::quicksort_asc(int* arr, int left, int right, const int* a, int (*pt2Func)(const void*, const void *)) 
{
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];
      /* partition */
      while (i <= j) {
		    while(pt2Func(&a[arr[i]], &a[pivot]) < 0) 
                  i++;
            while(pt2Func(&a[arr[j]], &a[pivot]) > 0) 
                  j--;
            if (i <= j) 
			{
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      }
      /* recursion */
      if (left < j)
            quicksort_asc(arr, left, j, a, pt2Func);
      if (i < right)
            quicksort_asc(arr, i, right, a, pt2Func);
}

//===============================================================================
//===============================================================================
//descending order for int array
//===============================================================================
//===============================================================================

void sort_data::quicksort_desc(int* arr, int left, int right, const int* a, int (*pt2Func)(const void*, const void *)) 
{
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];
      /* partition */
      while (i <= j) {
            while(pt2Func(&a[arr[i]], &a[pivot]) > 0) 
                  i++;
            while(pt2Func(&a[arr[j]], &a[pivot]) < 0) 
                  j--;
            if (i <= j) 
			{
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      }
      /* recursion */
      if (left < j)
            quicksort_desc(arr, left, j, a, pt2Func);
      if (i < right)
            quicksort_desc(arr, i, right, a, pt2Func);
}
//===============================================================================
//===============================================================================
//ascending order for double array
//===============================================================================
//===============================================================================

void sort_data::quicksort_asc(int* arr, int left, int right, const double* a, int (*pt2Func)(const void*, const void *)) 
{
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];
      /* partition */
      while (i <= j) {
		    while(pt2Func(&a[arr[i]], &a[pivot]) < 0) 
                  i++;
            while(pt2Func(&a[arr[j]], &a[pivot]) > 0) 
                  j--;
            if (i <= j) 
			{
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      }
      /* recursion */
      if (left < j)
            quicksort_asc(arr, left, j, a, pt2Func);
      if (i < right)
            quicksort_asc(arr, i, right, a, pt2Func);
}

//===============================================================================
//===============================================================================
//descending order for double array
//===============================================================================
//===============================================================================
void sort_data::quicksort_desc(int* arr, int left, int right, const double* a, int (*pt2Func)(const void*, const void *)) 
{
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];
      /* partition */
      while (i <= j) {
		    while(pt2Func(&a[arr[i]], &a[pivot]) > 0) 
                  i++;
            while(pt2Func(&a[arr[j]], &a[pivot]) < 0) 
                  j--;
            if (i <= j) 
			{
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      }
      /* recursion */
      if (left < j)
            quicksort_desc(arr, left, j, a, pt2Func);
      if (i < right)
            quicksort_desc(arr, i, right, a, pt2Func);
}

//===============================================================================
//===============================================================================
//ascending order for char* array
//===============================================================================
//===============================================================================

void sort_data::quicksort_asc(int* arr, int left, int right, char** a, int (*pt2Func)(const void*, const void *)) 
{
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];
      /* partition */
      while (i <= j) {
            while(pt2Func(a[arr[i]], a[pivot]) < 0) 
                  i++;
            while(pt2Func(a[arr[j]], a[pivot]) > 0) 
                  j--;
            if (i <= j) 
			{
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      }
      /* recursion */
      if (left < j)
            quicksort_asc(arr, left, j, a, pt2Func);
      if (i < right)
            quicksort_asc(arr, i, right, a, pt2Func);
}

//===============================================================================
//===============================================================================
//descending order for char* array
//===============================================================================
//===============================================================================
void sort_data::quicksort_desc(int* arr, int left, int right, char** a, int (*pt2Func)(const void*, const void *)) 
{
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];
      /* partition */
	  while (i <= j) {
		  while(pt2Func(a[arr[i]], a[pivot]) > 0) 
			  i++;
		  while(pt2Func(a[arr[j]], a[pivot]) < 0) 
			  j--;
		  if (i <= j) 
		  {
			  tmp = arr[i];
			  arr[i] = arr[j];
			  arr[j] = tmp;
			  i++;
			  j--;
		  }
	  }
      /* recursion */
      if (left < j)
            quicksort_desc(arr, left, j, a, pt2Func);
      if (i < right)
            quicksort_desc(arr, i, right, a, pt2Func);
}
//==================================================================

