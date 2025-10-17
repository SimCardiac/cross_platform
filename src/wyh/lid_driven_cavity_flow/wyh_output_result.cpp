#include "wyh_output_result.h"
#include <fstream>
#include <iostream>
#include <stdio.h> 
#include <string>

using namespace std;

//=======================================================================
void wyh_output_result::output_result1D(double *result, int mesh_size, char *result_name_char, int judge_num)
{
	int o_num = mesh_size;

	char text_result_name[500];
	sprintf(text_result_name, "%s%d%s", result_name_char, judge_num, ".txt");
	ofstream outfile_result;
	outfile_result.open(text_result_name);

	for (int k = 0; k < o_num; k++){
		outfile_result << result[k] << endl;
	}
	outfile_result.close();
}

//=======================================================================
void wyh_output_result::output_result2D(double *result, int *mesh_size, char *result_name_char, int judge_num)
{
	int m_num = mesh_size[0];
	int n_num = mesh_size[1];

	char text_result_name[500];
	sprintf(text_result_name, "%s%d%s", result_name_char, judge_num, ".txt");
	ofstream outfile_result;
	outfile_result.open(text_result_name);

	for (int i = 0; i < n_num; i++){
		for (int j = 0; j < m_num; j++){
			outfile_result << result[i*m_num + j] << "\t";
		}
		outfile_result << endl;
	}
	outfile_result.close();
}

//=======================================================================
void wyh_output_result::output_result3D(double *result, int *mesh_size, char *result_name_char, int judge_num)
{
	int m_num = mesh_size[0];
	int n_num = mesh_size[1];
	int o_num = mesh_size[2];

	char text_result_name[500];
	sprintf(text_result_name, "%s%d%s", result_name_char, judge_num, ".txt");
	ofstream outfile_result;
	outfile_result.open(text_result_name);

	for (int k = 0; k < o_num; k++){
		for (int i = 0; i < n_num; i++){
			for (int j = 0; j < m_num; j++){
				outfile_result << result[k*n_num*m_num + i*m_num + j] << endl;
			}
		}
	}
	outfile_result.close();
}