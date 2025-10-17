#ifndef WYH_output_result_H
#define WYH_output_result_H

class wyh_output_result
{
private:
	;
public:
	void output_result1D(double *result, int mesh_size, char *result_name, int judge_num);
	void output_result2D(double *result, int *mesh_size, char *result_name, int judge_num);
	void output_result3D(double *result, int *mesh_size, char *result_name, int judge_num);
};


#endif