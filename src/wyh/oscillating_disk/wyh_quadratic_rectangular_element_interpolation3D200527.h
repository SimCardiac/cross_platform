#ifndef wyh_quadratic_rectangular_element_interpolation3D200527H
#define wyh_quadratic_rectangular_element_interpolation3D200527H

#include <petsc.h>
#include <mpi.h>

class wyh_quadratic_rectangular_element_interpolation3D200527
{
public:
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int root = 0;
	double pi = 3.141592653;

	PetscInt istart;
	PetscInt iend;

	//与实际问题相关+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//变量---------------------------------------
	
	int m_num = 32;
	int n_num = 32;
	int o_num = 32;
	int all_mesh_num = m_num*n_num*o_num;

	int Lagrange_node_num = 4385;

	//空间
	double x_L = 0.0;	double x_U = 1.0;
	double y_L = 0.0;	double y_U = 1.0;
	double z_L = 0.0;	double z_U = 1.0;
	
	double hx = (x_U - x_L) / m_num;
	double hy = (y_U - y_L) / n_num;
	double hz = (z_U - z_L) / o_num;

	//时间
	double dt = 0.125*hx;

	//坐标
	Vec Lagrange_node_x;
	Vec Lagrange_node_y;
	Vec Lagrange_node_z;

	Vec Lagrange_velocity_u;
	Vec Lagrange_velocity_v;
	Vec Lagrange_velocity_w;
	
	Mat Double_quadratic_rectangular_element_matrix3D;

	//函数---------------------------------------

	//更新位置
	void update_Lagrange_node3D(MPI_Comm PETSC_COMM_WORLD,
		Vec velocity_u_at_p, Vec velocity_v_at_p, Vec velocity_w_at_p);

	//组装矩阵
	void assembly_rectangular_element_matrix3D(Mat Double_quadratic_rectangular_element_matrix3D);

	//设置变量大小
	void set_variable_size(MPI_Comm PETSC_COMM_WORLD);

	//Destroy
	void destroy_variable();

};
#endif
