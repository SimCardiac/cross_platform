#include "wyh_quadratic_rectangular_element_interpolation3D200527.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include <petsc.h>
#include <mpi.h>

using namespace std;

//设置变量大小-----------------------------------
void wyh_quadratic_rectangular_element_interpolation3D200527::set_variable_size(MPI_Comm PETSC_COMM_WORLD)
{
	//拉格朗日坐标
	VecCreate(PETSC_COMM_WORLD, &Lagrange_node_x);
	VecSetSizes(Lagrange_node_x, PETSC_DECIDE, Lagrange_node_num);
	VecSetFromOptions(Lagrange_node_x);
	VecDuplicate(Lagrange_node_x, &Lagrange_node_y);
	VecDuplicate(Lagrange_node_x, &Lagrange_node_z);
	VecDuplicate(Lagrange_node_x, &Lagrange_velocity_u);
	VecDuplicate(Lagrange_node_x, &Lagrange_velocity_v);
	VecDuplicate(Lagrange_node_x, &Lagrange_velocity_w);
	
	MatCreate(PETSC_COMM_WORLD, &Double_quadratic_rectangular_element_matrix3D);
	MatSetSizes(Double_quadratic_rectangular_element_matrix3D, PETSC_DECIDE, PETSC_DECIDE, Lagrange_node_num, all_mesh_num);
	MatSetFromOptions(Double_quadratic_rectangular_element_matrix3D);
}

//-----------------------------------------------
void wyh_quadratic_rectangular_element_interpolation3D200527::update_Lagrange_node3D(MPI_Comm PETSC_COMM_WORLD,
	Vec velocity_u_at_p, Vec velocity_v_at_p, Vec velocity_w_at_p)
{
	//组装矩阵
	//assembly_Double_quadratic_rectangular_element_matrix3D(Double_quadratic_rectangular_element_matrix3D);
	assembly_rectangular_element_matrix3D(Double_quadratic_rectangular_element_matrix3D);

	//求解拉格朗日点速度
	MatMult(Double_quadratic_rectangular_element_matrix3D, velocity_u_at_p, Lagrange_velocity_u);
	MatMult(Double_quadratic_rectangular_element_matrix3D, velocity_v_at_p, Lagrange_velocity_v);
	MatMult(Double_quadratic_rectangular_element_matrix3D, velocity_w_at_p, Lagrange_velocity_w);
	
	//更新位置
	VecAXPY(Lagrange_node_x, dt, Lagrange_velocity_u);
	VecAXPY(Lagrange_node_y, dt, Lagrange_velocity_v);
	VecAXPY(Lagrange_node_z, dt, Lagrange_velocity_w);
}

//组装拉格朗日插值矩阵Euler to Lagrange
void wyh_quadratic_rectangular_element_interpolation3D200527::assembly_rectangular_element_matrix3D(Mat Double_quadratic_rectangular_element_matrix3D)
{
	//边界
	double real_x_L = x_L + hx / 2.0;
	double real_y_L = y_L + hy / 2.0;
	double real_z_L = z_L + hz / 2.0;
	double half_hx = hx / 2.0;
	double half_hy = hy / 2.0;
	double half_hz = hz / 2.0;

	//
	double *Lagrange_node_x_temp;
	double *Lagrange_node_y_temp;
	double *Lagrange_node_z_temp;
	VecGetArray(Lagrange_node_x, &Lagrange_node_x_temp);
	VecGetArray(Lagrange_node_y, &Lagrange_node_y_temp);
	VecGetArray(Lagrange_node_z, &Lagrange_node_z_temp);

	//组装
	MatMPIAIJSetPreallocation(Double_quadratic_rectangular_element_matrix3D, 8, NULL, 8, NULL);
	MatSeqAIJSetPreallocation(Double_quadratic_rectangular_element_matrix3D, 8, NULL);
	MatZeroEntries(Double_quadratic_rectangular_element_matrix3D);

	MatGetOwnershipRange(Double_quadratic_rectangular_element_matrix3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		double node_x_temp = Lagrange_node_x_temp[Ii - istart];
		double node_y_temp = Lagrange_node_y_temp[Ii - istart];
		double node_z_temp = Lagrange_node_z_temp[Ii - istart];

		int i = floor((node_y_temp - real_y_L) / hy);
		int j = floor((node_x_temp - real_x_L) / hx);
		int k = floor((node_z_temp - real_z_L) / hz);

		double real_elem_x_L = real_x_L + j*hx;
		double real_elem_y_L = real_y_L + i*hy;
		double real_elem_z_L = real_z_L + k*hz;
		
		double x0 = real_elem_x_L + hx / 2.0;
		double y0 = real_elem_y_L + hy / 2.0;
		double z0 = real_elem_z_L + hz / 2.0;

		double node_x_redu_x0 = node_x_temp - x0;
		double node_y_redu_y0 = node_y_temp - y0;
		double node_z_redu_z0 = node_z_temp - z0;

		double insert_value_1 = -1.0*(node_x_redu_x0 + half_hx)*(node_y_redu_y0 + half_hy)*(node_z_redu_z0 - half_hz) / (hx*hy*hz);
		double insert_value_2 =  1.0*(node_x_redu_x0 - half_hx)*(node_y_redu_y0 + half_hy)*(node_z_redu_z0 - half_hz) / (hx*hy*hz);
		double insert_value_3 = -1.0*(node_x_redu_x0 - half_hx)*(node_y_redu_y0 - half_hy)*(node_z_redu_z0 - half_hz) / (hx*hy*hz);
		double insert_value_4 =  1.0*(node_x_redu_x0 + half_hx)*(node_y_redu_y0 - half_hy)*(node_z_redu_z0 - half_hz) / (hx*hy*hz);

		double insert_value_5 =  1.0*(node_x_redu_x0 + half_hx)*(node_y_redu_y0 + half_hy)*(node_z_redu_z0 + half_hz) / (hx*hy*hz);
		double insert_value_6 = -1.0*(node_x_redu_x0 - half_hx)*(node_y_redu_y0 + half_hy)*(node_z_redu_z0 + half_hz) / (hx*hy*hz);
		double insert_value_7 =  1.0*(node_x_redu_x0 - half_hx)*(node_y_redu_y0 - half_hy)*(node_z_redu_z0 + half_hz) / (hx*hy*hz);
		double insert_value_8 = -1.0*(node_x_redu_x0 + half_hx)*(node_y_redu_y0 - half_hy)*(node_z_redu_z0 + half_hz) / (hx*hy*hz);

		/*
		double x1 = real_elem_x_L + hx;			double y1 = real_elem_y_L + hy;			double z1 = real_elem_z_L;
		double x2 = real_elem_x_L;				double y2 = real_elem_y_L + hy;			double z2 = real_elem_z_L;
		double x3 = real_elem_x_L;				double y3 = real_elem_y_L;				double z3 = real_elem_z_L;
		double x4 = real_elem_x_L + hx;			double y4 = real_elem_y_L;				double z4 = real_elem_z_L;

		double x5 = real_elem_x_L + hx;			double y5 = real_elem_y_L + hy;			double z5 = real_elem_z_L + hz;
		double x6 = real_elem_x_L;				double y6 = real_elem_y_L + hy;			double z6 = real_elem_z_L + hz;
		double x7 = real_elem_x_L;				double y7 = real_elem_y_L;				double z7 = real_elem_z_L + hz;
		double x8 = real_elem_x_L + hx;			double y8 = real_elem_y_L;				double z8 = real_elem_z_L + hz;

		double insert_value_1 = (node_x_temp - x2)*(node_y_temp - y4)*(node_z_temp - z5) / (x1 - x2) / (y1 - y4) / (z1 - z5);
		double insert_value_2 = (node_x_temp - x1)*(node_y_temp - y4)*(node_z_temp - z5) / (x2 - x1) / (y2 - y4) / (z2 - z5);
		double insert_value_3 = (node_x_temp - x1)*(node_y_temp - y1)*(node_z_temp - z5) / (x3 - x1) / (y3 - y1) / (z3 - z5);
		double insert_value_4 = (node_x_temp - x2)*(node_y_temp - y1)*(node_z_temp - z5) / (x4 - x2) / (y4 - y1) / (z4 - z5);

		double insert_value_5 = (node_x_temp - x2)*(node_y_temp - y4)*(node_z_temp - z1) / (x5 - x2) / (y5 - y4) / (z5 - z1);
		double insert_value_6 = (node_x_temp - x1)*(node_y_temp - y4)*(node_z_temp - z1) / (x6 - x1) / (y6 - y4) / (z6 - z1);
		double insert_value_7 = (node_x_temp - x1)*(node_y_temp - y1)*(node_z_temp - z1) / (x7 - x1) / (y7 - y1) / (z7 - z1);
		double insert_value_8 = (node_x_temp - x2)*(node_y_temp - y1)*(node_z_temp - z1) / (x8 - x2) / (y8 - y1) / (z8 - z1);
		*/
		double insert_vector_value[8];

		insert_vector_value[0] = insert_value_1;	insert_vector_value[1] = insert_value_2;
		insert_vector_value[2] = insert_value_3;	insert_vector_value[3] = insert_value_4;
		insert_vector_value[4] = insert_value_5;	insert_vector_value[5] = insert_value_6;
		insert_vector_value[6] = insert_value_7;	insert_vector_value[7] = insert_value_8;

		//
		int P0_mark = (k + 0)*m_num*n_num + (i + 1)*m_num + j + 1;
		int P1_mark = (k + 0)*m_num*n_num + (i + 1)*m_num + j + 0;
		int P2_mark = (k + 0)*m_num*n_num + (i + 0)*m_num + j + 0;
		int P3_mark = (k + 0)*m_num*n_num + (i + 0)*m_num + j + 1;

		int P4_mark = (k + 1)*m_num*n_num + (i + 1)*m_num + j + 1;
		int P5_mark = (k + 1)*m_num*n_num + (i + 1)*m_num + j + 0;
		int P6_mark = (k + 1)*m_num*n_num + (i + 0)*m_num + j + 0;
		int P7_mark = (k + 1)*m_num*n_num + (i + 0)*m_num + j + 1;

		int Ij_vector[8];
		Ij_vector[0] = P0_mark;		Ij_vector[1] = P1_mark;		
		Ij_vector[2] = P2_mark;		Ij_vector[3] = P3_mark;
		Ij_vector[4] = P4_mark;		Ij_vector[5] = P5_mark;
		Ij_vector[6] = P6_mark;		Ij_vector[7] = P7_mark;

		MatSetValues(Double_quadratic_rectangular_element_matrix3D, 1, &Ii, 8, Ij_vector, insert_vector_value, INSERT_VALUES);

	}
	MatAssemblyBegin(Double_quadratic_rectangular_element_matrix3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Double_quadratic_rectangular_element_matrix3D, MAT_FINAL_ASSEMBLY);

	VecRestoreArray(Lagrange_node_x, &Lagrange_node_x_temp);
	VecRestoreArray(Lagrange_node_y, &Lagrange_node_y_temp);
	VecRestoreArray(Lagrange_node_z, &Lagrange_node_z_temp);

	//
	delete[] Lagrange_node_x_temp;
	delete[] Lagrange_node_y_temp;
	delete[] Lagrange_node_z_temp;
}

//*************************************************************************************************

//Destroy==========================================================================================

//*************************************************************************************************
void wyh_quadratic_rectangular_element_interpolation3D200527::destroy_variable()
{
	//
	VecDestroy(&Lagrange_node_x);		
	VecDestroy(&Lagrange_node_y);
	VecDestroy(&Lagrange_node_z);
	VecDestroy(&Lagrange_velocity_u);
	VecDestroy(&Lagrange_velocity_v);
	VecDestroy(&Lagrange_velocity_w);

	//释放变量
	MatDestroy(&Double_quadratic_rectangular_element_matrix3D);
}
