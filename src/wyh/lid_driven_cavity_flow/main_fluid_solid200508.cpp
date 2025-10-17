static char help[] = "A structured-grid fluid solver using DMDA+KSP.\n\n";
#include <petsc.h>
#include <mpi.h>

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include <fstream>
#include <string>

#include "wyh_output_result.h"
#include "wyh_solve_fluid_structure_interaction_system3D200630.h"
#include "wyh_quadratic_rectangular_element_interpolation3D200527.h"

using namespace std;

//转换=============================================================================================
void change_variable_PETSC_to_CPP(MPI_Comm PETSC_COMM_WORLD, Vec PETSC_vector_value, double *CPP_vector_value, int CPP_num, int Petsc_size, int Petsc_rank);
void change_variable_CPP_to_PETSC(double *CPP_vector_value, Vec PETSC_vector_value);

//初始网格=========================================================================================
void get_initial_value(double * CPP_Lagrange_node_x, double * CPP_Lagrange_node_y, double * CPP_Lagrange_node_z, int Lagrange_node_num);

//读入数据=========================================================================================
void input_data(
	double * CPP_velocity_u, double * CPP_velocity_v, double * CPP_velocity_w, double * CPP_pressure,
	double * CPP_FF00, double * CPP_FF01, double * CPP_FF02,
	double * CPP_FF10, double * CPP_FF11, double * CPP_FF12,
	double * CPP_FF20, double * CPP_FF21, double * CPP_FF22,
	double * CPP_FF_inv00, double * CPP_FF_inv01, double * CPP_FF_inv02,
	double * CPP_FF_inv10, double * CPP_FF_inv11, double * CPP_FF_inv12,
	double * CPP_FF_inv20, double * CPP_FF_inv21, double * CPP_FF_inv22,
	double * CPP_phi,
	double * CPP_u_NS_convection_ori, double * CPP_v_NS_convection_ori, double * CPP_w_NS_convection_ori,
	double * CPP_Lagrange_node_x, double * CPP_Lagrange_node_y, double * CPP_Lagrange_node_z,
	int all_mesh_num, int Lagrange_node_num, int judge_num);

//主程序===========================================================================================
int main(int argc, char* argv[])
{
	clock_t run_start, run_end;
	run_start = clock();

	//-------------------------------------------
	wyh_solve_fluid_structure_interaction_system3D200630 FS;
	wyh_quadratic_rectangular_element_interpolation3D200527 rect_inter;
	wyh_output_result output;
	
	//PTESc准备----------------------------------

	//初始化-----------------
	PetscMPIInt	Petsc_rank;
	PetscMPIInt Petsc_size;

	PetscInitialize(&argc, &argv, (char*)0, help);

	MPI_Comm_rank(PETSC_COMM_WORLD, &Petsc_rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &Petsc_size);

	//流固耦合变量准备--------
	FS.set_variable_size(PETSC_COMM_WORLD);
	rect_inter.set_variable_size(PETSC_COMM_WORLD);

	FS.prepare_assembly_matrix_and_vector();

	//输入参数-----------------------------------
	int input_num = atoi(argv[1]);
	int run_step = 10400;
	int judge_num = input_num*run_step;
	int end_num = (input_num + 1)*run_step;

	double t = FS.dt*judge_num;
	double t_end = FS.dt*end_num;

	if (Petsc_rank == FS.root){
		cout << "judge_num is " << judge_num << "; t is " << t << " at begining; " << endl;
	}
	if (Petsc_rank == FS.root){
		cout << "time is " << t << "-" << t_end << "; judge_num is " << judge_num << "; run_step is " << run_step << endl;
		cout << "mesh_num is " << FS.m_num << "*" << FS.n_num << "*" << FS.o_num << "; dt is " << FS.dt << endl;		
	}

	//用于输出的变量-------------------------------

	//CPP--------------------

	//欧拉点
	double * CPP_velocity_u = new double[FS.all_mesh_num];
	double * CPP_velocity_v = new double[FS.all_mesh_num];
	double * CPP_velocity_w = new double[FS.all_mesh_num];
	double * CPP_pressure = new double[FS.all_mesh_num];

	double * CPP_FF00 = new double[FS.all_mesh_num];
	double * CPP_FF01 = new double[FS.all_mesh_num];
	double * CPP_FF02 = new double[FS.all_mesh_num];
	double * CPP_FF10 = new double[FS.all_mesh_num];
	double * CPP_FF11 = new double[FS.all_mesh_num];
	double * CPP_FF12 = new double[FS.all_mesh_num];
	double * CPP_FF20 = new double[FS.all_mesh_num];
	double * CPP_FF21 = new double[FS.all_mesh_num];
	double * CPP_FF22 = new double[FS.all_mesh_num];

	double * CPP_FF_inv00 = new double[FS.all_mesh_num];
	double * CPP_FF_inv01 = new double[FS.all_mesh_num];
	double * CPP_FF_inv02 = new double[FS.all_mesh_num];
	double * CPP_FF_inv10 = new double[FS.all_mesh_num];
	double * CPP_FF_inv11 = new double[FS.all_mesh_num];
	double * CPP_FF_inv12 = new double[FS.all_mesh_num];
	double * CPP_FF_inv20 = new double[FS.all_mesh_num];
	double * CPP_FF_inv21 = new double[FS.all_mesh_num];
	double * CPP_FF_inv22 = new double[FS.all_mesh_num];

	double * CPP_phi = new double[FS.all_mesh_num];

	double * CPP_u_NS_convection_ori = new double[FS.all_mesh_num];
	double * CPP_v_NS_convection_ori = new double[FS.all_mesh_num];
	double * CPP_w_NS_convection_ori = new double[FS.all_mesh_num];

	//拉格朗日点
	double * CPP_Lagrange_node_x = new double[rect_inter.Lagrange_node_num];
	double * CPP_Lagrange_node_y = new double[rect_inter.Lagrange_node_num];
	double * CPP_Lagrange_node_z = new double[rect_inter.Lagrange_node_num];

	//输入网格坐标-------------------------------
	get_initial_value(CPP_Lagrange_node_x, CPP_Lagrange_node_y, CPP_Lagrange_node_z, rect_inter.Lagrange_node_num);
	change_variable_CPP_to_PETSC(CPP_Lagrange_node_x, rect_inter.Lagrange_node_x);
	change_variable_CPP_to_PETSC(CPP_Lagrange_node_y, rect_inter.Lagrange_node_y);
	change_variable_CPP_to_PETSC(CPP_Lagrange_node_z, rect_inter.Lagrange_node_z);

	//是否接受中间数据---------------------------
	if (judge_num != 0){
		input_data(
			CPP_velocity_u, CPP_velocity_v, CPP_velocity_w, CPP_pressure,
			CPP_FF00, CPP_FF01, CPP_FF02,
			CPP_FF10, CPP_FF11, CPP_FF12,
			CPP_FF20, CPP_FF21, CPP_FF22,
			CPP_FF_inv00, CPP_FF_inv01, CPP_FF_inv02,
			CPP_FF_inv10, CPP_FF_inv11, CPP_FF_inv12,
			CPP_FF_inv20, CPP_FF_inv21, CPP_FF_inv22,
			CPP_phi,
			CPP_u_NS_convection_ori, CPP_v_NS_convection_ori, CPP_w_NS_convection_ori,
			CPP_Lagrange_node_x, CPP_Lagrange_node_y, CPP_Lagrange_node_z,
			FS.all_mesh_num, rect_inter.Lagrange_node_num, judge_num);

		change_variable_CPP_to_PETSC(CPP_velocity_u, FS.velocity_u);
		change_variable_CPP_to_PETSC(CPP_velocity_v, FS.velocity_v);
		change_variable_CPP_to_PETSC(CPP_velocity_w, FS.velocity_w);
		change_variable_CPP_to_PETSC(CPP_pressure, FS.pressure);

		change_variable_CPP_to_PETSC(CPP_FF00, FS.FF00);	change_variable_CPP_to_PETSC(CPP_FF01, FS.FF01);	change_variable_CPP_to_PETSC(CPP_FF02, FS.FF02);
		change_variable_CPP_to_PETSC(CPP_FF10, FS.FF10);	change_variable_CPP_to_PETSC(CPP_FF11, FS.FF11);	change_variable_CPP_to_PETSC(CPP_FF12, FS.FF12);
		change_variable_CPP_to_PETSC(CPP_FF20, FS.FF20);	change_variable_CPP_to_PETSC(CPP_FF21, FS.FF21);	change_variable_CPP_to_PETSC(CPP_FF22, FS.FF22);
		
		change_variable_CPP_to_PETSC(CPP_FF_inv00, FS.FF_inv00);	change_variable_CPP_to_PETSC(CPP_FF_inv01, FS.FF_inv01);	change_variable_CPP_to_PETSC(CPP_FF_inv02, FS.FF_inv02);
		change_variable_CPP_to_PETSC(CPP_FF_inv10, FS.FF_inv10);	change_variable_CPP_to_PETSC(CPP_FF_inv11, FS.FF_inv11);	change_variable_CPP_to_PETSC(CPP_FF_inv12, FS.FF_inv12);
		change_variable_CPP_to_PETSC(CPP_FF_inv20, FS.FF_inv20);	change_variable_CPP_to_PETSC(CPP_FF_inv21, FS.FF_inv21);	change_variable_CPP_to_PETSC(CPP_FF_inv22, FS.FF_inv22);

		change_variable_CPP_to_PETSC(CPP_phi, FS.phi);

		change_variable_CPP_to_PETSC(CPP_Lagrange_node_x, rect_inter.Lagrange_node_x);
		change_variable_CPP_to_PETSC(CPP_Lagrange_node_y, rect_inter.Lagrange_node_y);
		change_variable_CPP_to_PETSC(CPP_Lagrange_node_z, rect_inter.Lagrange_node_z);
	}

	//-------------------------------------------
	if (Petsc_rank == FS.root){
		run_end = clock();
		cout << "begin to run and the CPU time is " << (double)(run_end - run_start) / CLOCKS_PER_SEC << "s" << endl;
	}

	while (t <= t_end && t<=10) {
		judge_num = judge_num + 1;
		t = t + FS.dt;
		FS.judge_num = judge_num;

		//求解流固耦合系统------------------------
		FS.solve_fluid_structure_interaction_system3D();
		
		//更细拉格朗日点
		rect_inter.update_Lagrange_node3D(PETSC_COMM_WORLD,
			FS.velocity_u_at_p, FS.velocity_v_at_p, FS.velocity_w_at_p);

		//out put--------------------------------
		if (judge_num == 1 || judge_num % 100 == 0){
			//转换-----------

			//欧拉点
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.velocity_u, CPP_velocity_u, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.velocity_v, CPP_velocity_v, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.velocity_w, CPP_velocity_w, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.pressure, CPP_pressure, FS.all_mesh_num, Petsc_size, Petsc_rank);

			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF00, CPP_FF00, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF01, CPP_FF01, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF02, CPP_FF02, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF10, CPP_FF10, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF11, CPP_FF11, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF12, CPP_FF12, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF20, CPP_FF20, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF21, CPP_FF21, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF22, CPP_FF22, FS.all_mesh_num, Petsc_size, Petsc_rank);

			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv00, CPP_FF_inv00, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv01, CPP_FF_inv01, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv02, CPP_FF_inv02, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv10, CPP_FF_inv10, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv11, CPP_FF_inv11, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv12, CPP_FF_inv12, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv20, CPP_FF_inv20, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv21, CPP_FF_inv21, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.FF_inv22, CPP_FF_inv22, FS.all_mesh_num, Petsc_size, Petsc_rank);

			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.phi, CPP_phi, FS.all_mesh_num, Petsc_size, Petsc_rank);

			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.u_NS_convection_ori, CPP_u_NS_convection_ori, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.v_NS_convection_ori, CPP_v_NS_convection_ori, FS.all_mesh_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, FS.w_NS_convection_ori, CPP_w_NS_convection_ori, FS.all_mesh_num, Petsc_size, Petsc_rank);

			//拉格朗日点
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, rect_inter.Lagrange_node_x, CPP_Lagrange_node_x, rect_inter.Lagrange_node_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, rect_inter.Lagrange_node_y, CPP_Lagrange_node_y, rect_inter.Lagrange_node_num, Petsc_size, Petsc_rank);
			change_variable_PETSC_to_CPP(PETSC_COMM_WORLD, rect_inter.Lagrange_node_z, CPP_Lagrange_node_z, rect_inter.Lagrange_node_num, Petsc_size, Petsc_rank);

			if (Petsc_rank == FS.root){
				run_end = clock();
				cout << "Time is :" << t << ", step is : " << judge_num << "; CPU time is " << (double)(run_end - run_start) / CLOCKS_PER_SEC << "s" << endl;

				int mesh_size[3] = { FS.n_num, FS.m_num, FS.o_num };
				int Lagrange_node_size[2] = { 1, rect_inter.Lagrange_node_num };

				//存储文件名---------------------

				//欧拉点	
				char velocity_u_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/velocity_u" };
				char velocity_v_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/velocity_v" };
				char velocity_w_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/velocity_w" };
				char pressure_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/pressure" };

				char FF00_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF00_" };
				char FF01_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF01_" };
				char FF02_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF02_" };
				char FF10_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF10_" };
				char FF11_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF11_" };
				char FF12_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF12_" };
				char FF20_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF20_" };
				char FF21_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF21_" };
				char FF22_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF22_" };

				char FF_inv00_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv00_" };
				char FF_inv01_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv01_" };
				char FF_inv02_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv02_" };
				char FF_inv10_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv10_" };
				char FF_inv11_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv11_" };
				char FF_inv12_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv12_" };
				char FF_inv20_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv20_" };
				char FF_inv21_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv21_" };
				char FF_inv22_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv22_" };

				char phi_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/phi" };

				char u_NS_convection_ori_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/u_NS_convection_ori" };
				char v_NS_convection_ori_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/v_NS_convection_ori" };
				char w_NS_convection_ori_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/w_NS_convection_ori" };

				//拉格朗日点
				char Lagrange_node_x_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/Lagrange_node_x" };
				char Lagrange_node_y_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/Lagrange_node_y" };
				char Lagrange_node_z_name[] = { "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/Lagrange_node_z" };

				//存储数据-----------------------

				//欧拉点
				output.output_result3D(CPP_velocity_u, mesh_size, velocity_u_name, judge_num);
				output.output_result3D(CPP_velocity_v, mesh_size, velocity_v_name, judge_num);
				output.output_result3D(CPP_velocity_w, mesh_size, velocity_w_name, judge_num);
				output.output_result3D(CPP_pressure, mesh_size, pressure_name, judge_num);

				output.output_result3D(CPP_FF00, mesh_size, FF00_name, judge_num);
				output.output_result3D(CPP_FF01, mesh_size, FF01_name, judge_num);
				output.output_result3D(CPP_FF02, mesh_size, FF02_name, judge_num);
				output.output_result3D(CPP_FF10, mesh_size, FF10_name, judge_num);
				output.output_result3D(CPP_FF11, mesh_size, FF11_name, judge_num);
				output.output_result3D(CPP_FF12, mesh_size, FF12_name, judge_num);
				output.output_result3D(CPP_FF20, mesh_size, FF20_name, judge_num);
				output.output_result3D(CPP_FF21, mesh_size, FF21_name, judge_num);
				output.output_result3D(CPP_FF22, mesh_size, FF22_name, judge_num);

				output.output_result3D(CPP_FF_inv00, mesh_size, FF_inv00_name, judge_num);
				output.output_result3D(CPP_FF_inv01, mesh_size, FF_inv01_name, judge_num);
				output.output_result3D(CPP_FF_inv02, mesh_size, FF_inv02_name, judge_num);
				output.output_result3D(CPP_FF_inv10, mesh_size, FF_inv10_name, judge_num);
				output.output_result3D(CPP_FF_inv11, mesh_size, FF_inv11_name, judge_num);
				output.output_result3D(CPP_FF_inv12, mesh_size, FF_inv12_name, judge_num);
				output.output_result3D(CPP_FF_inv20, mesh_size, FF_inv20_name, judge_num);
				output.output_result3D(CPP_FF_inv21, mesh_size, FF_inv21_name, judge_num);
				output.output_result3D(CPP_FF_inv22, mesh_size, FF_inv22_name, judge_num);

				output.output_result3D(CPP_phi, mesh_size, phi_name, judge_num);

				output.output_result3D(CPP_u_NS_convection_ori, mesh_size, u_NS_convection_ori_name, judge_num);
				output.output_result3D(CPP_v_NS_convection_ori, mesh_size, v_NS_convection_ori_name, judge_num);
				output.output_result3D(CPP_w_NS_convection_ori, mesh_size, w_NS_convection_ori_name, judge_num);

				//拉格朗日点
				output.output_result2D(CPP_Lagrange_node_x, Lagrange_node_size, Lagrange_node_x_name, judge_num);
				output.output_result2D(CPP_Lagrange_node_y, Lagrange_node_size, Lagrange_node_y_name, judge_num);
				output.output_result2D(CPP_Lagrange_node_z, Lagrange_node_size, Lagrange_node_z_name, judge_num);
			}
			MPI_Barrier(PETSC_COMM_WORLD);
		}
	}

	MPI_Barrier(PETSC_COMM_WORLD);
	FS.destroy_variable();

	PetscFinalize(); 
	
	run_end = clock();
	if (Petsc_rank == FS.root){
		cout << (double)(run_end - run_start) / CLOCKS_PER_SEC << "s" << endl;
	}

	//释放内存
	delete[] CPP_velocity_u;
	delete[] CPP_velocity_v;
	delete[] CPP_pressure;

	return 0;
}


//*************************************************************************************************

//变量转换=========================================================================================

//*************************************************************************************************

//change PETSC to CPP
void change_variable_PETSC_to_CPP(MPI_Comm PETSC_COMM_WORLD, Vec PETSC_vector_value, double *CPP_vector_value, int CPP_num, int Petsc_size, int Petsc_rank)
{
	PetscInt istart;
	PetscInt iend;
	int root = 0;

	//
	VecGetOwnershipRange(PETSC_vector_value, &istart, &iend);

	//确定开始、总量、位移量
	int matrix_send_count = iend - istart;
	double *CPP_vector_value_temp = new double[matrix_send_count];
	int *begin_table = new int[Petsc_size];
	int *rank_num_table = new int[Petsc_size];
	int *displace_table = new int[Petsc_size];

	int row_num;
	int row_num_temp;

	//rank num
	for (int i = 0; i < Petsc_size; i++){
		row_num = CPP_num / Petsc_size;
		if (CPP_num % Petsc_size == 0){
			rank_num_table[i] = row_num;
		}
		else{
			row_num_temp = CPP_num - row_num*Petsc_size;
			if (i < row_num_temp){
				rank_num_table[i] = row_num + 1;
			}
			else{
				rank_num_table[i] = row_num;
			}
		}
	}

	//displace
	displace_table[0] = 0;
	for (int i = 1; i < Petsc_size; i++){
		displace_table[i] = displace_table[i - 1] + rank_num_table[i - 1];
	}

	//begin
	begin_table[0] = 0;
	for (int i = 1; i < Petsc_size; i++){
		begin_table[i] = begin_table[i - 1] + rank_num_table[i - 1];
	}

	if (begin_table[Petsc_rank] != istart){
		cout << "error in PETSC to CPP" << endl;
	}

	//组装CPP
	double *PETSC_vector_value_temp;
	VecGetArray(PETSC_vector_value, &PETSC_vector_value_temp);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		CPP_vector_value_temp[Ii - istart] = PETSC_vector_value_temp[Ii - istart];
	}
	VecRestoreArray(PETSC_vector_value, &PETSC_vector_value_temp);
	MPI_Gatherv(CPP_vector_value_temp, matrix_send_count, MPI_DOUBLE, CPP_vector_value, rank_num_table, displace_table, MPI_DOUBLE, root, PETSC_COMM_WORLD);
	MPI_Bcast(CPP_vector_value, CPP_num, MPI_DOUBLE, root, PETSC_COMM_WORLD);
	MPI_Barrier(PETSC_COMM_WORLD);

	//
	delete[] CPP_vector_value_temp;
	delete[] begin_table;
	delete[] rank_num_table;
	delete[] displace_table;
	delete[] PETSC_vector_value_temp;

	CPP_vector_value_temp = NULL;
	begin_table = NULL;
	rank_num_table = NULL;
	displace_table = NULL;
	PETSC_vector_value_temp = NULL;
}

//change CPP to PETSC
void change_variable_CPP_to_PETSC(double *CPP_vector_value, Vec PETSC_vector_value)
{
	PetscInt istart;
	PetscInt iend;

	double *PETSC_vector_value_temp;

	VecGetArray(PETSC_vector_value, &PETSC_vector_value_temp);
	VecGetOwnershipRange(PETSC_vector_value, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PETSC_vector_value_temp[Ii - istart] = CPP_vector_value[Ii];
	}
	VecRestoreArray(PETSC_vector_value, &PETSC_vector_value_temp);

	delete[] PETSC_vector_value_temp;
	PETSC_vector_value_temp = NULL;
}

//*************************************************************************************************

//读入数据=========================================================================================

//*************************************************************************************************

//初值-------------------------------------------
void get_initial_value(double * CPP_Lagrange_node_x, double * CPP_Lagrange_node_y, double * CPP_Lagrange_node_z, int Lagrange_node_num)
{
	//文件名
	char input_Lagrange_node_name[500];
	sprintf(input_Lagrange_node_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/cricle_domain_nodes.txt");
	ifstream fin_Lagrange_node(input_Lagrange_node_name);

	//读入剖分单元
	for (int i = 0; i < Lagrange_node_num; i++) {
		for (int j = 0; j < 3; j++){
			if (j == 0){
				fin_Lagrange_node >> CPP_Lagrange_node_x[i];
				CPP_Lagrange_node_x[i] = CPP_Lagrange_node_x[i] + 0.6;
			}
			if (j == 1){
				fin_Lagrange_node >> CPP_Lagrange_node_y[i];
				CPP_Lagrange_node_y[i] = CPP_Lagrange_node_y[i] + 0.5;
			}
			if (j == 2){
				fin_Lagrange_node >> CPP_Lagrange_node_z[i];
				CPP_Lagrange_node_z[i] = CPP_Lagrange_node_z[i] + 0.5;
			}
		}
	}
}

void input_data(
	double * CPP_velocity_u, double * CPP_velocity_v, double * CPP_velocity_w, double * CPP_pressure,
	double * CPP_FF00, double * CPP_FF01, double * CPP_FF02,
	double * CPP_FF10, double * CPP_FF11, double * CPP_FF12,
	double * CPP_FF20, double * CPP_FF21, double * CPP_FF22,
	double * CPP_FF_inv00, double * CPP_FF_inv01, double * CPP_FF_inv02,
	double * CPP_FF_inv10, double * CPP_FF_inv11, double * CPP_FF_inv12,
	double * CPP_FF_inv20, double * CPP_FF_inv21, double * CPP_FF_inv22,
	double * CPP_phi,
	double * CPP_u_NS_convection_ori, double * CPP_v_NS_convection_ori, double * CPP_w_NS_convection_ori,
	double * CPP_Lagrange_node_x, double * CPP_Lagrange_node_y, double * CPP_Lagrange_node_z,
	int all_mesh_num, int Lagrange_node_num, int judge_num)
{
	//文件名
	char input_Lagrange_node_x_name[500];
	char input_Lagrange_node_y_name[500];
	char input_Lagrange_node_z_name[500];

	char input_velocity_u_name[500];	char input_velocity_v_name[500];	char input_velocity_w_name[500];
	char input_pressure_name[500];
	char input_phi_name[500];

	char input_u_NS_convection_ori_name[500];
	char input_v_NS_convection_ori_name[500];
	char input_w_NS_convection_ori_name[500];

	char input_FF00_name[500];		char input_FF01_name[500];		char input_FF02_name[500];
	char input_FF10_name[500];		char input_FF11_name[500];		char input_FF12_name[500];
	char input_FF20_name[500];		char input_FF21_name[500];		char input_FF22_name[500];

	char input_FF_inv00_name[500];		char input_FF_inv01_name[500];		char input_FF_inv02_name[500];
	char input_FF_inv10_name[500];		char input_FF_inv11_name[500];		char input_FF_inv12_name[500];
	char input_FF_inv20_name[500];		char input_FF_inv21_name[500];		char input_FF_inv22_name[500];

	sprintf(input_Lagrange_node_x_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/Lagrange_node_x%d%s", judge_num, ".txt");
	sprintf(input_Lagrange_node_y_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/Lagrange_node_y%d%s", judge_num, ".txt");
	sprintf(input_Lagrange_node_z_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/Lagrange_node_z%d%s", judge_num, ".txt");

	sprintf(input_velocity_u_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/velocity_u%d%s", judge_num, ".txt");
	sprintf(input_velocity_v_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/velocity_v%d%s", judge_num, ".txt");
	sprintf(input_velocity_w_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/velocity_w%d%s", judge_num, ".txt");

	sprintf(input_pressure_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/pressure%d%s", judge_num, ".txt");
	sprintf(input_phi_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/phi%d%s", judge_num, ".txt");

	sprintf(input_u_NS_convection_ori_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/u_NS_convection_ori%d%s", judge_num, ".txt");
	sprintf(input_v_NS_convection_ori_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/v_NS_convection_ori%d%s", judge_num, ".txt");
	sprintf(input_w_NS_convection_ori_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/w_NS_convection_ori%d%s", judge_num, ".txt");

	sprintf(input_FF00_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF00_%d%s", judge_num, ".txt");
	sprintf(input_FF10_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF10_%d%s", judge_num, ".txt");
	sprintf(input_FF20_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF20_%d%s", judge_num, ".txt");

	sprintf(input_FF01_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF01_%d%s", judge_num, ".txt");
	sprintf(input_FF11_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF11_%d%s", judge_num, ".txt");
	sprintf(input_FF21_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF21_%d%s", judge_num, ".txt");

	sprintf(input_FF02_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF02_%d%s", judge_num, ".txt");
	sprintf(input_FF12_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF12_%d%s", judge_num, ".txt");
	sprintf(input_FF22_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF22_%d%s", judge_num, ".txt");

	sprintf(input_FF_inv00_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv00_%d%s", judge_num, ".txt");
	sprintf(input_FF_inv10_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv10_%d%s", judge_num, ".txt");
	sprintf(input_FF_inv20_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv20_%d%s", judge_num, ".txt");

	sprintf(input_FF_inv01_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv01_%d%s", judge_num, ".txt");
	sprintf(input_FF_inv11_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv11_%d%s", judge_num, ".txt");
	sprintf(input_FF_inv21_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv21_%d%s", judge_num, ".txt");

	sprintf(input_FF_inv02_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv02_%d%s", judge_num, ".txt");
	sprintf(input_FF_inv12_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv12_%d%s", judge_num, ".txt");
	sprintf(input_FF_inv22_name, "/home/large/ibamr/wyh/wyh_PETSc/PETSC_fulid_solid_thirteen/fluid_solid/3D/lid_driven_cavity_flow/data/FF_inv22_%d%s", judge_num, ".txt");

	ifstream fin_Lagrange_node_x(input_Lagrange_node_x_name);
	ifstream fin_Lagrange_node_y(input_Lagrange_node_y_name);
	ifstream fin_Lagrange_node_z(input_Lagrange_node_z_name);

	ifstream fin_velocity_u(input_velocity_u_name);
	ifstream fin_velocity_v(input_velocity_v_name);
	ifstream fin_velocity_w(input_velocity_w_name);

	ifstream fin_pressure(input_pressure_name);
	ifstream fin_phi(input_phi_name);

	ifstream fin_u_NS_convection_ori(input_u_NS_convection_ori_name);
	ifstream fin_v_NS_convection_ori(input_v_NS_convection_ori_name);
	ifstream fin_w_NS_convection_ori(input_w_NS_convection_ori_name);

	ifstream fin_FF00(input_FF00_name);	ifstream fin_FF01(input_FF01_name);	ifstream fin_FF02(input_FF02_name);
	ifstream fin_FF10(input_FF10_name);	ifstream fin_FF11(input_FF11_name);	ifstream fin_FF12(input_FF12_name);
	ifstream fin_FF20(input_FF20_name);	ifstream fin_FF21(input_FF21_name);	ifstream fin_FF22(input_FF22_name);

	ifstream fin_FF_inv00(input_FF_inv00_name);	ifstream fin_FF_inv01(input_FF_inv01_name);	ifstream fin_FF_inv02(input_FF_inv02_name);
	ifstream fin_FF_inv10(input_FF_inv10_name);	ifstream fin_FF_inv11(input_FF_inv11_name);	ifstream fin_FF_inv12(input_FF_inv12_name);
	ifstream fin_FF_inv20(input_FF_inv20_name);	ifstream fin_FF_inv21(input_FF_inv21_name);	ifstream fin_FF_inv22(input_FF_inv22_name);

	//读入剖分单元
	for (int i = 0; i < all_mesh_num; i++) {
		fin_velocity_u >> CPP_velocity_u[i];
		fin_velocity_v >> CPP_velocity_v[i];
		fin_velocity_w >> CPP_velocity_w[i];

		fin_pressure >> CPP_pressure[i];

		fin_FF00 >> CPP_FF00[i];		fin_FF01 >> CPP_FF01[i];		fin_FF02 >> CPP_FF02[i];
		fin_FF10 >> CPP_FF10[i];		fin_FF11 >> CPP_FF11[i];		fin_FF12 >> CPP_FF12[i];
		fin_FF20 >> CPP_FF20[i];		fin_FF21 >> CPP_FF21[i];		fin_FF22 >> CPP_FF22[i];

		fin_FF_inv00 >> CPP_FF_inv00[i];		fin_FF_inv01 >> CPP_FF_inv01[i];		fin_FF_inv02 >> CPP_FF_inv02[i];
		fin_FF_inv10 >> CPP_FF_inv10[i];		fin_FF_inv11 >> CPP_FF_inv11[i];		fin_FF_inv12 >> CPP_FF_inv12[i];
		fin_FF_inv20 >> CPP_FF_inv20[i];		fin_FF_inv21 >> CPP_FF_inv21[i];		fin_FF_inv22 >> CPP_FF_inv22[i];

		fin_phi >> CPP_phi[i];

		fin_u_NS_convection_ori >> CPP_u_NS_convection_ori[i];
		fin_v_NS_convection_ori >> CPP_v_NS_convection_ori[i];
		fin_w_NS_convection_ori >> CPP_w_NS_convection_ori[i];
	}

	for (int i = 0; i < Lagrange_node_num; i++){
		fin_Lagrange_node_x >> CPP_Lagrange_node_x[i];
		fin_Lagrange_node_y >> CPP_Lagrange_node_y[i];
		fin_Lagrange_node_z >> CPP_Lagrange_node_z[i];
	}
}
