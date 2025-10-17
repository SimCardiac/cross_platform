#include "wyh_solve_fluid_structure_interaction_system3D200630.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include <petsc.h>
#include <mpi.h>

using namespace std;

//***************************************************************************************************************************************************

//设定变量大小并赋值==================================================================================================================================

//***************************************************************************************************************************************************

//变量大小==========================================================================================
void wyh_solve_fluid_structure_interaction_system3D200630::set_variable_size(MPI_Comm PETSC_COMM_WORLD)
{
	//与实际问题相关变量-------------------------

	//方程组系数
	MatCreate(PETSC_COMM_WORLD, &NS_coeff_u3D);
	MatSetSizes(NS_coeff_u3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num_take_out_x_L, all_mesh_num_take_out_x_L);
	MatSetFromOptions(NS_coeff_u3D);

	MatCreate(PETSC_COMM_WORLD, &NS_coeff_v3D);
	MatSetSizes(NS_coeff_v3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num_take_out_y_L, all_mesh_num_take_out_y_L);
	MatSetFromOptions(NS_coeff_v3D);

	MatCreate(PETSC_COMM_WORLD, &NS_coeff_w3D);
	MatSetSizes(NS_coeff_w3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num_take_out_z_L, all_mesh_num_take_out_z_L);
	MatSetFromOptions(NS_coeff_w3D);

	MatCreate(PETSC_COMM_WORLD, &NS_coeff_p3D);
	MatSetSizes(NS_coeff_p3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num);
	MatSetFromOptions(NS_coeff_p3D);

	//修正方程组右端项
	VecCreate(PETSC_COMM_WORLD, &modified_u_NSeq_vector);
	VecSetSizes(modified_u_NSeq_vector, PETSC_DECIDE, all_mesh_num);
	VecSetFromOptions(modified_u_NSeq_vector);

	//去边界的过度速度
	VecCreate(PETSC_COMM_WORLD, &velocity_u_star_take_out_x_L);
	VecSetSizes(velocity_u_star_take_out_x_L, PETSC_DECIDE, all_mesh_num_take_out_x_L);
	VecSetFromOptions(velocity_u_star_take_out_x_L);

	VecCreate(PETSC_COMM_WORLD, &velocity_v_star_take_out_y_L);
	VecSetSizes(velocity_v_star_take_out_y_L, PETSC_DECIDE, all_mesh_num_take_out_y_L);
	VecSetFromOptions(velocity_v_star_take_out_y_L);

	VecCreate(PETSC_COMM_WORLD, &velocity_w_star_take_out_z_L);
	VecSetSizes(velocity_w_star_take_out_z_L, PETSC_DECIDE, all_mesh_num_take_out_z_L);
	VecSetFromOptions(velocity_w_star_take_out_z_L);

	//去边界的NS方程右端项
	VecDuplicate(velocity_u_star_take_out_x_L, &u_NS_right_term_take_out_x_L);
	VecDuplicate(velocity_v_star_take_out_y_L, &v_NS_right_term_take_out_y_L);
	VecDuplicate(velocity_w_star_take_out_z_L, &w_NS_right_term_take_out_z_L);

	//去掉边界的矩阵
	MatCreate(PETSC_COMM_WORLD, &take_out_x_L_boundary_matrix3D);
	MatSetSizes(take_out_x_L_boundary_matrix3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num_take_out_x_L, all_mesh_num);
	MatSetFromOptions(take_out_x_L_boundary_matrix3D);

	MatCreate(PETSC_COMM_WORLD, &take_out_y_L_boundary_matrix3D);
	MatSetSizes(take_out_y_L_boundary_matrix3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num_take_out_y_L, all_mesh_num);
	MatSetFromOptions(take_out_y_L_boundary_matrix3D);

	MatCreate(PETSC_COMM_WORLD, &take_out_z_L_boundary_matrix3D);
	MatSetSizes(take_out_z_L_boundary_matrix3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num_take_out_z_L, all_mesh_num);
	MatSetFromOptions(take_out_z_L_boundary_matrix3D);

	//添加边界的矩阵
	MatCreate(PETSC_COMM_WORLD, &x_solution_out_x_L_extend_to_ori_size_matrix3D);
	MatSetSizes(x_solution_out_x_L_extend_to_ori_size_matrix3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num_take_out_x_L);
	MatSetFromOptions(x_solution_out_x_L_extend_to_ori_size_matrix3D);

	MatCreate(PETSC_COMM_WORLD, &y_solution_out_y_L_extend_to_ori_size_matrix3D);
	MatSetSizes(y_solution_out_y_L_extend_to_ori_size_matrix3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num_take_out_y_L);
	MatSetFromOptions(y_solution_out_y_L_extend_to_ori_size_matrix3D);

	MatCreate(PETSC_COMM_WORLD, &z_solution_out_z_L_extend_to_ori_size_matrix3D);
	MatSetSizes(z_solution_out_z_L_extend_to_ori_size_matrix3D, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num_take_out_z_L);
	MatSetFromOptions(z_solution_out_z_L_extend_to_ori_size_matrix3D);

	//全局变量-----------------------------------

	//速度...................
	VecCreate(PETSC_COMM_WORLD, &velocity_u);
	VecSetSizes(velocity_u, PETSC_DECIDE, all_mesh_num);
	VecSetFromOptions(velocity_u);
	VecDuplicate(velocity_u, &velocity_v);
	VecDuplicate(velocity_u, &velocity_w);

	VecDuplicate(velocity_u, &velocity_u_at_p);
	VecDuplicate(velocity_u, &velocity_v_at_p);
	VecDuplicate(velocity_u, &velocity_w_at_p);

	VecDuplicate(velocity_u, &u_NS_convection_ori);
	VecDuplicate(velocity_u, &v_NS_convection_ori);
	VecDuplicate(velocity_w, &w_NS_convection_ori);

	//压力...................
	VecDuplicate(velocity_u, &pressure);

	//phi....................
	VecDuplicate(velocity_u, &phi);

	//FF.....................
	VecDuplicate(velocity_u, &FF00);	VecDuplicate(velocity_u, &FF01);	VecDuplicate(velocity_u, &FF02);
	VecDuplicate(velocity_u, &FF10);	VecDuplicate(velocity_u, &FF11);	VecDuplicate(velocity_u, &FF12);
	VecDuplicate(velocity_u, &FF20);	VecDuplicate(velocity_u, &FF21);	VecDuplicate(velocity_u, &FF22);

	//形变梯度的逆...........
	VecDuplicate(FF00, &FF_inv00);	VecDuplicate(FF00, &FF_inv01);	VecDuplicate(FF00, &FF_inv02);
	VecDuplicate(FF00, &FF_inv10);	VecDuplicate(FF00, &FF_inv11);	VecDuplicate(FF00, &FF_inv12);
	VecDuplicate(FF00, &FF_inv20);	VecDuplicate(FF00, &FF_inv21);	VecDuplicate(FF00, &FF_inv22);

	//应力张量...............
	VecDuplicate(FF00, &tensor00);	VecDuplicate(FF00, &tensor01);	VecDuplicate(FF00, &tensor02);
	VecDuplicate(FF00, &tensor10);	VecDuplicate(FF00, &tensor11);	VecDuplicate(FF00, &tensor12);
	VecDuplicate(FF00, &tensor20);	VecDuplicate(FF00, &tensor21);	VecDuplicate(FF00, &tensor22);

	//应力...................
	VecDuplicate(FF00, &F_x_at_p);
	VecDuplicate(FF00, &F_y_at_p);
	VecDuplicate(FF00, &F_z_at_p);

	VecDuplicate(FF00, &tensor_body_x_at_p);
	VecDuplicate(FF00, &tensor_body_y_at_p);
	VecDuplicate(FF00, &tensor_body_z_at_p);

	VecDuplicate(FF00, &tensor_normal_x_at_p);
	VecDuplicate(FF00, &tensor_normal_y_at_p);
	VecDuplicate(FF00, &tensor_normal_z_at_p);

	//方程相关-----------------------------------

	//修正项
	VecDuplicate(velocity_u, &velocity_u_star);
	VecDuplicate(velocity_u, &velocity_v_star);
	VecDuplicate(velocity_u, &velocity_w_star);
	VecDuplicate(velocity_u, &modified_term_at_p);

	//NS方程右端项
	VecDuplicate(velocity_u, &u_NS_right_term);
	VecDuplicate(velocity_u, &v_NS_right_term);
	VecDuplicate(velocity_u, &w_NS_right_term);
	VecDuplicate(velocity_u, &modified_equation_right_term);

	//应力
	VecDuplicate(u_NS_right_term, &F_x_at_u);
	VecDuplicate(v_NS_right_term, &F_y_at_v);
	VecDuplicate(w_NS_right_term, &F_z_at_w);

	//压力
	VecDuplicate(u_NS_right_term, &d_pressure_x_at_u);
	VecDuplicate(v_NS_right_term, &d_pressure_y_at_v);
	VecDuplicate(w_NS_right_term, &d_pressure_z_at_w);

	//对流项
	VecDuplicate(u_NS_right_term, &u_NS_convection);
	VecDuplicate(v_NS_right_term, &v_NS_convection);
	VecDuplicate(w_NS_right_term, &w_NS_convection);

	//扩散项
	VecDuplicate(u_NS_right_term, &u_NS_diffusion);
	VecDuplicate(v_NS_right_term, &v_NS_diffusion);
	VecDuplicate(w_NS_right_term, &w_NS_diffusion);

	//法线向量-----------------------------------
	VecDuplicate(FF00, &phi_normal_x);
	VecDuplicate(FF00, &phi_normal_y);
	VecDuplicate(FF00, &phi_normal_z);

	VecDuplicate(FF00, &normal_norm);
	VecDuplicate(FF00, &square_phi_normal_x);
	VecDuplicate(FF00, &square_phi_normal_y);

	//区域划分-----------------------------------

	//内部区域
	VecDuplicate(FF00, &extract_interior_domain_vector);

	//边界区域
	VecDuplicate(FF00, &extract_boundary_domain_vector);
	
	//修正phi----------------
	VecDuplicate(FF00, &ones_vector);
	VecDuplicate(FF00, &one_redu_phi);

	VecDuplicate(FF00, &phi_boundary_mark_vector);
	VecDuplicate(FF00, &sign_boundary_mark);

	VecDuplicate(FF00, &modified_right_term);
	VecDuplicate(FF00, &modified_part1);
	VecDuplicate(FF00, &modified_part2);
	VecDuplicate(FF00, &modified_part_x);
	VecDuplicate(FF00, &modified_part_y);
	VecDuplicate(FF00, &modified_part_z);

	VecDuplicate(FF00, &d_phi_x);
	VecDuplicate(FF00, &d_phi_y);
	VecDuplicate(FF00, &d_phi_z);
	VecDuplicate(FF00, &d_phi_norm);
	
	//导数---------------------------------------

	//WENO相关
	VecDuplicate(phi, &d_uphi_x);
	VecDuplicate(phi, &d_vphi_y);

	VecDuplicate(u_NS_convection, &d_uu_x);
	VecDuplicate(u_NS_convection, &d_vu_y);
	VecDuplicate(u_NS_convection, &d_uv_x);
	VecDuplicate(u_NS_convection, &d_vv_y);
	VecDuplicate(u_NS_convection, &d_uw_x);
	VecDuplicate(u_NS_convection, &d_vw_y);

	//中间速度1阶导数
	VecDuplicate(velocity_u_at_p, &d_u_x_at_p);
	VecDuplicate(velocity_u_at_p, &d_u_y_at_p);
	VecDuplicate(velocity_u_at_p, &d_u_z_at_p);

	VecDuplicate(velocity_u_at_p, &d_v_x_at_p);
	VecDuplicate(velocity_u_at_p, &d_v_y_at_p);
	VecDuplicate(velocity_u_at_p, &d_v_z_at_p);

	VecDuplicate(velocity_u_at_p, &d_w_x_at_p);
	VecDuplicate(velocity_u_at_p, &d_w_y_at_p);
	VecDuplicate(velocity_u_at_p, &d_w_z_at_p);

	//过渡速度1阶导数
	VecDuplicate(velocity_u_star, &d_u_star_x);
	VecDuplicate(velocity_v_star, &d_v_star_y);

	//速度二阶导数
	VecDuplicate(u_NS_diffusion, &dd_u_xx);
	VecDuplicate(u_NS_diffusion, &dd_u_yy);

	VecDuplicate(v_NS_diffusion, &dd_v_xx);
	VecDuplicate(v_NS_diffusion, &dd_v_yy);

	VecDuplicate(w_NS_diffusion, &dd_w_xx);
	VecDuplicate(w_NS_diffusion, &dd_w_yy);

	//应力...................
	VecDuplicate(FF00, &d_tensor00_x);	VecDuplicate(FF00, &d_tensor01_x);	VecDuplicate(FF00, &d_tensor02_x);
	VecDuplicate(FF00, &d_tensor10_y);	VecDuplicate(FF00, &d_tensor11_y);	VecDuplicate(FF00, &d_tensor12_y);
	VecDuplicate(FF00, &d_tensor20_z);	VecDuplicate(FF00, &d_tensor21_z);	VecDuplicate(FF00, &d_tensor22_z);

	VecDuplicate(modified_term_at_p, &d_modified_term_x);
	VecDuplicate(modified_term_at_p, &d_modified_term_y);
	VecDuplicate(modified_term_at_p, &d_modified_term_z);

	//修正压力...............
	VecDuplicate(modified_term_at_p, &dd_modified_term_xx);
	VecDuplicate(modified_term_at_p, &dd_modified_term_yy);
	VecDuplicate(modified_term_at_p, &lap_modified_term);

	//过渡项-------------------------------------
	VecDuplicate(FF00, &right_term_part0);
	VecDuplicate(FF00, &right_term_part1);
	VecDuplicate(FF00, &right_term_part1_0);
	VecDuplicate(FF00, &right_term_part1_1);

	VecDuplicate(FF00, &mult_part01);
	VecDuplicate(FF00, &mult_part02);
	VecDuplicate(FF00, &tensor00_transition);	VecDuplicate(FF00, &tensor01_transition);	VecDuplicate(FF00, &tensor02_transition);
	VecDuplicate(FF00, &tensor10_transition);	VecDuplicate(FF00, &tensor11_transition);	VecDuplicate(FF00, &tensor12_transition);
	VecDuplicate(FF00, &tensor20_transition);	VecDuplicate(FF00, &tensor21_transition);	VecDuplicate(FF00, &tensor22_transition);

	//拓展变量-----------------------------------

	//速度...................

	//
	VecCreate(PETSC_COMM_WORLD, &velocity_u_star_at_u010000);
	VecSetSizes(velocity_u_star_at_u010000, PETSC_DECIDE, all_mesh_num100);
	VecSetFromOptions(velocity_u_star_at_u010000);

	VecCreate(PETSC_COMM_WORLD, &velocity_v_star_at_v000100);
	VecSetSizes(velocity_v_star_at_v000100, PETSC_DECIDE, all_mesh_num010);
	VecSetFromOptions(velocity_v_star_at_v000100);

	VecCreate(PETSC_COMM_WORLD, &velocity_w_star_at_w000001);
	VecSetSizes(velocity_w_star_at_w000001, PETSC_DECIDE, all_mesh_num001);
	VecSetFromOptions(velocity_w_star_at_w000001);

	//
	VecCreate(PETSC_COMM_WORLD, &velocity_u_at_u010000);
	VecSetSizes(velocity_u_at_u010000, PETSC_DECIDE, all_mesh_num100);
	VecSetFromOptions(velocity_u_at_u010000);

	VecCreate(PETSC_COMM_WORLD, &velocity_v_at_v000100);
	VecSetSizes(velocity_v_at_v000100, PETSC_DECIDE, all_mesh_num010);
	VecSetFromOptions(velocity_v_at_v000100);

	VecCreate(PETSC_COMM_WORLD, &velocity_w_at_w000001);
	VecSetSizes(velocity_w_at_w000001, PETSC_DECIDE, all_mesh_num001);
	VecSetFromOptions(velocity_w_at_w000001);

	//
	VecCreate(PETSC_COMM_WORLD, &velocity_u_at_u111111);
	VecSetSizes(velocity_u_at_u111111, PETSC_DECIDE, all_mesh_num222);
	VecSetFromOptions(velocity_u_at_u111111);
	VecDuplicate(velocity_u_at_u111111, &velocity_v_at_v111111);
	VecDuplicate(velocity_u_at_u111111, &velocity_w_at_w111111);

	//
	VecCreate(PETSC_COMM_WORLD, &velocity_u_at_u333333);
	VecSetSizes(velocity_u_at_u333333, PETSC_DECIDE, all_mesh_num666);
	VecSetFromOptions(velocity_u_at_u333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_u_at_v333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_u_at_w333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_v_at_u333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_v_at_v333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_v_at_w333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_w_at_u333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_w_at_v333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_w_at_w333333);

	//
	VecCreate(PETSC_COMM_WORLD, &velocity_u_at_u344333);
	VecSetSizes(velocity_u_at_u344333, PETSC_DECIDE, all_mesh_num776);
	VecSetFromOptions(velocity_u_at_u344333);
	VecDuplicate(velocity_u_at_u344333, &velocity_v_at_v433433);

	VecCreate(PETSC_COMM_WORLD, &velocity_u_at_u343343);
	VecSetSizes(velocity_u_at_u343343, PETSC_DECIDE, all_mesh_num767);
	VecSetFromOptions(velocity_u_at_u343343);
	VecDuplicate(velocity_u_at_u343343, &velocity_w_at_w433334);

	VecCreate(PETSC_COMM_WORLD, &velocity_v_at_v333443);
	VecSetSizes(velocity_v_at_v333443, PETSC_DECIDE, all_mesh_num677);
	VecSetFromOptions(velocity_v_at_v333443);
	VecDuplicate(velocity_v_at_v333443, &velocity_w_at_w334334);

	//
	VecCreate(PETSC_COMM_WORLD, &velocity_u_at_p111111);
	VecSetSizes(velocity_u_at_p111111, PETSC_DECIDE, all_mesh_num222);
	VecSetFromOptions(velocity_u_at_p111111);

	VecCreate(PETSC_COMM_WORLD, &velocity_v_at_p111111);
	VecSetSizes(velocity_v_at_p111111, PETSC_DECIDE, all_mesh_num222);
	VecSetFromOptions(velocity_v_at_p111111);

	VecCreate(PETSC_COMM_WORLD, &velocity_w_at_p111111);
	VecSetSizes(velocity_w_at_p111111, PETSC_DECIDE, all_mesh_num222);
	VecSetFromOptions(velocity_w_at_p111111);

	VecDuplicate(velocity_u_at_u333333, &velocity_u_at_p333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_v_at_p333333);
	VecDuplicate(velocity_u_at_u333333, &velocity_w_at_p333333);

	//压力...................
	VecCreate(PETSC_COMM_WORLD, &pressure_at_p100000);
	VecSetSizes(pressure_at_p100000, PETSC_DECIDE, all_mesh_num100);
	VecSetFromOptions(pressure_at_p100000);

	VecCreate(PETSC_COMM_WORLD, &pressure_at_p001000);
	VecSetSizes(pressure_at_p001000, PETSC_DECIDE, all_mesh_num010);
	VecSetFromOptions(pressure_at_p001000);

	VecCreate(PETSC_COMM_WORLD, &pressure_at_p000010);
	VecSetSizes(pressure_at_p000010, PETSC_DECIDE, all_mesh_num001);
	VecSetFromOptions(pressure_at_p000010);

	//phi....................
	VecCreate(PETSC_COMM_WORLD, &phi_at_p111111);
	VecSetSizes(phi_at_p111111, PETSC_DECIDE, all_mesh_num222);
	VecSetFromOptions(phi_at_p111111);

	//应力...................
	VecCreate(PETSC_COMM_WORLD, &F_x_at_p100000);
	VecSetSizes(F_x_at_p100000, PETSC_DECIDE, all_mesh_num100);
	VecSetFromOptions(F_x_at_p100000);

	VecCreate(PETSC_COMM_WORLD, &F_y_at_p001000);
	VecSetSizes(F_y_at_p001000, PETSC_DECIDE, all_mesh_num010);
	VecSetFromOptions(F_y_at_p001000);

	VecCreate(PETSC_COMM_WORLD, &F_z_at_p000010);
	VecSetSizes(F_z_at_p000010, PETSC_DECIDE, all_mesh_num001);
	VecSetFromOptions(F_z_at_p000010);

	VecCreate(PETSC_COMM_WORLD, &tensor00_111111);
	VecSetSizes(tensor00_111111, PETSC_DECIDE, all_mesh_num222);
	VecSetFromOptions(tensor00_111111);
	VecDuplicate(tensor00_111111, &tensor01_111111);
	VecDuplicate(tensor00_111111, &tensor02_111111);
	VecDuplicate(tensor00_111111, &tensor10_111111);
	VecDuplicate(tensor00_111111, &tensor11_111111);
	VecDuplicate(tensor00_111111, &tensor12_111111);
	VecDuplicate(tensor00_111111, &tensor20_111111);
	VecDuplicate(tensor00_111111, &tensor21_111111);
	VecDuplicate(tensor00_111111, &tensor22_111111);

	//修正项.................
	VecCreate(PETSC_COMM_WORLD, &modified_term_at_p100000);
	VecSetSizes(modified_term_at_p100000, PETSC_DECIDE, all_mesh_num100);
	VecSetFromOptions(modified_term_at_p100000);

	VecCreate(PETSC_COMM_WORLD, &modified_term_at_p001000);
	VecSetSizes(modified_term_at_p001000, PETSC_DECIDE, all_mesh_num010);
	VecSetFromOptions(modified_term_at_p001000);

	VecCreate(PETSC_COMM_WORLD, &modified_term_at_p000010);
	VecSetSizes(modified_term_at_p000010, PETSC_DECIDE, all_mesh_num001);
	VecSetFromOptions(modified_term_at_p000010);

	VecCreate(PETSC_COMM_WORLD, &modified_term_at_p111111);
	VecSetSizes(modified_term_at_p111111, PETSC_DECIDE, all_mesh_num222);
	VecSetFromOptions(modified_term_at_p111111);

	//WENO相关-----------------------------------

	//phi....................
	VecCreate(PETSC_COMM_WORLD, &variable333333);
	VecSetSizes(variable333333, PETSC_DECIDE, all_mesh_num666);
	VecSetFromOptions(variable333333);

	VecDuplicate(variable333333, &f_variable333333);

	VecDuplicate(variable333333, &uphi_f_P333333);		VecDuplicate(variable333333, &uphi_f_N333333);
	VecDuplicate(variable333333, &vphi_f_P333333);		VecDuplicate(variable333333, &vphi_f_N333333);
	VecDuplicate(variable333333, &wphi_f_P333333);		VecDuplicate(variable333333, &wphi_f_N333333);

	//速度
	VecDuplicate(variable333333, &uu_f_P);		VecDuplicate(variable333333, &uu_f_N);
	VecDuplicate(variable333333, &vu_f_P);		VecDuplicate(variable333333, &vu_f_N);
	VecDuplicate(variable333333, &wu_f_P);		VecDuplicate(variable333333, &wu_f_N);
	VecDuplicate(variable333333, &uv_f_P);		VecDuplicate(variable333333, &uv_f_N);
	VecDuplicate(variable333333, &vv_f_P);		VecDuplicate(variable333333, &vv_f_N);
	VecDuplicate(variable333333, &wv_f_P);		VecDuplicate(variable333333, &wv_f_N);
	VecDuplicate(variable333333, &uw_f_P);		VecDuplicate(variable333333, &uw_f_N);
	VecDuplicate(variable333333, &vw_f_P);		VecDuplicate(variable333333, &vw_f_N);
	VecDuplicate(variable333333, &ww_f_P);		VecDuplicate(variable333333, &ww_f_N);

	//WENO---------------------------------------

	//速度绝对值的最大值.....
	VecDuplicate(variable333333, &variable_abs);

	//WENO_at_x..............

	//分量、权重和导数
	VecCreate(PETSC_COMM_WORLD, &d_variable_x_temp);
	VecSetSizes(d_variable_x_temp, PETSC_DECIDE, all_mesh_num100);
	VecSetFromOptions(d_variable_x_temp);

	VecDuplicate(d_variable_x_temp, &OmegaP0_at_x);	VecDuplicate(d_variable_x_temp, &OmegaN0_at_x);
	VecDuplicate(d_variable_x_temp, &OmegaP1_at_x);	VecDuplicate(d_variable_x_temp, &OmegaN1_at_x);
	VecDuplicate(d_variable_x_temp, &OmegaP2_at_x);	VecDuplicate(d_variable_x_temp, &OmegaN2_at_x);

	VecDuplicate(d_variable_x_temp, &f_P0_at_x);		VecDuplicate(d_variable_x_temp, &f_N0_at_x);
	VecDuplicate(d_variable_x_temp, &f_P1_at_x);		VecDuplicate(d_variable_x_temp, &f_N1_at_x);
	VecDuplicate(d_variable_x_temp, &f_P2_at_x);		VecDuplicate(d_variable_x_temp, &f_N2_at_x);

	VecDuplicate(d_variable_x_temp, &FP_at_x);		VecDuplicate(d_variable_x_temp, &FN_at_x);
	VecDuplicate(d_variable_x_temp, &FP_at_x0);		VecDuplicate(d_variable_x_temp, &FN_at_x0);
	VecDuplicate(d_variable_x_temp, &FP_at_x1);		VecDuplicate(d_variable_x_temp, &FN_at_x1);
	VecDuplicate(d_variable_x_temp, &FP_at_x2);		VecDuplicate(d_variable_x_temp, &FN_at_x2);

	//过渡量
	VecDuplicate(OmegaP0_at_x, &ISP0_part1_value_at_x);		VecDuplicate(OmegaP0_at_x, &ISP0_part2_value_at_x);
	VecDuplicate(OmegaP0_at_x, &ISP1_part1_value_at_x);		VecDuplicate(OmegaP0_at_x, &ISP1_part2_value_at_x);
	VecDuplicate(OmegaP0_at_x, &ISP2_part1_value_at_x);		VecDuplicate(OmegaP0_at_x, &ISP2_part2_value_at_x);

	VecDuplicate(OmegaP0_at_x, &sum_alphaP_at_x);
	VecDuplicate(OmegaP0_at_x, &alphaP0_at_x);
	VecDuplicate(OmegaP0_at_x, &alphaP1_at_x);
	VecDuplicate(OmegaP0_at_x, &alphaP2_at_x);

	VecDuplicate(OmegaN0_at_x, &ISN0_part1_value_at_x);		VecDuplicate(OmegaN0_at_x, &ISN0_part2_value_at_x);
	VecDuplicate(OmegaN0_at_x, &ISN1_part1_value_at_x);		VecDuplicate(OmegaN0_at_x, &ISN1_part2_value_at_x);
	VecDuplicate(OmegaN0_at_x, &ISN2_part1_value_at_x);		VecDuplicate(OmegaN0_at_x, &ISN2_part2_value_at_x);

	VecDuplicate(OmegaN0_at_x, &sum_alphaN_at_x);
	VecDuplicate(OmegaN0_at_x, &alphaN0_at_x);
	VecDuplicate(OmegaN0_at_x, &alphaN1_at_x);
	VecDuplicate(OmegaN0_at_x, &alphaN2_at_x);

	//WENO_at_y..............

	//分量、权重和导数
	VecCreate(PETSC_COMM_WORLD, &d_variable_y_temp);
	VecSetSizes(d_variable_y_temp, PETSC_DECIDE, all_mesh_num010);
	VecSetFromOptions(d_variable_y_temp);

	VecDuplicate(d_variable_y_temp, &OmegaP0_at_y);	VecDuplicate(d_variable_y_temp, &OmegaN0_at_y);
	VecDuplicate(d_variable_y_temp, &OmegaP1_at_y);	VecDuplicate(d_variable_y_temp, &OmegaN1_at_y);
	VecDuplicate(d_variable_y_temp, &OmegaP2_at_y);	VecDuplicate(d_variable_y_temp, &OmegaN2_at_y);

	VecDuplicate(d_variable_y_temp, &f_P0_at_y);		VecDuplicate(d_variable_y_temp, &f_N0_at_y);
	VecDuplicate(d_variable_y_temp, &f_P1_at_y);		VecDuplicate(d_variable_y_temp, &f_N1_at_y);
	VecDuplicate(d_variable_y_temp, &f_P2_at_y);		VecDuplicate(d_variable_y_temp, &f_N2_at_y);

	VecDuplicate(d_variable_y_temp, &FP_at_y);		VecDuplicate(d_variable_y_temp, &FN_at_y);
	VecDuplicate(d_variable_y_temp, &FP_at_y0);		VecDuplicate(d_variable_y_temp, &FN_at_y0);
	VecDuplicate(d_variable_y_temp, &FP_at_y1);		VecDuplicate(d_variable_y_temp, &FN_at_y1);
	VecDuplicate(d_variable_y_temp, &FP_at_y2);		VecDuplicate(d_variable_y_temp, &FN_at_y2);

	//过渡量
	VecDuplicate(OmegaP0_at_y, &ISP0_part1_value_at_y);		VecDuplicate(OmegaP0_at_y, &ISP0_part2_value_at_y);
	VecDuplicate(OmegaP0_at_y, &ISP1_part1_value_at_y);		VecDuplicate(OmegaP0_at_y, &ISP1_part2_value_at_y);
	VecDuplicate(OmegaP0_at_y, &ISP2_part1_value_at_y);		VecDuplicate(OmegaP0_at_y, &ISP2_part2_value_at_y);

	VecDuplicate(OmegaP0_at_y, &sum_alphaP_at_y);
	VecDuplicate(OmegaP0_at_y, &alphaP0_at_y);
	VecDuplicate(OmegaP0_at_y, &alphaP1_at_y);
	VecDuplicate(OmegaP0_at_y, &alphaP2_at_y);

	VecDuplicate(OmegaN0_at_y, &ISN0_part1_value_at_y);		VecDuplicate(OmegaN0_at_y, &ISN0_part2_value_at_y);
	VecDuplicate(OmegaN0_at_y, &ISN1_part1_value_at_y);		VecDuplicate(OmegaN0_at_y, &ISN1_part2_value_at_y);
	VecDuplicate(OmegaN0_at_y, &ISN2_part1_value_at_y);		VecDuplicate(OmegaN0_at_y, &ISN2_part2_value_at_y);

	VecDuplicate(OmegaN0_at_y, &sum_alphaN_at_y);
	VecDuplicate(OmegaN0_at_y, &alphaN0_at_y);
	VecDuplicate(OmegaN0_at_y, &alphaN1_at_y);
	VecDuplicate(OmegaN0_at_y, &alphaN2_at_y);

	//WENO_at_z..............

	//分量、权重和导数
	VecCreate(PETSC_COMM_WORLD, &d_variable_z_temp);
	VecSetSizes(d_variable_z_temp, PETSC_DECIDE, all_mesh_num001);
	VecSetFromOptions(d_variable_z_temp);

	VecDuplicate(d_variable_z_temp, &OmegaP0_at_z);	VecDuplicate(d_variable_z_temp, &OmegaN0_at_z);
	VecDuplicate(d_variable_z_temp, &OmegaP1_at_z);	VecDuplicate(d_variable_z_temp, &OmegaN1_at_z);
	VecDuplicate(d_variable_z_temp, &OmegaP2_at_z);	VecDuplicate(d_variable_z_temp, &OmegaN2_at_z);

	VecDuplicate(d_variable_z_temp, &f_P0_at_z);		VecDuplicate(d_variable_z_temp, &f_N0_at_z);
	VecDuplicate(d_variable_z_temp, &f_P1_at_z);		VecDuplicate(d_variable_z_temp, &f_N1_at_z);
	VecDuplicate(d_variable_z_temp, &f_P2_at_z);		VecDuplicate(d_variable_z_temp, &f_N2_at_z);

	VecDuplicate(d_variable_z_temp, &FP_at_z);		VecDuplicate(d_variable_z_temp, &FN_at_z);
	VecDuplicate(d_variable_z_temp, &FP_at_z0);		VecDuplicate(d_variable_z_temp, &FN_at_z0);
	VecDuplicate(d_variable_z_temp, &FP_at_z1);		VecDuplicate(d_variable_z_temp, &FN_at_z1);
	VecDuplicate(d_variable_z_temp, &FP_at_z2);		VecDuplicate(d_variable_z_temp, &FN_at_z2);

	//过渡量
	VecDuplicate(OmegaP0_at_z, &ISP0_part1_value_at_z);		VecDuplicate(OmegaP0_at_z, &ISP0_part2_value_at_z);
	VecDuplicate(OmegaP0_at_z, &ISP1_part1_value_at_z);		VecDuplicate(OmegaP0_at_z, &ISP1_part2_value_at_z);
	VecDuplicate(OmegaP0_at_z, &ISP2_part1_value_at_z);		VecDuplicate(OmegaP0_at_z, &ISP2_part2_value_at_z);

	VecDuplicate(OmegaP0_at_z, &sum_alphaP_at_z);
	VecDuplicate(OmegaP0_at_z, &alphaP0_at_z);
	VecDuplicate(OmegaP0_at_z, &alphaP1_at_z);
	VecDuplicate(OmegaP0_at_z, &alphaP2_at_z);

	VecDuplicate(OmegaN0_at_z, &ISN0_part1_value_at_z);		VecDuplicate(OmegaN0_at_z, &ISN0_part2_value_at_z);
	VecDuplicate(OmegaN0_at_z, &ISN1_part1_value_at_z);		VecDuplicate(OmegaN0_at_z, &ISN1_part2_value_at_z);
	VecDuplicate(OmegaN0_at_z, &ISN2_part1_value_at_z);		VecDuplicate(OmegaN0_at_z, &ISN2_part2_value_at_z);

	VecDuplicate(OmegaN0_at_z, &sum_alphaN_at_z);
	VecDuplicate(OmegaN0_at_z, &alphaN0_at_z);
	VecDuplicate(OmegaN0_at_z, &alphaN1_at_z);
	VecDuplicate(OmegaN0_at_z, &alphaN2_at_z);

	//TVD变量------------------------------------

	//phi....................
	VecDuplicate(phi, &TVD_phi1);
	VecDuplicate(phi, &TVD_phi2);
	VecDuplicate(phi, &TVD_phi_LSeq_right_term0);
	VecDuplicate(phi, &TVD_phi_LSeq_right_term1);
	VecDuplicate(phi, &TVD_phi_LSeq_right_term2);

	//FF.....................
	VecDuplicate(FF00, &TVD_FF00_1);		VecDuplicate(FF00, &TVD_FF00_2);
	VecDuplicate(FF00, &TVD_FF01_1);		VecDuplicate(FF00, &TVD_FF01_2);
	VecDuplicate(FF00, &TVD_FF02_1);		VecDuplicate(FF00, &TVD_FF02_2);
	VecDuplicate(FF00, &TVD_FF10_1);		VecDuplicate(FF00, &TVD_FF10_2);
	VecDuplicate(FF00, &TVD_FF11_1);		VecDuplicate(FF00, &TVD_FF11_2);
	VecDuplicate(FF00, &TVD_FF12_1);		VecDuplicate(FF00, &TVD_FF12_2);
	VecDuplicate(FF00, &TVD_FF20_1);		VecDuplicate(FF00, &TVD_FF20_2);
	VecDuplicate(FF00, &TVD_FF21_1);		VecDuplicate(FF00, &TVD_FF21_2);
	VecDuplicate(FF00, &TVD_FF22_1);		VecDuplicate(FF00, &TVD_FF22_2);

	VecDuplicate(FF00, &TVD_FF00_right_term0);	VecDuplicate(FF00, &TVD_FF00_right_term1);	VecDuplicate(FF00, &TVD_FF00_right_term2);
	VecDuplicate(FF00, &TVD_FF01_right_term0);	VecDuplicate(FF00, &TVD_FF01_right_term1);	VecDuplicate(FF00, &TVD_FF01_right_term2);
	VecDuplicate(FF00, &TVD_FF02_right_term0);	VecDuplicate(FF00, &TVD_FF02_right_term1);	VecDuplicate(FF00, &TVD_FF02_right_term2);
	VecDuplicate(FF00, &TVD_FF10_right_term0);	VecDuplicate(FF00, &TVD_FF10_right_term1);	VecDuplicate(FF00, &TVD_FF10_right_term2);
	VecDuplicate(FF00, &TVD_FF11_right_term0);	VecDuplicate(FF00, &TVD_FF11_right_term1);	VecDuplicate(FF00, &TVD_FF11_right_term2);
	VecDuplicate(FF00, &TVD_FF12_right_term0);	VecDuplicate(FF00, &TVD_FF12_right_term1);	VecDuplicate(FF00, &TVD_FF12_right_term2);
	VecDuplicate(FF00, &TVD_FF20_right_term0);	VecDuplicate(FF00, &TVD_FF20_right_term1);	VecDuplicate(FF00, &TVD_FF20_right_term2);
	VecDuplicate(FF00, &TVD_FF21_right_term0);	VecDuplicate(FF00, &TVD_FF21_right_term1);	VecDuplicate(FF00, &TVD_FF21_right_term2);
	VecDuplicate(FF00, &TVD_FF22_right_term0);	VecDuplicate(FF00, &TVD_FF22_right_term1);	VecDuplicate(FF00, &TVD_FF22_right_term2);

	//FF_inv.................
	VecDuplicate(FF_inv00, &TVD_FF_inv00_1);		VecDuplicate(FF_inv00, &TVD_FF_inv00_2);
	VecDuplicate(FF_inv00, &TVD_FF_inv01_1);		VecDuplicate(FF_inv00, &TVD_FF_inv01_2);
	VecDuplicate(FF_inv00, &TVD_FF_inv02_1);		VecDuplicate(FF_inv00, &TVD_FF_inv02_2);
	VecDuplicate(FF_inv00, &TVD_FF_inv10_1);		VecDuplicate(FF_inv00, &TVD_FF_inv10_2);
	VecDuplicate(FF_inv00, &TVD_FF_inv11_1);		VecDuplicate(FF_inv00, &TVD_FF_inv11_2);
	VecDuplicate(FF_inv00, &TVD_FF_inv12_1);		VecDuplicate(FF_inv00, &TVD_FF_inv12_2);
	VecDuplicate(FF_inv00, &TVD_FF_inv20_1);		VecDuplicate(FF_inv00, &TVD_FF_inv20_2);
	VecDuplicate(FF_inv00, &TVD_FF_inv21_1);		VecDuplicate(FF_inv00, &TVD_FF_inv21_2);
	VecDuplicate(FF_inv00, &TVD_FF_inv22_1);		VecDuplicate(FF_inv00, &TVD_FF_inv22_2);

	VecDuplicate(FF_inv00, &TVD_FF_inv00_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv00_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv00_right_term2);
	VecDuplicate(FF_inv00, &TVD_FF_inv01_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv01_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv01_right_term2);
	VecDuplicate(FF_inv00, &TVD_FF_inv02_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv02_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv02_right_term2);
	VecDuplicate(FF_inv00, &TVD_FF_inv10_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv10_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv10_right_term2);
	VecDuplicate(FF_inv00, &TVD_FF_inv11_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv11_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv11_right_term2);
	VecDuplicate(FF_inv00, &TVD_FF_inv12_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv12_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv12_right_term2);
	VecDuplicate(FF_inv00, &TVD_FF_inv20_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv20_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv20_right_term2);
	VecDuplicate(FF_inv00, &TVD_FF_inv21_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv21_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv21_right_term2);
	VecDuplicate(FF_inv00, &TVD_FF_inv22_right_term0);	VecDuplicate(FF_inv00, &TVD_FF_inv22_right_term1);	VecDuplicate(FF_inv00, &TVD_FF_inv22_right_term2);

	//矩阵变量=====================================================================================

	//导数相关矩阵-------------------------------
	MatCreate(PETSC_COMM_WORLD, &dx_matrix100000_or_010000_to_000000);
	MatSetSizes(dx_matrix100000_or_010000_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num100);
	MatSetFromOptions(dx_matrix100000_or_010000_to_000000);

	MatCreate(PETSC_COMM_WORLD, &dy_matrix001000_or_000100_to_000000);
	MatSetSizes(dy_matrix001000_or_000100_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num010);
	MatSetFromOptions(dy_matrix001000_or_000100_to_000000);

	MatCreate(PETSC_COMM_WORLD, &dz_matrix000010_or_000001_to_000000);
	MatSetSizes(dz_matrix000010_or_000001_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num001);
	MatSetFromOptions(dz_matrix000010_or_000001_to_000000);

	MatCreate(PETSC_COMM_WORLD, &dx_matrix111111_to_000000);
	MatSetSizes(dx_matrix111111_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num222);
	MatSetFromOptions(dx_matrix111111_to_000000);

	MatCreate(PETSC_COMM_WORLD, &dy_matrix111111_to_000000);
	MatSetSizes(dy_matrix111111_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num222);
	MatSetFromOptions(dy_matrix111111_to_000000);

	MatCreate(PETSC_COMM_WORLD, &dz_matrix111111_to_000000);
	MatSetSizes(dz_matrix111111_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num222);
	MatSetFromOptions(dz_matrix111111_to_000000);

	MatCreate(PETSC_COMM_WORLD, &ddx_matrix111111_to_000000);
	MatSetSizes(ddx_matrix111111_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num222);
	MatSetFromOptions(ddx_matrix111111_to_000000);

	MatCreate(PETSC_COMM_WORLD, &ddy_matrix111111_to_000000);
	MatSetSizes(ddy_matrix111111_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num222);
	MatSetFromOptions(ddy_matrix111111_to_000000);

	MatCreate(PETSC_COMM_WORLD, &ddz_matrix111111_to_000000);
	MatSetSizes(ddz_matrix111111_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num222);
	MatSetFromOptions(ddz_matrix111111_to_000000);

	//平均数据矩阵-------------------------------
	MatCreate(PETSC_COMM_WORLD, &average_matrix100000_or_010000_to_000000);
	MatSetSizes(average_matrix100000_or_010000_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num100);
	MatSetFromOptions(average_matrix100000_or_010000_to_000000);

	MatCreate(PETSC_COMM_WORLD, &average_matrix001000_or_000100_to_000000);
	MatSetSizes(average_matrix001000_or_000100_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num010);
	MatSetFromOptions(average_matrix001000_or_000100_to_000000);

	MatCreate(PETSC_COMM_WORLD, &average_matrix000010_or_000001_to_000000);
	MatSetSizes(average_matrix000010_or_000001_to_000000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num, all_mesh_num001);
	MatSetFromOptions(average_matrix000010_or_000001_to_000000);

	MatCreate(PETSC_COMM_WORLD, &average_matrix344333_or_433433_to_333333);
	MatSetSizes(average_matrix344333_or_433433_to_333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num776);
	MatSetFromOptions(average_matrix344333_or_433433_to_333333);

	MatCreate(PETSC_COMM_WORLD, &average_matrix343343_or_433334_to_333333);
	MatSetSizes(average_matrix343343_or_433334_to_333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num767);
	MatSetFromOptions(average_matrix343343_or_433334_to_333333);

	MatCreate(PETSC_COMM_WORLD, &average_matrix333443_or_334334_to_333333);
	MatSetSizes(average_matrix333443_or_334334_to_333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num677);
	MatSetFromOptions(average_matrix333443_or_334334_to_333333);

	//拓展矩阵-----------------------------------

	//速度

	//nonstagger_extend_u_matrix111111
	MatCreate(PETSC_COMM_WORLD, &nonstagger_extend_u_matrix111111);
	MatSetSizes(nonstagger_extend_u_matrix111111, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num222, all_mesh_num);
	MatSetFromOptions(nonstagger_extend_u_matrix111111);

	//nonstagger_extend_v_matrix111111
	MatCreate(PETSC_COMM_WORLD, &nonstagger_extend_v_matrix111111);
	MatSetSizes(nonstagger_extend_v_matrix111111, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num222, all_mesh_num);
	MatSetFromOptions(nonstagger_extend_v_matrix111111);

	//nonstagger_extend_w_matrix111111
	MatCreate(PETSC_COMM_WORLD, &nonstagger_extend_w_matrix111111);
	MatSetSizes(nonstagger_extend_w_matrix111111, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num222, all_mesh_num);
	MatSetFromOptions(nonstagger_extend_w_matrix111111);

	//nonstagger_extend_u_matrix333333
	MatCreate(PETSC_COMM_WORLD, &nonstagger_extend_u_matrix333333);
	MatSetSizes(nonstagger_extend_u_matrix333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num);
	MatSetFromOptions(nonstagger_extend_u_matrix333333);

	//nonstagger_extend_v_matrix333333
	MatCreate(PETSC_COMM_WORLD, &nonstagger_extend_v_matrix333333);
	MatSetSizes(nonstagger_extend_v_matrix333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num);
	MatSetFromOptions(nonstagger_extend_v_matrix333333);

	//nonstagger_extend_w_matrix333333
	MatCreate(PETSC_COMM_WORLD, &nonstagger_extend_w_matrix333333);
	MatSetSizes(nonstagger_extend_w_matrix333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num);
	MatSetFromOptions(nonstagger_extend_w_matrix333333);

	//extend_u_matrix010000
	MatCreate(PETSC_COMM_WORLD, &extend_u_matrix010000);
	MatSetSizes(extend_u_matrix010000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num);
	MatSetFromOptions(extend_u_matrix010000);

	//extend_v_matrix000100
	MatCreate(PETSC_COMM_WORLD, &extend_v_matrix000100);
	MatSetSizes(extend_v_matrix000100, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num);
	MatSetFromOptions(extend_v_matrix000100);

	//extend_w_matrix000001
	MatCreate(PETSC_COMM_WORLD, &extend_w_matrix000001);
	MatSetSizes(extend_w_matrix000001, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num);
	MatSetFromOptions(extend_w_matrix000001);

	//extend_u_matrix111111 and extend_v_matrix111111 and extend_w_matrix111111
	MatCreate(PETSC_COMM_WORLD, &extend_u_matrix111111);
	MatSetSizes(extend_u_matrix111111, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num222, all_mesh_num);
	MatSetFromOptions(extend_u_matrix111111);

	MatCreate(PETSC_COMM_WORLD, &extend_v_matrix111111);
	MatSetSizes(extend_v_matrix111111, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num222, all_mesh_num);
	MatSetFromOptions(extend_v_matrix111111);

	MatCreate(PETSC_COMM_WORLD, &extend_w_matrix111111);
	MatSetSizes(extend_w_matrix111111, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num222, all_mesh_num);
	MatSetFromOptions(extend_w_matrix111111);

	//extend_u_matrix333333 and extend_v_matrix333333 and extend_w_matrix333333
	MatCreate(PETSC_COMM_WORLD, &extend_u_matrix333333);
	MatSetSizes(extend_u_matrix333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num);
	MatSetFromOptions(extend_u_matrix333333);

	MatCreate(PETSC_COMM_WORLD, &extend_v_matrix333333);
	MatSetSizes(extend_v_matrix333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num);
	MatSetFromOptions(extend_v_matrix333333);

	MatCreate(PETSC_COMM_WORLD, &extend_w_matrix333333);
	MatSetSizes(extend_w_matrix333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num);
	MatSetFromOptions(extend_w_matrix333333);

	//extend_u_matrix344333 and extend_u_matrix343343
	MatCreate(PETSC_COMM_WORLD, &extend_u_matrix344333);
	MatSetSizes(extend_u_matrix344333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num776, all_mesh_num);
	MatSetFromOptions(extend_u_matrix344333);

	MatCreate(PETSC_COMM_WORLD, &extend_u_matrix343343);
	MatSetSizes(extend_u_matrix343343, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num767, all_mesh_num);
	MatSetFromOptions(extend_u_matrix343343);

	//extend_v_matrix433433 and extend_v_matrix333443
	MatCreate(PETSC_COMM_WORLD, &extend_v_matrix433433);
	MatSetSizes(extend_v_matrix433433, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num776, all_mesh_num);
	MatSetFromOptions(extend_v_matrix433433);

	MatCreate(PETSC_COMM_WORLD, &extend_v_matrix333443);
	MatSetSizes(extend_v_matrix333443, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num677, all_mesh_num);
	MatSetFromOptions(extend_v_matrix333443);

	//extend_w_matrix433334 and extend_w_matrix334334
	MatCreate(PETSC_COMM_WORLD, &extend_w_matrix433334);
	MatSetSizes(extend_w_matrix433334, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num767, all_mesh_num);
	MatSetFromOptions(extend_w_matrix433334);

	MatCreate(PETSC_COMM_WORLD, &extend_w_matrix334334);
	MatSetSizes(extend_w_matrix334334, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num677, all_mesh_num);
	MatSetFromOptions(extend_w_matrix334334);

	//u_auxiliary_vector010000
	VecCreate(PETSC_COMM_WORLD, &u_auxiliary_vector010000);
	VecSetSizes(u_auxiliary_vector010000, PETSC_DECIDE, all_mesh_num100);
	VecSetFromOptions(u_auxiliary_vector010000);

	//v_auxiliary_vector000100
	//VecCreate(PETSC_COMM_WORLD, &v_auxiliary_vector000100);
	//VecSetSizes(v_auxiliary_vector000100, PETSC_DECIDE, all_mesh_num010);
	//VecSetFromOptions(v_auxiliary_vector000100);

	//w_auxiliary_vector000100
	//VecCreate(PETSC_COMM_WORLD, &w_auxiliary_vector000001);
	//VecSetSizes(w_auxiliary_vector000001, PETSC_DECIDE, all_mesh_num001);
	//VecSetFromOptions(w_auxiliary_vector000001);

	//u_auxiliary_vector111111
	VecCreate(PETSC_COMM_WORLD, &u_auxiliary_vector111111);
	VecSetSizes(u_auxiliary_vector111111, PETSC_DECIDE, all_mesh_num222);
	VecSetFromOptions(u_auxiliary_vector111111);

	//v_auxiliary_vector111111
	//VecCreate(PETSC_COMM_WORLD, &v_auxiliary_vector111111);
	//VecSetSizes(v_auxiliary_vector111111, PETSC_DECIDE, all_mesh_num222);
	//VecSetFromOptions(v_auxiliary_vector111111);

	//w_auxiliary_vector111111
	//VecCreate(PETSC_COMM_WORLD, &w_auxiliary_vector111111);
	//VecSetSizes(w_auxiliary_vector111111, PETSC_DECIDE, all_mesh_num222);
	//VecSetFromOptions(w_auxiliary_vector111111);

	//u_auxiliary_vector333333
	VecCreate(PETSC_COMM_WORLD, &u_auxiliary_vector333333);
	VecSetSizes(u_auxiliary_vector333333, PETSC_DECIDE, all_mesh_num666);
	VecSetFromOptions(u_auxiliary_vector333333);

	//v_auxiliary_vector333333
	//VecCreate(PETSC_COMM_WORLD, &v_auxiliary_vector333333);
	//VecSetSizes(v_auxiliary_vector333333, PETSC_DECIDE, all_mesh_num666);
	//VecSetFromOptions(v_auxiliary_vector333333);

	//w_auxiliary_vector333333
	//VecCreate(PETSC_COMM_WORLD, &w_auxiliary_vector333333);
	//VecSetSizes(w_auxiliary_vector333333, PETSC_DECIDE, all_mesh_num666);
	//VecSetFromOptions(w_auxiliary_vector333333);

	//u_auxiliary_vector344333 and u_auxiliary_vector343343
	VecCreate(PETSC_COMM_WORLD, &u_auxiliary_vector344333);
	VecSetSizes(u_auxiliary_vector344333, PETSC_DECIDE, all_mesh_num776);
	VecSetFromOptions(u_auxiliary_vector344333);

	VecCreate(PETSC_COMM_WORLD, &u_auxiliary_vector343343);
	VecSetSizes(u_auxiliary_vector343343, PETSC_DECIDE, all_mesh_num767);
	VecSetFromOptions(u_auxiliary_vector343343);

	//v_auxiliary_vector433433 and v_auxiliary_vector333443
	//VecCreate(PETSC_COMM_WORLD, &v_auxiliary_vector433433);
	//VecSetSizes(v_auxiliary_vector433433, PETSC_DECIDE, all_mesh_num776);
	//VecSetFromOptions(v_auxiliary_vector433433);

	//VecCreate(PETSC_COMM_WORLD, &v_auxiliary_vector333443);
	//VecSetSizes(v_auxiliary_vector333443, PETSC_DECIDE, all_mesh_num677);
	//VecSetFromOptions(v_auxiliary_vector333443);

	//w_auxiliary_vector433334 and w_auxiliary_vector334334
	//VecCreate(PETSC_COMM_WORLD, &w_auxiliary_vector433334);
	//VecSetSizes(w_auxiliary_vector433334, PETSC_DECIDE, all_mesh_num767);
	//VecSetFromOptions(w_auxiliary_vector433334);

	//VecCreate(PETSC_COMM_WORLD, &w_auxiliary_vector334334);
	//VecSetSizes(w_auxiliary_vector334334, PETSC_DECIDE, all_mesh_num677);
	//VecSetFromOptions(w_auxiliary_vector334334);

	//pressure...............
	MatCreate(PETSC_COMM_WORLD, &extend_pressure_matrix100000);
	MatSetSizes(extend_pressure_matrix100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num);
	MatSetFromOptions(extend_pressure_matrix100000);

	MatCreate(PETSC_COMM_WORLD, &extend_pressure_matrix001000);
	MatSetSizes(extend_pressure_matrix001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num);
	MatSetFromOptions(extend_pressure_matrix001000);

	MatCreate(PETSC_COMM_WORLD, &extend_pressure_matrix000010);
	MatSetSizes(extend_pressure_matrix000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num);
	MatSetFromOptions(extend_pressure_matrix000010);

	MatCreate(PETSC_COMM_WORLD, &extend_pressure_matrix111111);
	MatSetSizes(extend_pressure_matrix111111, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num222, all_mesh_num);
	MatSetFromOptions(extend_pressure_matrix111111);

	//VecDuplicate(u_auxiliary_vector010000, &pressure_auxiliary_vector100000);
	//VecDuplicate(v_auxiliary_vector000100, &pressure_auxiliary_vector001000);
	//VecDuplicate(w_auxiliary_vector000001, &pressure_auxiliary_vector000010);
	//VecDuplicate(u_auxiliary_vector111111, &pressure_auxiliary_vector111111);

	//phi....................
	MatCreate(PETSC_COMM_WORLD, &extend_phi_matrix100000);
	MatSetSizes(extend_phi_matrix100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num);
	MatSetFromOptions(extend_phi_matrix100000);

	MatCreate(PETSC_COMM_WORLD, &extend_phi_matrix001000);
	MatSetSizes(extend_phi_matrix001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num);
	MatSetFromOptions(extend_phi_matrix001000);

	MatCreate(PETSC_COMM_WORLD, &extend_phi_matrix000010);
	MatSetSizes(extend_phi_matrix000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num);
	MatSetFromOptions(extend_phi_matrix000010);

	MatCreate(PETSC_COMM_WORLD, &extend_phi_matrix111111);
	MatSetSizes(extend_phi_matrix111111, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num222, all_mesh_num);
	MatSetFromOptions(extend_phi_matrix111111);

	MatCreate(PETSC_COMM_WORLD, &extend_phi_matrix333333);
	MatSetSizes(extend_phi_matrix333333, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num666, all_mesh_num);
	MatSetFromOptions(extend_phi_matrix333333);

	//VecDuplicate(u_auxiliary_vector010000, &phi_auxiliary_vector100000);
	//VecDuplicate(v_auxiliary_vector000100, &phi_auxiliary_vector001000);
	//VecDuplicate(w_auxiliary_vector000001, &phi_auxiliary_vector000010);
	//VecDuplicate(u_auxiliary_vector111111, &phi_auxiliary_vector111111);
	//VecDuplicate(u_auxiliary_vector333333, &phi_auxiliary_vector333333);

	//WENO 相关矩阵------------------------------

	//x

	//f_P
	MatCreate(PETSC_COMM_WORLD, &f_P0_x_matrix333333_to_100000);
	MatSetSizes(f_P0_x_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(f_P0_x_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &f_P1_x_matrix333333_to_100000);
	MatSetSizes(f_P1_x_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(f_P1_x_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &f_P2_x_matrix333333_to_100000);
	MatSetSizes(f_P2_x_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(f_P2_x_matrix333333_to_100000);

	//f_N
	MatCreate(PETSC_COMM_WORLD, &f_N0_x_matrix333333_to_100000);
	MatSetSizes(f_N0_x_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(f_N0_x_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &f_N1_x_matrix333333_to_100000);
	MatSetSizes(f_N1_x_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(f_N1_x_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &f_N2_x_matrix333333_to_100000);
	MatSetSizes(f_N2_x_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(f_N2_x_matrix333333_to_100000);

	//ISP part1
	MatCreate(PETSC_COMM_WORLD, &ISP0_x_part1_matrix333333_to_100000);
	MatSetSizes(ISP0_x_part1_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISP0_x_part1_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &ISP1_x_part1_matrix333333_to_100000);
	MatSetSizes(ISP1_x_part1_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISP1_x_part1_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &ISP2_x_part1_matrix333333_to_100000);
	MatSetSizes(ISP2_x_part1_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISP2_x_part1_matrix333333_to_100000);

	//ISP part2
	MatCreate(PETSC_COMM_WORLD, &ISP0_x_part2_matrix333333_to_100000);
	MatSetSizes(ISP0_x_part2_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISP0_x_part2_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &ISP1_x_part2_matrix333333_to_100000);
	MatSetSizes(ISP1_x_part2_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISP1_x_part2_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &ISP2_x_part2_matrix333333_to_100000);
	MatSetSizes(ISP2_x_part2_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISP2_x_part2_matrix333333_to_100000);

	//ISN part1
	MatCreate(PETSC_COMM_WORLD, &ISN0_x_part1_matrix333333_to_100000);
	MatSetSizes(ISN0_x_part1_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISN0_x_part1_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &ISN1_x_part1_matrix333333_to_100000);
	MatSetSizes(ISN1_x_part1_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISN1_x_part1_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &ISN2_x_part1_matrix333333_to_100000);
	MatSetSizes(ISN2_x_part1_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISN2_x_part1_matrix333333_to_100000);

	//ISN part2
	MatCreate(PETSC_COMM_WORLD, &ISN0_x_part2_matrix333333_to_100000);
	MatSetSizes(ISN0_x_part2_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISN0_x_part2_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &ISN1_x_part2_matrix333333_to_100000);
	MatSetSizes(ISN1_x_part2_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISN1_x_part2_matrix333333_to_100000);

	MatCreate(PETSC_COMM_WORLD, &ISN2_x_part2_matrix333333_to_100000);
	MatSetSizes(ISN2_x_part2_matrix333333_to_100000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num100, all_mesh_num666);
	MatSetFromOptions(ISN2_x_part2_matrix333333_to_100000);

	//y

	//f_P
	MatCreate(PETSC_COMM_WORLD, &f_P0_y_matrix333333_to_001000);
	MatSetSizes(f_P0_y_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(f_P0_y_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &f_P1_y_matrix333333_to_001000);
	MatSetSizes(f_P1_y_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(f_P1_y_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &f_P2_y_matrix333333_to_001000);
	MatSetSizes(f_P2_y_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(f_P2_y_matrix333333_to_001000);

	//f_N
	MatCreate(PETSC_COMM_WORLD, &f_N0_y_matrix333333_to_001000);
	MatSetSizes(f_N0_y_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(f_N0_y_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &f_N1_y_matrix333333_to_001000);
	MatSetSizes(f_N1_y_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(f_N1_y_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &f_N2_y_matrix333333_to_001000);
	MatSetSizes(f_N2_y_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(f_N2_y_matrix333333_to_001000);

	//ISP part1
	MatCreate(PETSC_COMM_WORLD, &ISP0_y_part1_matrix333333_to_001000);
	MatSetSizes(ISP0_y_part1_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISP0_y_part1_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &ISP1_y_part1_matrix333333_to_001000);
	MatSetSizes(ISP1_y_part1_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISP1_y_part1_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &ISP2_y_part1_matrix333333_to_001000);
	MatSetSizes(ISP2_y_part1_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISP2_y_part1_matrix333333_to_001000);

	//ISP part2
	MatCreate(PETSC_COMM_WORLD, &ISP0_y_part2_matrix333333_to_001000);
	MatSetSizes(ISP0_y_part2_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISP0_y_part2_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &ISP1_y_part2_matrix333333_to_001000);
	MatSetSizes(ISP1_y_part2_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISP1_y_part2_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &ISP2_y_part2_matrix333333_to_001000);
	MatSetSizes(ISP2_y_part2_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISP2_y_part2_matrix333333_to_001000);

	//ISN part1
	MatCreate(PETSC_COMM_WORLD, &ISN0_y_part1_matrix333333_to_001000);
	MatSetSizes(ISN0_y_part1_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISN0_y_part1_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &ISN1_y_part1_matrix333333_to_001000);
	MatSetSizes(ISN1_y_part1_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISN1_y_part1_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &ISN2_y_part1_matrix333333_to_001000);
	MatSetSizes(ISN2_y_part1_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISN2_y_part1_matrix333333_to_001000);

	//ISN part2
	MatCreate(PETSC_COMM_WORLD, &ISN0_y_part2_matrix333333_to_001000);
	MatSetSizes(ISN0_y_part2_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISN0_y_part2_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &ISN1_y_part2_matrix333333_to_001000);
	MatSetSizes(ISN1_y_part2_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISN1_y_part2_matrix333333_to_001000);

	MatCreate(PETSC_COMM_WORLD, &ISN2_y_part2_matrix333333_to_001000);
	MatSetSizes(ISN2_y_part2_matrix333333_to_001000, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num010, all_mesh_num666);
	MatSetFromOptions(ISN2_y_part2_matrix333333_to_001000);

	//z

	//f_P
	MatCreate(PETSC_COMM_WORLD, &f_P0_z_matrix333333_to_000010);
	MatSetSizes(f_P0_z_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(f_P0_z_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &f_P1_z_matrix333333_to_000010);
	MatSetSizes(f_P1_z_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(f_P1_z_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &f_P2_z_matrix333333_to_000010);
	MatSetSizes(f_P2_z_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(f_P2_z_matrix333333_to_000010);

	//f_N
	MatCreate(PETSC_COMM_WORLD, &f_N0_z_matrix333333_to_000010);
	MatSetSizes(f_N0_z_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(f_N0_z_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &f_N1_z_matrix333333_to_000010);
	MatSetSizes(f_N1_z_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(f_N1_z_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &f_N2_z_matrix333333_to_000010);
	MatSetSizes(f_N2_z_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(f_N2_z_matrix333333_to_000010);

	//ISP part1
	MatCreate(PETSC_COMM_WORLD, &ISP0_z_part1_matrix333333_to_000010);
	MatSetSizes(ISP0_z_part1_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISP0_z_part1_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &ISP1_z_part1_matrix333333_to_000010);
	MatSetSizes(ISP1_z_part1_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISP1_z_part1_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &ISP2_z_part1_matrix333333_to_000010);
	MatSetSizes(ISP2_z_part1_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISP2_z_part1_matrix333333_to_000010);

	//ISP part2
	MatCreate(PETSC_COMM_WORLD, &ISP0_z_part2_matrix333333_to_000010);
	MatSetSizes(ISP0_z_part2_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISP0_z_part2_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &ISP1_z_part2_matrix333333_to_000010);
	MatSetSizes(ISP1_z_part2_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISP1_z_part2_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &ISP2_z_part2_matrix333333_to_000010);
	MatSetSizes(ISP2_z_part2_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISP2_z_part2_matrix333333_to_000010);

	//ISN part1
	MatCreate(PETSC_COMM_WORLD, &ISN0_z_part1_matrix333333_to_000010);
	MatSetSizes(ISN0_z_part1_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISN0_z_part1_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &ISN1_z_part1_matrix333333_to_000010);
	MatSetSizes(ISN1_z_part1_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISN1_z_part1_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &ISN2_z_part1_matrix333333_to_000010);
	MatSetSizes(ISN2_z_part1_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISN2_z_part1_matrix333333_to_000010);

	//ISN part2
	MatCreate(PETSC_COMM_WORLD, &ISN0_z_part2_matrix333333_to_000010);
	MatSetSizes(ISN0_z_part2_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISN0_z_part2_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &ISN1_z_part2_matrix333333_to_000010);
	MatSetSizes(ISN1_z_part2_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISN1_z_part2_matrix333333_to_000010);

	MatCreate(PETSC_COMM_WORLD, &ISN2_z_part2_matrix333333_to_000010);
	MatSetSizes(ISN2_z_part2_matrix333333_to_000010, PETSC_DECIDE, PETSC_DECIDE, all_mesh_num001, all_mesh_num666);
	MatSetFromOptions(ISN2_z_part2_matrix333333_to_000010);
}

//变量赋值==========================================================================================
void wyh_solve_fluid_structure_interaction_system3D200630::prepare_assembly_matrix_and_vector()
{
	//基于实际问题的函数--------------------------

	//初始变量
	get_initial_condition();

	//方程组系数
	assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_u3D(NS_coeff_u3D);
	assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_v3D(NS_coeff_v3D);
	assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_w3D(NS_coeff_w3D);
	assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_p3D(NS_coeff_p3D);

	//KSP_NSeq_u系数矩阵
	KSPCreate(PETSC_COMM_WORLD, &KSP_NSeq_u);
	KSPSetOperators(KSP_NSeq_u, NS_coeff_u3D, NS_coeff_u3D);
	KSPSetUp(KSP_NSeq_u);

	KSPSetFromOptions(KSP_NSeq_u);
	KSPSetUp(KSP_NSeq_u);
	KSPSetType(KSP_NSeq_u, KSPBCGSL);
	KSPSetTolerances(KSP_NSeq_u, PETSC_DEFAULT, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT);

	//KSP_NSeq_v系数矩阵
	KSPCreate(PETSC_COMM_WORLD, &KSP_NSeq_v);
	KSPSetOperators(KSP_NSeq_v, NS_coeff_v3D, NS_coeff_v3D);
	KSPSetUp(KSP_NSeq_v);

	KSPSetFromOptions(KSP_NSeq_v);
	KSPSetUp(KSP_NSeq_v);
	KSPSetType(KSP_NSeq_v, KSPBCGSL);
	KSPSetTolerances(KSP_NSeq_v, PETSC_DEFAULT, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT);

	//KSP_NSeq_w系数矩阵
	KSPCreate(PETSC_COMM_WORLD, &KSP_NSeq_w);
	KSPSetOperators(KSP_NSeq_w, NS_coeff_w3D, NS_coeff_w3D);
	KSPSetUp(KSP_NSeq_w);

	KSPSetFromOptions(KSP_NSeq_w);
	KSPSetUp(KSP_NSeq_w);
	KSPSetType(KSP_NSeq_w, KSPBCGSL);
	KSPSetTolerances(KSP_NSeq_w, PETSC_DEFAULT, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT);

	//KSP_NSeq_p系数矩阵
	KSPCreate(PETSC_COMM_WORLD, &KSP_NSeq_p);
	KSPSetOperators(KSP_NSeq_p, NS_coeff_p3D, NS_coeff_p3D);
	KSPSetUp(KSP_NSeq_p);

	KSPSetFromOptions(KSP_NSeq_p);
	KSPSetUp(KSP_NSeq_p);
	KSPSetType(KSP_NSeq_p, KSPBCGSL);
	KSPSetTolerances(KSP_NSeq_p, PETSC_DEFAULT, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT);

	//修正方程右端项
	assembly_modified_u_NSeq_vector_for_lid_driven_cavity_flow_boundary_condition3D();

	//去掉边界的向量
	assembly_take_out_x_L_boundary_matrix3D();
	assembly_take_out_y_L_boundary_matrix3D();
	assembly_take_out_z_L_boundary_matrix3D();

	//添加边界的向量
	assembly_x_solution_out_x_L_extend_to_ori_size_matrix3D();
	assembly_y_solution_out_y_L_extend_to_ori_size_matrix3D();
	assembly_z_solution_out_z_L_extend_to_ori_size_matrix3D();
	
	//全局变量------------------------------------

	//WENO相关计算矩阵........

	//Cell boundary FP
	assembly_CellBoundaryFP_x_matrix333333_to_100000();
	assembly_CellBoundaryFP_y_matrix333333_to_001000();
	assembly_CellBoundaryFP_z_matrix333333_to_000010();

	//Cell boundary FN
	assembly_CellBoundaryFN_x_matrix333333_to_100000();
	assembly_CellBoundaryFN_y_matrix333333_to_001000();
	assembly_CellBoundaryFN_z_matrix333333_to_000010();

	//part 1 of the Weight FP
	assembly_part1_WeightFP_x_matrix333333_to_100000();
	assembly_part1_WeightFP_y_matrix333333_to_001000();
	assembly_part1_WeightFP_z_matrix333333_to_000010();

	//part 2 of the Weight FP
	assembly_part2_WeightFP_x_matrix333333_to_100000();
	assembly_part2_WeightFP_y_matrix333333_to_001000();
	assembly_part2_WeightFP_z_matrix333333_to_000010();

	//part 1 of the Weight FN
	assembly_part1_WeightFN_x_matrix333333_to_100000();
	assembly_part1_WeightFN_y_matrix333333_to_001000();
	assembly_part1_WeightFN_z_matrix333333_to_000010();

	//part 2 of the Weight FN
	assembly_part2_WeightFN_x_matrix333333_to_100000();
	assembly_part2_WeightFN_y_matrix333333_to_001000();
	assembly_part2_WeightFN_z_matrix333333_to_000010();

	//求导相关矩阵............

	assembly_dx_matrix100000_or_010000_to_000000();
	assembly_dy_matrix001000_or_000100_to_000000();
	assembly_dz_matrix000010_or_000001_to_000000();
	assembly_dx_matrix111111_to_000000();
	assembly_dy_matrix111111_to_000000();
	assembly_dz_matrix111111_to_000000();

	//2阶中心差分
	assembly_ddx_matrix111111_to_000000();
	assembly_ddy_matrix111111_to_000000();
	assembly_ddz_matrix111111_to_000000();

	//平均数据................

	//将(100000或010000类型数据求两点平均000000)
	assembly_average_matrix100000_or_010000_to_000000();

	//将(001000或000100类型数据求两点平均000000)
	assembly_average_matrix001000_or_000100_to_000000();

	//将(000010或000001类型数据求两点平均000000)
	assembly_average_matrix000010_or_000001_to_000000();

	//组装由交错网格处速度获得非交错网格速度
	assembly_average_matrix344333_or_433433_to_333333();
	assembly_average_matrix343343_or_433334_to_333333();
	assembly_average_matrix333443_or_334334_to_333333();

	//拓展数据点..............

	//速度
	assembly_nonstagger_extend_velocity_u_matrix3D(nonstagger_extend_u_matrix111111, extend_num111111, velocity_boundary);
	assembly_nonstagger_extend_velocity_v_matrix3D(nonstagger_extend_v_matrix111111, extend_num111111, velocity_boundary);
	assembly_nonstagger_extend_velocity_w_matrix3D(nonstagger_extend_w_matrix111111, extend_num111111, velocity_boundary);

	assembly_nonstagger_extend_velocity_u_matrix3D(nonstagger_extend_u_matrix333333, extend_num333333, velocity_boundary);
	assembly_nonstagger_extend_velocity_v_matrix3D(nonstagger_extend_v_matrix333333, extend_num333333, velocity_boundary);
	assembly_nonstagger_extend_velocity_w_matrix3D(nonstagger_extend_w_matrix333333, extend_num333333, velocity_boundary);

	assembly_stagger_extend_velocity_u_matrix3D(extend_u_matrix010000, extend_num010000, velocity_boundary);
	assembly_stagger_extend_velocity_u_matrix3D(extend_u_matrix111111, extend_num111111, velocity_boundary);
	assembly_stagger_extend_velocity_u_matrix3D(extend_u_matrix333333, extend_num333333, velocity_boundary);
	assembly_stagger_extend_velocity_u_matrix3D(extend_u_matrix344333, extend_num344333, velocity_boundary);
	assembly_stagger_extend_velocity_u_matrix3D(extend_u_matrix343343, extend_num343343, velocity_boundary);

	assembly_stagger_extend_velocity_v_matrix3D(extend_v_matrix000100, extend_num000100, velocity_boundary);
	assembly_stagger_extend_velocity_v_matrix3D(extend_v_matrix111111, extend_num111111, velocity_boundary);
	assembly_stagger_extend_velocity_v_matrix3D(extend_v_matrix333333, extend_num333333, velocity_boundary);
	assembly_stagger_extend_velocity_v_matrix3D(extend_v_matrix433433, extend_num433433, velocity_boundary);
	assembly_stagger_extend_velocity_v_matrix3D(extend_v_matrix333443, extend_num333443, velocity_boundary);

	assembly_stagger_extend_velocity_w_matrix3D(extend_w_matrix000001, extend_num000001, velocity_boundary);
	assembly_stagger_extend_velocity_w_matrix3D(extend_w_matrix111111, extend_num111111, velocity_boundary);
	assembly_stagger_extend_velocity_w_matrix3D(extend_w_matrix333333, extend_num333333, velocity_boundary);
	assembly_stagger_extend_velocity_w_matrix3D(extend_w_matrix433334, extend_num433334, velocity_boundary);
	assembly_stagger_extend_velocity_w_matrix3D(extend_w_matrix334334, extend_num334334, velocity_boundary);

	assembly_extend_velocity_u_auxiliary_vector3D(u_auxiliary_vector010000, extend_num010000, velocity_boundary);
	assembly_extend_velocity_u_auxiliary_vector3D(u_auxiliary_vector111111, extend_num111111, velocity_boundary);
	assembly_extend_velocity_u_auxiliary_vector3D(u_auxiliary_vector333333, extend_num333333, velocity_boundary);
	assembly_extend_velocity_u_auxiliary_vector3D(u_auxiliary_vector344333, extend_num344333, velocity_boundary);
	assembly_extend_velocity_u_auxiliary_vector3D(u_auxiliary_vector343343, extend_num343343, velocity_boundary);

	//assembly_extend_velocity_v_auxiliary_vector3D(v_auxiliary_vector000100, extend_num000100, velocity_boundary);
	//assembly_extend_velocity_v_auxiliary_vector3D(v_auxiliary_vector111111, extend_num111111, velocity_boundary);
	//assembly_extend_velocity_v_auxiliary_vector3D(v_auxiliary_vector333333, extend_num333333, velocity_boundary);
	//assembly_extend_velocity_v_auxiliary_vector3D(v_auxiliary_vector433433, extend_num433433, velocity_boundary);
	//assembly_extend_velocity_v_auxiliary_vector3D(v_auxiliary_vector333443, extend_num333443, velocity_boundary);

	//assembly_extend_velocity_w_auxiliary_vector3D(w_auxiliary_vector000001, extend_num000001, velocity_boundary);
	//assembly_extend_velocity_w_auxiliary_vector3D(w_auxiliary_vector111111, extend_num111111, velocity_boundary);
	//assembly_extend_velocity_w_auxiliary_vector3D(w_auxiliary_vector333333, extend_num333333, velocity_boundary);
	//assembly_extend_velocity_w_auxiliary_vector3D(w_auxiliary_vector433334, extend_num433334, velocity_boundary);
	//assembly_extend_velocity_w_auxiliary_vector3D(w_auxiliary_vector334334, extend_num334334, velocity_boundary);

	//压力
	assembly_nonstagger_extend_pressure_matrix3D(extend_pressure_matrix100000, extend_num100000, pressure_boundary);
	assembly_nonstagger_extend_pressure_matrix3D(extend_pressure_matrix001000, extend_num001000, pressure_boundary);
	assembly_nonstagger_extend_pressure_matrix3D(extend_pressure_matrix000010, extend_num000010, pressure_boundary);
	assembly_nonstagger_extend_pressure_matrix3D(extend_pressure_matrix111111, extend_num111111, pressure_boundary);

	//assembly_extend_pressure_auxiliary_vector3D(pressure_auxiliary_vector100000, extend_num100000, pressure_boundary);
	//assembly_extend_pressure_auxiliary_vector3D(pressure_auxiliary_vector001000, extend_num001000, pressure_boundary);
	//assembly_extend_pressure_auxiliary_vector3D(pressure_auxiliary_vector000010, extend_num000010, pressure_boundary);
	//assembly_extend_pressure_auxiliary_vector3D(pressure_auxiliary_vector111111, extend_num111111, pressure_boundary);
	
	//phi
	assembly_nonstagger_extend_phi_matrix3D(extend_phi_matrix100000, extend_num100000, phi_boundary);
	assembly_nonstagger_extend_phi_matrix3D(extend_phi_matrix001000, extend_num001000, phi_boundary);
	assembly_nonstagger_extend_phi_matrix3D(extend_phi_matrix000010, extend_num000010, phi_boundary);
	assembly_nonstagger_extend_phi_matrix3D(extend_phi_matrix111111, extend_num111111, phi_boundary);
	assembly_nonstagger_extend_phi_matrix3D(extend_phi_matrix333333, extend_num333333, phi_boundary);

	//assembly_extend_phi_auxiliary_vector3D(phi_auxiliary_vector100000, extend_num100000, phi_boundary);
	//assembly_extend_phi_auxiliary_vector3D(phi_auxiliary_vector001000, extend_num001000, phi_boundary);
	//assembly_extend_phi_auxiliary_vector3D(phi_auxiliary_vector000010, extend_num000010, phi_boundary);
	//assembly_extend_phi_auxiliary_vector3D(phi_auxiliary_vector111111, extend_num111111, phi_boundary);
	//assembly_extend_phi_auxiliary_vector3D(phi_auxiliary_vector333333, extend_num333333, phi_boundary);
}

//***************************************************************************************************************************************************

//实际问题的计算======================================================================================================================================

//***************************************************************************************************************************************************

//计算=============================================================================================

//更新速度与压力---------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::solve_fluid_structure_interaction_system3D()
{
	//固体-------------------

	//输入速度
	get_velocity_at_p3D();

	//更新phi
	update_phi3D();
	modified_phi_based_on_intermediate_step();

	//内部区域
	assembly_extract_interior_domain_vector3D(extract_interior_domain_vector);

	//速度导数
	get_d_velocity_at_p3D();

	//更新FF
	update_FF3D();

	//更新FF_inv
	update_FF_inv3D();

	//流固耦合---------------

	//NS方程右端项
	get_fluid_solid_NSeq_right_term3D();
	VecAXPY(u_NS_right_term, 1.0, modified_u_NSeq_vector);
	MatMult(take_out_x_L_boundary_matrix3D, u_NS_right_term, u_NS_right_term_take_out_x_L);
	MatMult(take_out_y_L_boundary_matrix3D, v_NS_right_term, v_NS_right_term_take_out_y_L);
	MatMult(take_out_z_L_boundary_matrix3D, w_NS_right_term, w_NS_right_term_take_out_z_L);

	//求解NS方程
	KSPSolve(KSP_NSeq_u, u_NS_right_term_take_out_x_L, velocity_u_star_take_out_x_L);
	KSPSolve(KSP_NSeq_v, v_NS_right_term_take_out_y_L, velocity_v_star_take_out_y_L);
	KSPSolve(KSP_NSeq_w, w_NS_right_term_take_out_z_L, velocity_w_star_take_out_z_L);

	MatMult(x_solution_out_x_L_extend_to_ori_size_matrix3D, velocity_u_star_take_out_x_L, velocity_u_star);
	MatMult(y_solution_out_y_L_extend_to_ori_size_matrix3D, velocity_v_star_take_out_y_L, velocity_v_star);
	MatMult(z_solution_out_z_L_extend_to_ori_size_matrix3D, velocity_w_star_take_out_z_L, velocity_w_star);

	//修正项方程
	get_modified_equation_right_term3D();
	KSPSolve(KSP_NSeq_p, modified_equation_right_term, modified_term_at_p);

	//更新速度和压力
	modified_velocity3D();
	modified_pressure3D();
}

//本构模型计算应力--------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_NS_elastic_force3D()
{
	//赋值-------------------

	//tensor
	VecWAXPY(tensor00_transition, -1.0, FF_inv00, FF00);	VecWAXPY(tensor01_transition, -1.0, FF_inv10, FF01);	VecWAXPY(tensor02_transition, -1.0, FF_inv20, FF02);
	VecWAXPY(tensor10_transition, -1.0, FF_inv01, FF10);	VecWAXPY(tensor11_transition, -1.0, FF_inv11, FF11);	VecWAXPY(tensor12_transition, -1.0, FF_inv21, FF12);
	VecWAXPY(tensor20_transition, -1.0, FF_inv02, FF20);	VecWAXPY(tensor21_transition, -1.0, FF_inv12, FF21);	VecWAXPY(tensor22_transition, -1.0, FF_inv22, FF22);

	//tensor00
	VecPointwiseMult(mult_part01, tensor00_transition, FF00);
	VecPointwiseMult(mult_part02, tensor01_transition, FF01);
	VecPointwiseMult(tensor00, tensor02_transition, FF02);
	VecAXPBYPCZ(tensor00, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//tensor01
	VecPointwiseMult(mult_part01, tensor00_transition, FF10);
	VecPointwiseMult(mult_part02, tensor01_transition, FF11);
	VecPointwiseMult(tensor01, tensor02_transition, FF12);
	VecAXPBYPCZ(tensor01, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//tensor02
	VecPointwiseMult(mult_part01, tensor00_transition, FF20);
	VecPointwiseMult(mult_part02, tensor01_transition, FF21);
	VecPointwiseMult(tensor02, tensor02_transition, FF22);
	VecAXPBYPCZ(tensor02, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//tensor10
	VecPointwiseMult(mult_part01, tensor10_transition, FF00);
	VecPointwiseMult(mult_part02, tensor11_transition, FF01);
	VecPointwiseMult(tensor10, tensor12_transition, FF02);
	VecAXPBYPCZ(tensor10, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//tensor11
	VecPointwiseMult(mult_part01, tensor10_transition, FF10);
	VecPointwiseMult(mult_part02, tensor11_transition, FF11);
	VecPointwiseMult(tensor11, tensor12_transition, FF12);
	VecAXPBYPCZ(tensor11, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//tensor12
	VecPointwiseMult(mult_part01, tensor10_transition, FF20);
	VecPointwiseMult(mult_part02, tensor11_transition, FF21);
	VecPointwiseMult(tensor12, tensor12_transition, FF22);
	VecAXPBYPCZ(tensor12, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//tensor20
	VecPointwiseMult(mult_part01, tensor20_transition, FF00);
	VecPointwiseMult(mult_part02, tensor21_transition, FF01);
	VecPointwiseMult(tensor20, tensor22_transition, FF02);
	VecAXPBYPCZ(tensor20, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//tensor21
	VecPointwiseMult(mult_part01, tensor20_transition, FF10);
	VecPointwiseMult(mult_part02, tensor21_transition, FF11);
	VecPointwiseMult(tensor21, tensor22_transition, FF12);
	VecAXPBYPCZ(tensor21, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//tensor22
	VecPointwiseMult(mult_part01, tensor20_transition, FF20);
	VecPointwiseMult(mult_part02, tensor21_transition, FF21);
	VecPointwiseMult(tensor22, tensor22_transition, FF22);
	VecAXPBYPCZ(tensor22, mu_s, mu_s, mu_s, mult_part01, mult_part02);

	//体力
	MatMult(extend_phi_matrix111111, tensor00, tensor00_111111);
	MatMult(extend_phi_matrix111111, tensor01, tensor01_111111);
	MatMult(extend_phi_matrix111111, tensor02, tensor02_111111);
	MatMult(extend_phi_matrix111111, tensor10, tensor10_111111);
	MatMult(extend_phi_matrix111111, tensor11, tensor11_111111);
	MatMult(extend_phi_matrix111111, tensor12, tensor12_111111);
	MatMult(extend_phi_matrix111111, tensor20, tensor20_111111);
	MatMult(extend_phi_matrix111111, tensor21, tensor21_111111);
	MatMult(extend_phi_matrix111111, tensor22, tensor22_111111);

	MatMult(dx_matrix111111_to_000000, tensor00_111111, d_tensor00_x);
	MatMult(dx_matrix111111_to_000000, tensor01_111111, d_tensor01_x);
	MatMult(dx_matrix111111_to_000000, tensor02_111111, d_tensor02_x);

	MatMult(dy_matrix111111_to_000000, tensor10_111111, d_tensor10_y);
	MatMult(dy_matrix111111_to_000000, tensor11_111111, d_tensor11_y);
	MatMult(dy_matrix111111_to_000000, tensor12_111111, d_tensor12_y);

	MatMult(dz_matrix111111_to_000000, tensor20_111111, d_tensor20_z);
	MatMult(dz_matrix111111_to_000000, tensor21_111111, d_tensor21_z);
	MatMult(dz_matrix111111_to_000000, tensor22_111111, d_tensor22_z);

	VecWAXPY(tensor_body_x_at_p, 1.0, d_tensor00_x, d_tensor10_y);
	VecAXPY(tensor_body_x_at_p, 1.0, d_tensor20_z);

	VecWAXPY(tensor_body_y_at_p, 1.0, d_tensor01_x, d_tensor11_y);
	VecAXPY(tensor_body_y_at_p, 1.0, d_tensor21_z);

	VecWAXPY(tensor_body_z_at_p, 1.0, d_tensor02_x, d_tensor12_y);
	VecAXPY(tensor_body_z_at_p, 1.0, d_tensor22_z);
	/*
	//面力
	solve_phi_normal();

	VecPointwiseMult(mult_part01, tensor00, phi_normal_x);
	VecPointwiseMult(mult_part02, tensor01, phi_normal_y);
	VecPointwiseMult(tensor_normal_x_at_p, tensor02, phi_normal_z);
	VecAXPBYPCZ(tensor_normal_x_at_p, 1.0, 1.0, 1.0, mult_part01, mult_part02);

	VecPointwiseMult(mult_part01, tensor10, phi_normal_x);
	VecPointwiseMult(mult_part02, tensor11, phi_normal_y);
	VecPointwiseMult(tensor_normal_y_at_p, tensor12, phi_normal_z);
	VecAXPBYPCZ(tensor_normal_y_at_p, 1.0, 1.0, 1.0, mult_part01, mult_part02);

	VecPointwiseMult(mult_part01, tensor20, phi_normal_x);
	VecPointwiseMult(mult_part02, tensor21, phi_normal_y);
	VecPointwiseMult(tensor_normal_z_at_p, tensor22, phi_normal_z);
	VecAXPBYPCZ(tensor_normal_z_at_p, 1.0, 1.0, 1.0, mult_part01, mult_part02);
	
	//
	VecWAXPY(F_x_at_p, 1.0, tensor_normal_x_at_p, tensor_body_x_at_p);
	VecWAXPY(F_y_at_p, 1.0, tensor_normal_y_at_p, tensor_body_y_at_p);
	VecWAXPY(F_z_at_p, 1.0, tensor_normal_z_at_p, tensor_body_z_at_p);
	*/

	VecCopy(tensor_body_x_at_p, F_x_at_p);
	VecCopy(tensor_body_y_at_p, F_y_at_p);
	VecCopy(tensor_body_z_at_p, F_z_at_p);

	MatMult(extend_phi_matrix100000, F_x_at_p, F_x_at_p100000);
	MatMult(extend_phi_matrix001000, F_y_at_p, F_y_at_p001000);
	MatMult(extend_phi_matrix000010, F_z_at_p, F_z_at_p000010);
	//VecAXPY(F_x_at_p100000, 1.0, phi_auxiliary_vector100000);
	//VecAXPY(F_y_at_p001000, 1.0, phi_auxiliary_vector001000);
	//VecAXPY(F_z_at_p000010, 1.0, phi_auxiliary_vector000010);

	MatMult(average_matrix100000_or_010000_to_000000, F_x_at_p100000, F_x_at_u);
	MatMult(average_matrix001000_or_000100_to_000000, F_y_at_p001000, F_y_at_v);
	MatMult(average_matrix000010_or_000001_to_000000, F_z_at_p000010, F_z_at_w);
}

//初始值-----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_initial_condition()
{
	//速度
	VecSet(velocity_u, 0.0);
	VecSet(velocity_v, 0.0);
	VecSet(velocity_w, 0.0);

	//压力
	VecSet(pressure, 0.0);

	//phi
	VecSet(phi, 0.0);

	VecGetOwnershipRange(phi, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++) {
		PetscInt k = Ii / (m_num*n_num);
		PetscInt Ii02 = Ii - k*m_num*n_num;
		PetscInt i = Ii02 / m_num;
		PetscInt j = Ii02 - i*m_num;
		
		double node_x_temp = x_L + hx / 2.0 + j*hx;
		double node_y_temp = y_L + hy / 2.0 + i*hy;
		double node_z_temp = z_L + hz / 2.0 + k*hz;

		double phi_temp = sqrt(pow(node_x_temp - p_center[0], 2) + pow(node_y_temp - p_center[1], 2) + pow(node_z_temp - p_center[2], 2)) - r;
		
		if (phi_temp > epsilon){
			phi_temp = 0.0;
		}
		else{
			if (phi_temp < -epsilon){
				phi_temp = 1.0;
			}
			else{
				phi_temp = 0.5 - phi_temp / 2.0 / epsilon - 1.0 / 2.0 / pi*sin(pi*phi_temp / epsilon);
			}
		}
		VecSetValues(phi, 1, &Ii, &phi_temp, INSERT_VALUES);
	}
	VecAssemblyBegin(phi);
	VecAssemblyEnd(phi);

	//FF
	VecCopy(phi, FF00);
	VecSet(FF01, 0.0);
	VecSet(FF02, 0.0);
	VecSet(FF10, 0.0);
	VecCopy(phi, FF11);
	VecSet(FF12, 0.0);
	VecSet(FF20, 0.0);
	VecSet(FF21, 0.0);
	VecCopy(phi, FF22);

	//FF_inv
	VecCopy(phi, FF_inv00);
	VecSet(FF_inv01, 0.0);
	VecSet(FF_inv02, 0.0);
	VecSet(FF_inv10, 0.0);
	VecCopy(phi, FF_inv11);
	VecSet(FF_inv12, 0.0);
	VecSet(FF_inv20, 0.0);
	VecSet(FF_inv21, 0.0);
	VecCopy(phi, FF_inv22);
}

//方程组系数========================================================================================

//流体速度系数------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_u3D(Mat NS_coeff_u3D)
{
	int m_num_at_coeff = m_num - 1;
	int n_num_at_coeff = n_num;
	int o_num_at_coeff = o_num;

	MatSetFromOptions(NS_coeff_u3D);
	MatMPIAIJSetPreallocation(NS_coeff_u3D, 7, NULL, 7, NULL);
	MatSeqAIJSetPreallocation(NS_coeff_u3D, 7, NULL);
	MatZeroEntries(NS_coeff_u3D);

	MatGetOwnershipRange(NS_coeff_u3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++) {
		PetscInt k = Ii / (m_num_at_coeff*n_num_at_coeff);
		PetscInt Ii02 = Ii - k*m_num_at_coeff*n_num_at_coeff;
		PetscInt i = Ii02 / m_num_at_coeff;
		PetscInt j = Ii02 - i*m_num_at_coeff;

		PetscInt Ij;
		double insert_value;

		insert_value = rho*1.0 + mu_f * dt / 2.0 * (2.0 / hx / hx + 2.0 / hy / hy + 2.0 / hz / hz);
		MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ii, &insert_value, ADD_VALUES);

		//i
		if (i > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i - 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + 0 * m_num_at_coeff + j;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (i < n_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (n_num_at_coeff - 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		//j
		if (j > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j - 1;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (j < m_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j + 1;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		//k
		if (k > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			Ij = (k - 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			insert_value = -1.0*insert_value;
			Ij = 0 * m_num_at_coeff*n_num_at_coeff + (i + 0) * m_num_at_coeff + j;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (k < o_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			Ij = (k + 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			insert_value = -1.0*insert_value;
			Ij = (o_num_at_coeff - 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_u3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

	}
	MatAssemblyBegin(NS_coeff_u3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(NS_coeff_u3D, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_v3D(Mat NS_coeff_v3D)
{
	int m_num_at_coeff = m_num;
	int n_num_at_coeff = n_num - 1;
	int o_num_at_coeff = o_num;

	MatSetFromOptions(NS_coeff_v3D);
	MatMPIAIJSetPreallocation(NS_coeff_v3D, 7, NULL, 7, NULL);
	MatSeqAIJSetPreallocation(NS_coeff_v3D, 7, NULL);
	MatZeroEntries(NS_coeff_v3D);

	MatGetOwnershipRange(NS_coeff_v3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++) {
		PetscInt k = Ii / (m_num_at_coeff*n_num_at_coeff);
		PetscInt Ii02 = Ii - k*m_num_at_coeff*n_num_at_coeff;
		PetscInt i = Ii02 / m_num_at_coeff;
		PetscInt j = Ii02 - i*m_num_at_coeff;

		PetscInt Ij;
		double insert_value;

		insert_value = rho*1.0 + mu_f * dt / 2.0 * (2.0 / hx / hx + 2.0 / hy / hy + 2.0 / hz / hz);
		MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ii, &insert_value, ADD_VALUES);

		//i
		if (i > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i - 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (i < n_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		//j
		if (j > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j - 1;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + 0;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (j < m_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j + 1;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + m_num_at_coeff - 1;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		//k
		if (k > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			Ij = (k - 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			insert_value = -1.0*insert_value;
			Ij = 0 * m_num_at_coeff*n_num_at_coeff + (i + 0) * m_num_at_coeff + j;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (k < o_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			Ij = (k + 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			insert_value = -1.0*insert_value;
			Ij = (o_num_at_coeff - 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_v3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

	}
	MatAssemblyBegin(NS_coeff_v3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(NS_coeff_v3D, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_w3D(Mat NS_coeff_w3D)
{
	int m_num_at_coeff = m_num;
	int n_num_at_coeff = n_num;
	int o_num_at_coeff = o_num - 1;

	MatSetFromOptions(NS_coeff_w3D);
	MatMPIAIJSetPreallocation(NS_coeff_w3D, 7, NULL, 7, NULL);
	MatSeqAIJSetPreallocation(NS_coeff_w3D, 7, NULL);
	MatZeroEntries(NS_coeff_w3D);

	MatGetOwnershipRange(NS_coeff_w3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++) {
		PetscInt k = Ii / (m_num_at_coeff*n_num_at_coeff);
		PetscInt Ii02 = Ii - k*m_num_at_coeff*n_num_at_coeff;
		PetscInt i = Ii02 / m_num_at_coeff;
		PetscInt j = Ii02 - i*m_num_at_coeff;

		PetscInt Ij;
		double insert_value;

		insert_value = rho*1.0 + mu_f * dt / 2.0 * (2.0 / hx / hx + 2.0 / hy / hy + 2.0 / hz / hz);
		MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ii, &insert_value, ADD_VALUES);

		//i
		if (i > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i - 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + 0 * m_num_at_coeff + j;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (i < n_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hy / hy;
			insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (n_num_at_coeff - 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		//j
		if (j > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j - 1;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + 0;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (j < m_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j + 1;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = -1.0 * mu_f * dt / 2.0 / hx / hx;
			insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + m_num_at_coeff - 1;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		//k
		if (k > 0){
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			Ij = (k - 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (k < o_num_at_coeff - 1){
			insert_value = -1.0 * mu_f * dt / 2.0 / hz / hz;
			Ij = (k + 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_w3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

	}
	MatAssemblyBegin(NS_coeff_w3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(NS_coeff_w3D, MAT_FINAL_ASSEMBLY);
}

//流体压力系数------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_p3D(Mat NS_coeff_p3D)
{
	int m_num_at_coeff = m_num;
	int n_num_at_coeff = n_num;
	int o_num_at_coeff = o_num;

	MatSetFromOptions(NS_coeff_p3D);
	MatMPIAIJSetPreallocation(NS_coeff_p3D, 7, NULL, 7, NULL);
	MatSeqAIJSetPreallocation(NS_coeff_p3D, 7, NULL);
	MatZeroEntries(NS_coeff_p3D);

	MatGetOwnershipRange(NS_coeff_p3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++) {
		PetscInt k = Ii / (m_num_at_coeff*n_num_at_coeff);
		PetscInt Ii02 = Ii - k*m_num_at_coeff*n_num_at_coeff;
		PetscInt i = Ii02 / m_num_at_coeff;
		PetscInt j = Ii02 - i*m_num_at_coeff;

		PetscInt Ij;
		double insert_value;

		insert_value = -dt * (2.0 / hx / hx + 2.0 / hy / hy + 2.0 / hz / hz);
		MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ii, &insert_value, ADD_VALUES);

		//i
		if (i > 0){
			insert_value = dt*1.0 / hy / hy;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i - 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = dt*1.0 / hy / hy;
			//insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + 0 * m_num_at_coeff + j;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (i < n_num_at_coeff - 1){
			insert_value = dt*1.0 / hy / hy;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = dt*1.0 / hy / hy;
			//insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (n_num_at_coeff - 1)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		//j
		if (j > 0){
			insert_value = dt*1.0 / hx / hx;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j - 1;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = dt*1.0 / hx / hx;
			//insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + 0;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (j < m_num_at_coeff - 1){
			insert_value = dt*1.0 / hx / hx;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j + 1;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = dt*1.0 / hx / hx;
			//insert_value = -1.0*insert_value;
			Ij = (k + 0)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + m_num_at_coeff - 1;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		//k
		if (k > 0){
			insert_value = dt*1.0 / hz / hz;
			Ij = (k - 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = dt*1.0 / hz / hz;
			//insert_value = -1.0*insert_value;
			Ij = 0 * m_num_at_coeff*n_num_at_coeff + (i + 0) * m_num_at_coeff + j;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

		if (k < o_num_at_coeff - 1){
			insert_value = dt*1.0 / hz / hz;
			Ij = (k + 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}
		else{
			insert_value = dt*1.0 / hz / hz;
			//insert_value = -1.0*insert_value;
			Ij = (o_num_at_coeff - 1)*m_num_at_coeff*n_num_at_coeff + (i + 0)*m_num_at_coeff + j;
			MatSetValues(NS_coeff_p3D, 1, &Ii, 1, &Ij, &insert_value, ADD_VALUES);
		}

	}
	MatAssemblyBegin(NS_coeff_p3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(NS_coeff_p3D, MAT_FINAL_ASSEMBLY);
}

//处理实际问题的特殊向量与矩阵=======================================================================

//处理NS方程右端项--------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_modified_u_NSeq_vector_for_lid_driven_cavity_flow_boundary_condition3D()
{
	//组装
	VecGetOwnershipRange(modified_u_NSeq_vector, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++) {
		PetscInt k = Ii / (m_num*n_num);
		PetscInt Ii02 = Ii - k * (m_num*n_num);
		PetscInt i = Ii02 / m_num;
		PetscInt j = Ii02 - i*m_num;

		if (k == o_num - 1){
			double insert_value = mu_f*dt / 2.0*(2.0 / hz / hz);
			VecSetValues(modified_u_NSeq_vector, 1, &Ii, &insert_value, INSERT_VALUES);
		}
	}
	VecAssemblyBegin(modified_u_NSeq_vector);
	VecAssemblyEnd(modified_u_NSeq_vector);
}

//去掉边界数据------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_take_out_x_L_boundary_matrix3D()
{
	int m_num_take_out_boundary = m_num - 1;
	int n_num_take_out_boundary = n_num;
	int o_num_take_out_boundary = o_num;

	//组装
	MatMPIAIJSetPreallocation(take_out_x_L_boundary_matrix3D, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(take_out_x_L_boundary_matrix3D, 1, NULL);
	MatZeroEntries(take_out_x_L_boundary_matrix3D);

	MatGetOwnershipRange(take_out_x_L_boundary_matrix3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (m_num_take_out_boundary*n_num_take_out_boundary);
		PetscInt Ii02 = Ii - k * m_num_take_out_boundary*n_num_take_out_boundary;
		PetscInt i = Ii02 / m_num_take_out_boundary;
		PetscInt j = Ii02 - i*m_num_take_out_boundary;

		PetscInt Ij;
		double insert_value = 1.0;

		Ij = k*m_num*n_num + i*m_num + (j + 1);

		MatSetValues(take_out_x_L_boundary_matrix3D, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
	}
	MatAssemblyBegin(take_out_x_L_boundary_matrix3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(take_out_x_L_boundary_matrix3D, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_take_out_y_L_boundary_matrix3D()
{
	int m_num_take_out_boundary = m_num;
	int n_num_take_out_boundary = n_num - 1;
	int o_num_take_out_boundary = o_num;

	//组装
	MatMPIAIJSetPreallocation(take_out_y_L_boundary_matrix3D, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(take_out_y_L_boundary_matrix3D, 1, NULL);
	MatZeroEntries(take_out_y_L_boundary_matrix3D);

	MatGetOwnershipRange(take_out_y_L_boundary_matrix3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (m_num_take_out_boundary*n_num_take_out_boundary);
		PetscInt Ii02 = Ii - k * m_num_take_out_boundary*n_num_take_out_boundary;
		PetscInt i = Ii02 / m_num_take_out_boundary;
		PetscInt j = Ii02 - i*m_num_take_out_boundary;

		PetscInt Ij;
		double insert_value = 1.0;

		Ij = k*m_num*n_num + (i + 1)*m_num + j;

		MatSetValues(take_out_y_L_boundary_matrix3D, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
	}
	MatAssemblyBegin(take_out_y_L_boundary_matrix3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(take_out_y_L_boundary_matrix3D, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_take_out_z_L_boundary_matrix3D()
{
	int m_num_take_out_boundary = m_num;
	int n_num_take_out_boundary = n_num;
	int o_num_take_out_boundary = o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(take_out_z_L_boundary_matrix3D, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(take_out_z_L_boundary_matrix3D, 1, NULL);
	MatZeroEntries(take_out_z_L_boundary_matrix3D);

	MatGetOwnershipRange(take_out_z_L_boundary_matrix3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (m_num_take_out_boundary*n_num_take_out_boundary);
		PetscInt Ii02 = Ii - k * m_num_take_out_boundary*n_num_take_out_boundary;
		PetscInt i = Ii02 / m_num_take_out_boundary;
		PetscInt j = Ii02 - i*m_num_take_out_boundary;

		PetscInt Ij;
		double insert_value = 1.0;

		Ij = (k + 1)*m_num*n_num + i*m_num + j;

		MatSetValues(take_out_z_L_boundary_matrix3D, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
	}
	MatAssemblyBegin(take_out_z_L_boundary_matrix3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(take_out_z_L_boundary_matrix3D, MAT_FINAL_ASSEMBLY);
}

//添加边界数据------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_x_solution_out_x_L_extend_to_ori_size_matrix3D()
{
	//组装
	MatMPIAIJSetPreallocation(x_solution_out_x_L_extend_to_ori_size_matrix3D, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(x_solution_out_x_L_extend_to_ori_size_matrix3D, 1, NULL);
	MatZeroEntries(x_solution_out_x_L_extend_to_ori_size_matrix3D);

	MatGetOwnershipRange(x_solution_out_x_L_extend_to_ori_size_matrix3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (m_num*n_num);
		PetscInt Ii02 = Ii - k * m_num*n_num;
		PetscInt i = Ii02 / m_num;
		PetscInt j = Ii02 - i * m_num;

		PetscInt Ij;

		double insert_value = 1.0;
		if (j > 0){
			Ij = k*(m_num - 1)*n_num + i*(m_num - 1) + (j - 1);
			MatSetValues(x_solution_out_x_L_extend_to_ori_size_matrix3D, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}
	MatAssemblyBegin(x_solution_out_x_L_extend_to_ori_size_matrix3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(x_solution_out_x_L_extend_to_ori_size_matrix3D, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_y_solution_out_y_L_extend_to_ori_size_matrix3D()
{
	//组装
	MatMPIAIJSetPreallocation(y_solution_out_y_L_extend_to_ori_size_matrix3D, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(y_solution_out_y_L_extend_to_ori_size_matrix3D, 1, NULL);
	MatZeroEntries(y_solution_out_y_L_extend_to_ori_size_matrix3D);

	MatGetOwnershipRange(y_solution_out_y_L_extend_to_ori_size_matrix3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (m_num*n_num);
		PetscInt Ii02 = Ii - k * m_num*n_num;
		PetscInt i = Ii02 / m_num;
		PetscInt j = Ii02 - i * m_num;

		PetscInt Ij;

		double insert_value = 1.0;
		if (i > 0){
			Ij = k*m_num*(n_num - 1) + (i - 1)*m_num + j;
			MatSetValues(y_solution_out_y_L_extend_to_ori_size_matrix3D, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}
	MatAssemblyBegin(y_solution_out_y_L_extend_to_ori_size_matrix3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(y_solution_out_y_L_extend_to_ori_size_matrix3D, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_z_solution_out_z_L_extend_to_ori_size_matrix3D()
{
	//组装
	MatMPIAIJSetPreallocation(z_solution_out_z_L_extend_to_ori_size_matrix3D, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(z_solution_out_z_L_extend_to_ori_size_matrix3D, 1, NULL);
	MatZeroEntries(z_solution_out_z_L_extend_to_ori_size_matrix3D);

	MatGetOwnershipRange(z_solution_out_z_L_extend_to_ori_size_matrix3D, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (m_num*n_num);
		PetscInt Ii02 = Ii - k * m_num*n_num;
		PetscInt i = Ii02 / m_num;
		PetscInt j = Ii02 - i * m_num;

		PetscInt Ij;

		double insert_value = 1.0;
		if (k > 0){
			Ij = (k - 1)*m_num*n_num + i*m_num + j;
			MatSetValues(z_solution_out_z_L_extend_to_ori_size_matrix3D, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}
	MatAssemblyBegin(z_solution_out_z_L_extend_to_ori_size_matrix3D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(z_solution_out_z_L_extend_to_ori_size_matrix3D, MAT_FINAL_ASSEMBLY);
}

//提取边界区域的矩阵------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_extract_boundary_domain_vector3D(Vec extract_boundary_domain_vector)
{
	VecSet(extract_boundary_domain_vector, 0.0);

	double *phi_temp;
	VecGetArray(phi, &phi_temp);

	VecGetOwnershipRange(extract_boundary_domain_vector, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		if (phi_temp[Ii - istart] >= 0.01 && phi_temp[Ii - istart] <= 0.99){
			double insert_value = 1.0;
			VecSetValues(extract_boundary_domain_vector, 1, &Ii, &insert_value, INSERT_VALUES);
		}
	}
	VecAssemblyBegin(extract_boundary_domain_vector);
	VecAssemblyEnd(extract_boundary_domain_vector);
	VecRestoreArray(phi, &phi_temp);

	delete[] phi_temp;
	phi_temp = NULL;
}

//提取内部区域的矩阵------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_extract_interior_domain_vector3D(Vec extract_interior_domain_vector)
{
	VecSet(extract_interior_domain_vector, 0.0);

	double *phi_temp;
	VecGetArray(phi, &phi_temp);

	VecGetOwnershipRange(extract_interior_domain_vector, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		if (phi_temp[Ii - istart] > 0.01){
			double insert_value = 1.0;
			VecSetValues(extract_interior_domain_vector, 1, &Ii, &insert_value, INSERT_VALUES);
		}
	}
	VecAssemblyBegin(extract_interior_domain_vector);
	VecAssemblyEnd(extract_interior_domain_vector);
	VecRestoreArray(phi, &phi_temp);

	delete[] phi_temp;
	phi_temp = NULL;
}

//***************************************************************************************************************************************************

//求解N-S方程=========================================================================================================================================

//***************************************************************************************************************************************************

//修正速度及压力项==================================================================================

//修正速度----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::modified_velocity3D()
{
	//赋值-------------------

	//拓展修正项
	MatMult(extend_pressure_matrix100000, modified_term_at_p, modified_term_at_p100000);
	MatMult(extend_pressure_matrix001000, modified_term_at_p, modified_term_at_p001000);
	MatMult(extend_pressure_matrix000010, modified_term_at_p, modified_term_at_p000010);
	//VecAXPY(modified_term_at_p100000, 1.0, pressure_auxiliary_vector100000);
	//VecAXPY(modified_term_at_p001000, 1.0, pressure_auxiliary_vector001000);
	//VecAXPY(modified_term_at_p000010, 1.0, pressure_auxiliary_vector000010);

	//修正项导数
	MatMult(dx_matrix100000_or_010000_to_000000, modified_term_at_p100000, d_modified_term_x);
	MatMult(dy_matrix001000_or_000100_to_000000, modified_term_at_p001000, d_modified_term_y);
	MatMult(dz_matrix000010_or_000001_to_000000, modified_term_at_p000010, d_modified_term_z);

	//修正速度
	VecWAXPY(velocity_u, -1.0*dt / rho, d_modified_term_x, velocity_u_star);
	VecWAXPY(velocity_v, -1.0*dt / rho, d_modified_term_y, velocity_v_star);
	VecWAXPY(velocity_w, -1.0*dt / rho, d_modified_term_z, velocity_w_star);
}

//修正压力----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::modified_pressure3D()
{
	//赋值-------------------

	//拓展修正项
	MatMult(extend_pressure_matrix111111, modified_term_at_p, modified_term_at_p111111);
	//VecAXPY(modified_term_at_p111111, 1.0, pressure_auxiliary_vector111111);

	//计算拉普拉斯项
	MatMult(ddx_matrix111111_to_000000, modified_term_at_p111111, dd_modified_term_xx);
	MatMult(ddy_matrix111111_to_000000, modified_term_at_p111111, dd_modified_term_yy);
	MatMult(ddz_matrix111111_to_000000, modified_term_at_p111111, lap_modified_term);

	VecAXPBYPCZ(lap_modified_term, 1.0, 1.0, 1.0, dd_modified_term_xx, dd_modified_term_yy);

	//修正压力
	VecCopy(modified_term_at_p, pressure);
	VecAXPY(pressure, -1.0*dt*mu_f / rho / 2.0, lap_modified_term);
}

//求解NS方程所涉及的右端项===========================================================================

//修正方程的右端项--------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_modified_equation_right_term3D()
{
	//赋值-------------------

	//右端项初始化
	VecSet(modified_equation_right_term, 0.0);

	//拓展过渡速度
	MatMult(extend_u_matrix010000, velocity_u_star, velocity_u_star_at_u010000);
	MatMult(extend_v_matrix000100, velocity_v_star, velocity_v_star_at_v000100);
	MatMult(extend_w_matrix000001, velocity_w_star, velocity_w_star_at_w000001);
	VecAXPY(velocity_u_star_at_u010000, 1.0, u_auxiliary_vector010000);
	//VecAXPY(velocity_v_star_at_v001000, 1.0, v_auxiliary_vector000100);
	//VecAXPY(velocity_w_star_at_w000001, 1.0, w_auxiliary_vector000001);

	//过渡速度导数
	MatMult(dx_matrix100000_or_010000_to_000000, velocity_u_star_at_u010000, d_u_star_x);
	MatMult(dy_matrix001000_or_000100_to_000000, velocity_v_star_at_v000100, d_v_star_y);
	MatMult(dz_matrix000010_or_000001_to_000000, velocity_w_star_at_w000001, modified_equation_right_term);

	//右端项
	VecAXPBYPCZ(modified_equation_right_term, rho, rho, rho, d_u_star_x, d_v_star_y);
}

//NS右端项---------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_fluid_solid_NSeq_right_term3D()
{
	//赋值-------------------

	//计算固体应力
	get_NS_elastic_force3D();

	//压力
	//get_d_pressure_at_stagger3D(PETSC_COMM_WORLD, d_pressure_x_at_u, d_pressure_y_at_v, d_pressure_z_at_w);

	//对流项
	get_NS_diffusion3D();

	//扩散项
	if (judge_num == 1){
		get_NS_convection3D();
		VecCopy(u_NS_convection, u_NS_convection_ori);
		VecCopy(v_NS_convection, v_NS_convection_ori);
		VecCopy(w_NS_right_term, w_NS_convection_ori);
	}
	else{
		get_NS_convection3D();
		VecCopy(u_NS_convection, u_NS_right_term);
		VecCopy(v_NS_convection, v_NS_right_term);
		VecCopy(w_NS_convection, w_NS_right_term);
		VecScale(u_NS_convection, 3.0 / 2.0);
		VecScale(v_NS_convection, 3.0 / 2.0);
		VecScale(w_NS_convection, 3.0 / 2.0);
		VecAXPY(u_NS_convection, -1.0 / 2.0, u_NS_convection_ori);
		VecAXPY(v_NS_convection, -1.0 / 2.0, v_NS_convection_ori);
		VecAXPY(w_NS_convection, -1.0 / 2.0, w_NS_convection_ori);
		VecCopy(u_NS_right_term, u_NS_convection_ori);
		VecCopy(v_NS_right_term, v_NS_convection_ori);
		VecCopy(w_NS_right_term, w_NS_convection_ori);
	}

	//方程右端项初始化
	VecSet(u_NS_right_term, 0.0);
	VecSet(v_NS_right_term, 0.0);
	VecSet(w_NS_right_term, 0.0);

	//获得N-S方程右端项
	VecAXPY(u_NS_right_term, -1.0*rho*dt, u_NS_convection);
	VecAXPY(v_NS_right_term, -1.0*rho*dt, v_NS_convection);
	VecAXPY(w_NS_right_term, -1.0*rho*dt, w_NS_convection);

	VecAXPY(u_NS_right_term, dt*mu_f / 2.0, u_NS_diffusion);
	VecAXPY(v_NS_right_term, dt*mu_f / 2.0, v_NS_diffusion);
	VecAXPY(w_NS_right_term, dt*mu_f / 2.0, w_NS_diffusion);

	VecAXPY(u_NS_right_term, dt, F_x_at_u);
	VecAXPY(v_NS_right_term, dt, F_y_at_v);
	VecAXPY(w_NS_right_term, dt, F_z_at_w);

	//VecAXPY(u_NS_right_term, -1.0*dt, d_pressure_x_at_u);
	//VecAXPY(v_NS_right_term, -1.0*dt, d_pressure_y_at_v);
	//VecAXPY(w_NS_right_term, -1.0*dt, d_pressure_z_at_w);

	VecAXPY(u_NS_right_term, rho, velocity_u);
	VecAXPY(v_NS_right_term, rho, velocity_v);
	VecAXPY(w_NS_right_term, rho, velocity_w);
}

//压力的梯度--------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_d_pressure_at_stagger3D()
{
	//赋值-------------------

	//拓展压力
	MatMult(extend_pressure_matrix100000, pressure, pressure_at_p100000);
	MatMult(extend_pressure_matrix001000, pressure, pressure_at_p001000);
	MatMult(extend_pressure_matrix000010, pressure, pressure_at_p000010);
	//VecAXPY(pressure_at_p100000, 1.0, pressure_auxiliary_vector100000);
	//VecAXPY(pressure_at_p001000, 1.0, pressure_auxiliary_vector001000);
	//VecAXPY(pressure_at_p000010, 1.0, pressure_auxiliary_vector000010);

	//计算压力的梯度
	MatMult(dx_matrix100000_or_010000_to_000000, pressure_at_p100000, d_pressure_x_at_u);
	MatMult(dy_matrix001000_or_000100_to_000000, pressure_at_p001000, d_pressure_y_at_v);
	MatMult(dz_matrix000010_or_000001_to_000000, pressure_at_p000010, d_pressure_z_at_w);
}

//扩散项-----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_NS_diffusion3D()
{
	//赋值-------------------

	//拓展速度
	MatMult(extend_u_matrix111111, velocity_u, velocity_u_at_u111111);
	MatMult(extend_v_matrix111111, velocity_v, velocity_v_at_v111111);
	MatMult(extend_w_matrix111111, velocity_w, velocity_w_at_w111111);
	VecAXPY(velocity_u_at_u111111, 1.0, u_auxiliary_vector111111);
	//VecAXPY(velocity_v_at_v111111, 1.0, v_auxiliary_vector111111);
	//VecAXPY(velocity_w_at_w111111, 1.0, w_auxiliary_vector111111);

	//计算扩散项
	MatMult(ddx_matrix111111_to_000000, velocity_u_at_u111111, dd_u_xx);
	MatMult(ddy_matrix111111_to_000000, velocity_u_at_u111111, dd_u_yy);
	MatMult(ddz_matrix111111_to_000000, velocity_u_at_u111111, u_NS_diffusion);
	MatMult(ddx_matrix111111_to_000000, velocity_v_at_v111111, dd_v_xx);
	MatMult(ddy_matrix111111_to_000000, velocity_v_at_v111111, dd_v_yy);
	MatMult(ddz_matrix111111_to_000000, velocity_v_at_v111111, v_NS_diffusion);
	MatMult(ddx_matrix111111_to_000000, velocity_w_at_w111111, dd_w_xx);
	MatMult(ddy_matrix111111_to_000000, velocity_w_at_w111111, dd_w_yy);
	MatMult(ddz_matrix111111_to_000000, velocity_w_at_w111111, w_NS_diffusion);

	VecAXPBYPCZ(u_NS_diffusion, 1.0, 1.0, 1.0, dd_u_xx, dd_u_yy);
	VecAXPBYPCZ(v_NS_diffusion, 1.0, 1.0, 1.0, dd_v_xx, dd_v_yy);
	VecAXPBYPCZ(w_NS_diffusion, 1.0, 1.0, 1.0, dd_w_xx, dd_w_yy);
}

//对流项-----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_NS_convection3D()
{
	//赋值-------------------

	//拓展项
	MatMult(extend_u_matrix333333, velocity_u, velocity_u_at_u333333);
	MatMult(extend_u_matrix344333, velocity_u, velocity_u_at_u344333);
	MatMult(extend_u_matrix343343, velocity_u, velocity_u_at_u343343);
	VecAXPY(velocity_u_at_u333333, 1.0, u_auxiliary_vector333333);
	VecAXPY(velocity_u_at_u344333, 1.0, u_auxiliary_vector344333);
	VecAXPY(velocity_u_at_u343343, 1.0, u_auxiliary_vector343343);

	MatMult(extend_v_matrix333333, velocity_v, velocity_v_at_v333333);
	MatMult(extend_v_matrix433433, velocity_v, velocity_v_at_v433433);
	MatMult(extend_v_matrix333443, velocity_v, velocity_v_at_v333443);
	//VecAXPY(velocity_v_at_v333333, 1.0, v_auxiliary_vector333333);
	//VecAXPY(velocity_v_at_v433433, 1.0, v_auxiliary_vector433433);
	//VecAXPY(velocity_v_at_v333443, 1.0, v_auxiliary_vector333443);

	MatMult(extend_w_matrix333333, velocity_w, velocity_w_at_w333333);
	MatMult(extend_w_matrix433334, velocity_w, velocity_w_at_w433334);
	MatMult(extend_w_matrix334334, velocity_w, velocity_w_at_w334334);
	//VecAXPY(velocity_w_at_w333333, 1.0, w_auxiliary_vector333333);
	//VecAXPY(velocity_w_at_w433334, 1.0, w_auxiliary_vector433334);
	//VecAXPY(velocity_w_at_w334334, 1.0, w_auxiliary_vector334334);

	//速度平均值
	MatMult(average_matrix344333_or_433433_to_333333, velocity_u_at_u344333, velocity_u_at_v333333);
	MatMult(average_matrix343343_or_433334_to_333333, velocity_u_at_u343343, velocity_u_at_w333333);

	MatMult(average_matrix344333_or_433433_to_333333, velocity_v_at_v433433, velocity_v_at_u333333);
	MatMult(average_matrix333443_or_334334_to_333333, velocity_v_at_v333443, velocity_v_at_w333333);

	MatMult(average_matrix343343_or_433334_to_333333, velocity_w_at_w433334, velocity_w_at_u333333);
	MatMult(average_matrix333443_or_334334_to_333333, velocity_w_at_w334334, velocity_w_at_v333333);

	//获得PN
	double max_abs_velocity_u;
	double max_abs_velocity_v;
	double max_abs_velocity_w;

	//u
	max_abs_velocity_u = get_max_abs(velocity_u_at_u333333);
	VecPointwiseMult(f_variable333333, velocity_u_at_u333333, velocity_u_at_u333333);
	VecWAXPY(uu_f_P, 2.0*max_abs_velocity_u, velocity_u_at_u333333, f_variable333333);
	VecWAXPY(uu_f_N, -2.0*max_abs_velocity_u, velocity_u_at_u333333, f_variable333333);
	VecScale(uu_f_P, 0.5);
	VecScale(uu_f_N, 0.5);

	max_abs_velocity_v = get_max_abs(velocity_v_at_u333333);
	VecPointwiseMult(f_variable333333, velocity_v_at_u333333, velocity_u_at_u333333);
	VecWAXPY(vu_f_P, 1.0*max_abs_velocity_v, velocity_u_at_u333333, f_variable333333);
	VecWAXPY(vu_f_N, -1.0*max_abs_velocity_v, velocity_u_at_u333333, f_variable333333);
	VecScale(vu_f_P, 0.5);
	VecScale(vu_f_N, 0.5);

	max_abs_velocity_w = get_max_abs(velocity_w_at_u333333);
	VecPointwiseMult(f_variable333333, velocity_w_at_u333333, velocity_u_at_u333333);
	VecWAXPY(wu_f_P, 1.0*max_abs_velocity_w, velocity_u_at_u333333, f_variable333333);
	VecWAXPY(wu_f_N, -1.0*max_abs_velocity_w, velocity_u_at_u333333, f_variable333333);
	VecScale(wu_f_P, 0.5);
	VecScale(wu_f_N, 0.5);

	//v
	max_abs_velocity_u = get_max_abs(velocity_u_at_v333333);
	VecPointwiseMult(f_variable333333, velocity_u_at_v333333, velocity_v_at_v333333);
	VecWAXPY(uv_f_P, 1.0*max_abs_velocity_u, velocity_v_at_v333333, f_variable333333);
	VecWAXPY(uv_f_N, -1.0*max_abs_velocity_u, velocity_v_at_v333333, f_variable333333);
	VecScale(uv_f_P, 0.5);
	VecScale(uv_f_N, 0.5);

	max_abs_velocity_v = get_max_abs(velocity_v_at_v333333);
	VecPointwiseMult(f_variable333333, velocity_v_at_v333333, velocity_v_at_v333333);
	VecWAXPY(vv_f_P, 2.0*max_abs_velocity_v, velocity_v_at_v333333, f_variable333333);
	VecWAXPY(vv_f_N, -2.0*max_abs_velocity_v, velocity_v_at_v333333, f_variable333333);
	VecScale(vv_f_P, 0.5);
	VecScale(vv_f_N, 0.5);

	max_abs_velocity_w = get_max_abs(velocity_w_at_v333333);
	VecPointwiseMult(f_variable333333, velocity_w_at_v333333, velocity_v_at_v333333);
	VecWAXPY(wv_f_P, 1.0*max_abs_velocity_w, velocity_v_at_v333333, f_variable333333);
	VecWAXPY(wv_f_N, -1.0*max_abs_velocity_w, velocity_v_at_v333333, f_variable333333);
	VecScale(wv_f_P, 0.5);
	VecScale(wv_f_N, 0.5);

	//w
	max_abs_velocity_u = get_max_abs(velocity_u_at_w333333);
	VecPointwiseMult(f_variable333333, velocity_u_at_w333333, velocity_w_at_w333333);
	VecWAXPY(uw_f_P, 1.0*max_abs_velocity_u, velocity_w_at_w333333, f_variable333333);
	VecWAXPY(uw_f_N, -1.0*max_abs_velocity_u, velocity_w_at_w333333, f_variable333333);
	VecScale(uw_f_P, 0.5);
	VecScale(uw_f_N, 0.5);

	max_abs_velocity_v = get_max_abs(velocity_v_at_w333333);
	VecPointwiseMult(f_variable333333, velocity_v_at_w333333, velocity_w_at_w333333);
	VecWAXPY(vw_f_P, 1.0*max_abs_velocity_v, velocity_w_at_w333333, f_variable333333);
	VecWAXPY(vw_f_N, -1.0*max_abs_velocity_v, velocity_w_at_w333333, f_variable333333);
	VecScale(vw_f_P, 0.5);
	VecScale(vw_f_N, 0.5);

	max_abs_velocity_w = get_max_abs(velocity_w_at_w333333);
	VecPointwiseMult(f_variable333333, velocity_w_at_w333333, velocity_w_at_w333333);
	VecWAXPY(ww_f_P, 2.0*max_abs_velocity_w, velocity_w_at_w333333, f_variable333333);
	VecWAXPY(ww_f_N, -2.0*max_abs_velocity_w, velocity_w_at_w333333, f_variable333333);
	VecScale(ww_f_P, 0.5);
	VecScale(ww_f_N, 0.5);

	//导数项
	get_WENO_d_variable_x_value(uu_f_P, uu_f_N, d_uu_x);
	get_WENO_d_variable_y_value(vu_f_P, vu_f_N, d_vu_y);
	get_WENO_d_variable_z_value(wu_f_P, wu_f_N, u_NS_convection);

	get_WENO_d_variable_x_value(uv_f_P, uv_f_N, d_uv_x);
	get_WENO_d_variable_y_value(vv_f_P, vv_f_N, d_vv_y);
	get_WENO_d_variable_z_value(wv_f_P, wv_f_N, v_NS_convection);

	get_WENO_d_variable_x_value(uw_f_P, uw_f_N, d_uw_x);
	get_WENO_d_variable_y_value(vw_f_P, vw_f_N, d_vw_y);
	get_WENO_d_variable_z_value(ww_f_P, ww_f_N, w_NS_convection);

	//对流项
	VecAXPBYPCZ(u_NS_convection, 1.0, 1.0, 1.0, d_uu_x, d_vu_y);
	VecAXPBYPCZ(v_NS_convection, 1.0, 1.0, 1.0, d_uv_x, d_vv_y);
	VecAXPBYPCZ(w_NS_convection, 1.0, 1.0, 1.0, d_uw_x, d_vw_y);
}

//***************************************************************************************************************************************************

//更新FF与phi========================================================================================================================================

//***************************************************************************************************************************************************

//转换速度==========================================================================================

//获得中心速度------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_velocity_at_p3D()
{
	MatMult(extend_u_matrix010000, velocity_u, velocity_u_at_u010000);
	MatMult(extend_v_matrix000100, velocity_v, velocity_v_at_v000100);
	MatMult(extend_w_matrix000001, velocity_w, velocity_w_at_w000001);
	VecAXPY(velocity_u_at_u010000, 1.0, u_auxiliary_vector010000);
	//	VecAXPY(velocity_v_at_v000100, 1.0, v_auxiliary_vector000100);
	//	VecAXPY(velocity_w_at_w000001, 1.0, w_auxiliary_vector000001);

	MatMult(average_matrix100000_or_010000_to_000000, velocity_u_at_u010000, velocity_u_at_p);
	MatMult(average_matrix001000_or_000100_to_000000, velocity_v_at_v000100, velocity_v_at_p);
	MatMult(average_matrix000010_or_000001_to_000000, velocity_w_at_w000001, velocity_w_at_p);
}

//速度导数----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_d_velocity_at_p3D()
{
	//拓展向量
	MatMult(nonstagger_extend_u_matrix111111, velocity_u_at_p, velocity_u_at_p111111);
	MatMult(nonstagger_extend_v_matrix111111, velocity_v_at_p, velocity_v_at_p111111);
	MatMult(nonstagger_extend_w_matrix111111, velocity_w_at_p, velocity_w_at_p111111);

	VecAXPY(velocity_u_at_p111111, 1.0, u_auxiliary_vector111111);
	//VecAXPY(velocity_v_at_p111111, 1.0, v_auxiliary_vector111111);
	//VecAXPY(velocity_w_at_p111111, 1.0, w_auxiliary_vector111111);

	//求导
	MatMult(dx_matrix111111_to_000000, velocity_u_at_p111111, d_u_x_at_p);
	MatMult(dy_matrix111111_to_000000, velocity_u_at_p111111, d_u_y_at_p);
	MatMult(dz_matrix111111_to_000000, velocity_u_at_p111111, d_u_z_at_p);
	MatMult(dx_matrix111111_to_000000, velocity_v_at_p111111, d_v_x_at_p);
	MatMult(dy_matrix111111_to_000000, velocity_v_at_p111111, d_v_y_at_p);
	MatMult(dz_matrix111111_to_000000, velocity_v_at_p111111, d_v_z_at_p);
	MatMult(dx_matrix111111_to_000000, velocity_w_at_p111111, d_w_x_at_p);
	MatMult(dy_matrix111111_to_000000, velocity_w_at_p111111, d_w_y_at_p);
	MatMult(dz_matrix111111_to_000000, velocity_w_at_p111111, d_w_z_at_p);
}

//求解level输运方程=================================================================================

//phi输运方程右端项-------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_variable_LSeq_right_term3D(
	Vec variable, Vec variable_LSeq_right_term)
{
	//防止variable与v叠放后引起歧义，在此用phi代替variable

	//拓展向量
	MatMult(nonstagger_extend_u_matrix333333, velocity_u_at_p, velocity_u_at_p333333);
	MatMult(nonstagger_extend_v_matrix333333, velocity_v_at_p, velocity_v_at_p333333);
	MatMult(nonstagger_extend_w_matrix333333, velocity_w_at_p, velocity_w_at_p333333);
	MatMult(extend_phi_matrix333333, variable, variable333333);
	VecAXPY(velocity_u_at_p333333, 1.0, u_auxiliary_vector333333);
	//VecAXPY(velocity_v_at_p333333, 1.0, v_auxiliary_vector333333);
	//VecAXPY(velocity_w_at_p333333, 1.0, w_auxiliary_vector333333);
	//VecAXPY(variable333333, 1.0, phi_auxiliary_vector333333);

	//获得PN
	double max_abs_velocity_u;
	double max_abs_velocity_v;
	double max_abs_velocity_w;

	max_abs_velocity_u = get_max_abs(velocity_u_at_p333333);
	max_abs_velocity_v = get_max_abs(velocity_v_at_p333333);
	max_abs_velocity_w = get_max_abs(velocity_w_at_p333333);

	VecPointwiseMult(f_variable333333, velocity_u_at_p333333, variable333333);
	VecWAXPY(uphi_f_P333333, 1.0*max_abs_velocity_u, variable333333, f_variable333333);
	VecWAXPY(uphi_f_N333333, -1.0*max_abs_velocity_u, variable333333, f_variable333333);
	VecScale(uphi_f_P333333, 0.5);
	VecScale(uphi_f_N333333, 0.5);

	VecPointwiseMult(f_variable333333, velocity_v_at_p333333, variable333333);
	VecWAXPY(vphi_f_P333333, 1.0*max_abs_velocity_v, variable333333, f_variable333333);
	VecWAXPY(vphi_f_N333333, -1.0*max_abs_velocity_v, variable333333, f_variable333333);
	VecScale(vphi_f_P333333, 0.5);
	VecScale(vphi_f_N333333, 0.5);

	VecPointwiseMult(f_variable333333, velocity_w_at_p333333, variable333333);
	VecWAXPY(wphi_f_P333333, 1.0*max_abs_velocity_w, variable333333, f_variable333333);
	VecWAXPY(wphi_f_N333333, -1.0*max_abs_velocity_w, variable333333, f_variable333333);
	VecScale(wphi_f_P333333, 0.5);
	VecScale(wphi_f_N333333, 0.5);

	//求导
	get_WENO_d_variable_x_value(uphi_f_P333333, uphi_f_N333333, d_uphi_x);
	get_WENO_d_variable_y_value(vphi_f_P333333, vphi_f_N333333, d_vphi_y);
	get_WENO_d_variable_z_value(wphi_f_P333333, wphi_f_N333333, variable_LSeq_right_term);

	//获得右端项
	VecAXPBYPCZ(variable_LSeq_right_term, -1.0, -1.0, -1.0, d_uphi_x, d_vphi_y);
}

//FF输运方程右端项--------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_FF_LSeq_right_term3D(
	Vec variable00, Vec variable01, Vec variable02,
	Vec variable10, Vec variable11, Vec variable12,
	Vec variable20, Vec variable21, Vec variable22,
	Vec variable_right_term00, Vec variable_right_term01, Vec variable_right_term02,
	Vec variable_right_term10, Vec variable_right_term11, Vec variable_right_term12,
	Vec variable_right_term20, Vec variable_right_term21, Vec variable_right_term22)
{
	//00
	get_variable_LSeq_right_term3D(variable00, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_u_x_at_p, variable00);
	VecPointwiseMult(right_term_part1_1, d_u_y_at_p, variable10);
	VecPointwiseMult(right_term_part1, d_u_z_at_p, variable20);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term00, 1.0, right_term_part0, right_term_part1);

	//01
	get_variable_LSeq_right_term3D(variable01, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_u_x_at_p, variable01);
	VecPointwiseMult(right_term_part1_1, d_u_y_at_p, variable11);
	VecPointwiseMult(right_term_part1, d_u_z_at_p, variable21);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term01, 1.0, right_term_part0, right_term_part1);

	//02
	get_variable_LSeq_right_term3D(variable02, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_u_x_at_p, variable02);
	VecPointwiseMult(right_term_part1_1, d_u_y_at_p, variable12);
	VecPointwiseMult(right_term_part1, d_u_z_at_p, variable22);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term02, 1.0, right_term_part0, right_term_part1);

	//10
	get_variable_LSeq_right_term3D(variable10, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_v_x_at_p, variable00);
	VecPointwiseMult(right_term_part1_1, d_v_y_at_p, variable10);
	VecPointwiseMult(right_term_part1, d_v_z_at_p, variable20);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term10, 1.0, right_term_part0, right_term_part1);

	//11
	get_variable_LSeq_right_term3D(variable11, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_v_x_at_p, variable01);
	VecPointwiseMult(right_term_part1_1, d_v_y_at_p, variable11);
	VecPointwiseMult(right_term_part1, d_v_z_at_p, variable21);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term11, 1.0, right_term_part0, right_term_part1);

	//12
	get_variable_LSeq_right_term3D(variable12, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_v_x_at_p, variable02);
	VecPointwiseMult(right_term_part1_1, d_v_y_at_p, variable12);
	VecPointwiseMult(right_term_part1, d_v_z_at_p, variable22);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term12, 1.0, right_term_part0, right_term_part1);

	//20
	get_variable_LSeq_right_term3D(variable20, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_w_x_at_p, variable00);
	VecPointwiseMult(right_term_part1_1, d_w_y_at_p, variable10);
	VecPointwiseMult(right_term_part1, d_w_z_at_p, variable20);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term20, 1.0, right_term_part0, right_term_part1);

	//21
	get_variable_LSeq_right_term3D(variable21, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_w_x_at_p, variable01);
	VecPointwiseMult(right_term_part1_1, d_w_y_at_p, variable11);
	VecPointwiseMult(right_term_part1, d_w_z_at_p, variable21);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term21, 1.0, right_term_part0, right_term_part1);

	//22
	get_variable_LSeq_right_term3D(variable22, right_term_part0);

	VecPointwiseMult(right_term_part1_0, d_w_x_at_p, variable02);
	VecPointwiseMult(right_term_part1_1, d_w_y_at_p, variable12);
	VecPointwiseMult(right_term_part1, d_w_z_at_p, variable22);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);

	VecWAXPY(variable_right_term22, 1.0, right_term_part0, right_term_part1);
}

//FF_inv输运方程右端项---------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_FF_inv_LSeq_right_term3D(
	Vec variable00, Vec variable01, Vec variable02,
	Vec variable10, Vec variable11, Vec variable12,
	Vec variable20, Vec variable21, Vec variable22,
	Vec variable_right_term00, Vec variable_right_term01, Vec variable_right_term02,
	Vec variable_right_term10, Vec variable_right_term11, Vec variable_right_term12,
	Vec variable_right_term20, Vec variable_right_term21, Vec variable_right_term22)
{
	//00
	get_variable_LSeq_right_term3D(variable00, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable00, d_u_x_at_p);
	VecPointwiseMult(right_term_part1_1, variable01, d_v_x_at_p);
	VecPointwiseMult(right_term_part1, variable02, d_w_x_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term00, 1.0, right_term_part0, right_term_part1);

	//01
	get_variable_LSeq_right_term3D(variable01, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable00, d_u_y_at_p);
	VecPointwiseMult(right_term_part1_1, variable01, d_v_y_at_p);
	VecPointwiseMult(right_term_part1, variable02, d_w_y_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term01, 1.0, right_term_part0, right_term_part1);

	//02
	get_variable_LSeq_right_term3D(variable02, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable00, d_u_z_at_p);
	VecPointwiseMult(right_term_part1_1, variable01, d_v_z_at_p);
	VecPointwiseMult(right_term_part1, variable02, d_w_z_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term02, 1.0, right_term_part0, right_term_part1);

	//10
	get_variable_LSeq_right_term3D(variable10, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable10, d_u_x_at_p);
	VecPointwiseMult(right_term_part1_1, variable11, d_v_x_at_p);
	VecPointwiseMult(right_term_part1, variable12, d_w_x_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term10, 1.0, right_term_part0, right_term_part1);

	//11
	get_variable_LSeq_right_term3D(variable11, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable10, d_u_y_at_p);
	VecPointwiseMult(right_term_part1_1, variable11, d_v_y_at_p);
	VecPointwiseMult(right_term_part1, variable12, d_w_y_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term11, 1.0, right_term_part0, right_term_part1);

	//12
	get_variable_LSeq_right_term3D(variable12, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable10, d_u_z_at_p);
	VecPointwiseMult(right_term_part1_1, variable11, d_v_z_at_p);
	VecPointwiseMult(right_term_part1, variable12, d_w_z_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term12, 1.0, right_term_part0, right_term_part1);

	//20
	get_variable_LSeq_right_term3D(variable20, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable20, d_u_x_at_p);
	VecPointwiseMult(right_term_part1_1, variable21, d_v_x_at_p);
	VecPointwiseMult(right_term_part1, variable22, d_w_x_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term20, 1.0, right_term_part0, right_term_part1);

	//21
	get_variable_LSeq_right_term3D(variable21, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable20, d_u_y_at_p);
	VecPointwiseMult(right_term_part1_1, variable21, d_v_y_at_p);
	VecPointwiseMult(right_term_part1, variable22, d_w_y_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term21, 1.0, right_term_part0, right_term_part1);

	//22
	get_variable_LSeq_right_term3D(variable22, right_term_part0);

	VecPointwiseMult(right_term_part1_0, variable20, d_u_z_at_p);
	VecPointwiseMult(right_term_part1_1, variable21, d_v_z_at_p);
	VecPointwiseMult(right_term_part1, variable22, d_w_z_at_p);
	VecAXPBYPCZ(right_term_part1, 1.0, 1.0, 1.0, right_term_part1_0, right_term_part1_1);
	VecScale(right_term_part1, -1.0);

	VecWAXPY(variable_right_term22, 1.0, right_term_part0, right_term_part1);
}

//更新phi----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::update_phi3D()
{
	//TVD01
	get_variable_LSeq_right_term3D(phi, TVD_phi_LSeq_right_term0);
	VecWAXPY(TVD_phi1, dt, TVD_phi_LSeq_right_term0, phi);

	//TVD02
	get_variable_LSeq_right_term3D(TVD_phi1, TVD_phi_LSeq_right_term1);
	VecWAXPY(TVD_phi2, -3.0 / 4.0 * dt, TVD_phi_LSeq_right_term0, TVD_phi1);
	VecAXPY(TVD_phi2, 1.0 / 4.0 * dt, TVD_phi_LSeq_right_term1);

	//TVD03
	get_variable_LSeq_right_term3D(TVD_phi2, TVD_phi_LSeq_right_term2);
	VecWAXPY(phi, -1.0 / 12.0 * dt, TVD_phi_LSeq_right_term0, TVD_phi2);
	VecAXPY(phi, -1.0 / 12.0 * dt, TVD_phi_LSeq_right_term1);
	VecAXPY(phi, 8.0 / 12.0 * dt, TVD_phi_LSeq_right_term2);
}

//更新FF-----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::update_FF3D()
{
	//TVD0
	get_FF_LSeq_right_term3D(
		FF00, FF01, FF02,
		FF10, FF11, FF12,
		FF20, FF21, FF22,
		TVD_FF00_right_term0, TVD_FF01_right_term0, TVD_FF02_right_term0,
		TVD_FF10_right_term0, TVD_FF11_right_term0, TVD_FF12_right_term0,
		TVD_FF20_right_term0, TVD_FF21_right_term0, TVD_FF22_right_term0);
	VecWAXPY(TVD_FF00_1, dt, TVD_FF00_right_term0, FF00);
	VecWAXPY(TVD_FF01_1, dt, TVD_FF01_right_term0, FF01);
	VecWAXPY(TVD_FF02_1, dt, TVD_FF02_right_term0, FF02);
	VecWAXPY(TVD_FF10_1, dt, TVD_FF10_right_term0, FF10);
	VecWAXPY(TVD_FF11_1, dt, TVD_FF11_right_term0, FF11);
	VecWAXPY(TVD_FF12_1, dt, TVD_FF12_right_term0, FF12);
	VecWAXPY(TVD_FF20_1, dt, TVD_FF20_right_term0, FF20);
	VecWAXPY(TVD_FF21_1, dt, TVD_FF21_right_term0, FF21);
	VecWAXPY(TVD_FF22_1, dt, TVD_FF22_right_term0, FF22);

	//TVD1
	get_FF_LSeq_right_term3D(
		TVD_FF00_1, TVD_FF01_1, TVD_FF02_1,
		TVD_FF10_1, TVD_FF11_1, TVD_FF12_1,
		TVD_FF20_1, TVD_FF21_1, TVD_FF22_1,
		TVD_FF00_right_term1, TVD_FF01_right_term1, TVD_FF02_right_term1,
		TVD_FF10_right_term1, TVD_FF11_right_term1, TVD_FF12_right_term1,
		TVD_FF20_right_term1, TVD_FF21_right_term1, TVD_FF22_right_term1);
	VecWAXPY(TVD_FF00_2, -3.0 / 4.0 * dt, TVD_FF00_right_term0, TVD_FF00_1);
	VecWAXPY(TVD_FF01_2, -3.0 / 4.0 * dt, TVD_FF01_right_term0, TVD_FF01_1);
	VecWAXPY(TVD_FF02_2, -3.0 / 4.0 * dt, TVD_FF02_right_term0, TVD_FF02_1);
	VecWAXPY(TVD_FF10_2, -3.0 / 4.0 * dt, TVD_FF10_right_term0, TVD_FF10_1);
	VecWAXPY(TVD_FF11_2, -3.0 / 4.0 * dt, TVD_FF11_right_term0, TVD_FF11_1);
	VecWAXPY(TVD_FF12_2, -3.0 / 4.0 * dt, TVD_FF12_right_term0, TVD_FF12_1);
	VecWAXPY(TVD_FF20_2, -3.0 / 4.0 * dt, TVD_FF20_right_term0, TVD_FF20_1);
	VecWAXPY(TVD_FF21_2, -3.0 / 4.0 * dt, TVD_FF21_right_term0, TVD_FF21_1);
	VecWAXPY(TVD_FF22_2, -3.0 / 4.0 * dt, TVD_FF22_right_term0, TVD_FF22_1);

	VecAXPY(TVD_FF00_2, 1.0 / 4.0 * dt, TVD_FF00_right_term1);
	VecAXPY(TVD_FF01_2, 1.0 / 4.0 * dt, TVD_FF01_right_term1);
	VecAXPY(TVD_FF02_2, 1.0 / 4.0 * dt, TVD_FF02_right_term1);
	VecAXPY(TVD_FF10_2, 1.0 / 4.0 * dt, TVD_FF10_right_term1);
	VecAXPY(TVD_FF11_2, 1.0 / 4.0 * dt, TVD_FF11_right_term1);
	VecAXPY(TVD_FF12_2, 1.0 / 4.0 * dt, TVD_FF12_right_term1);
	VecAXPY(TVD_FF20_2, 1.0 / 4.0 * dt, TVD_FF20_right_term1);
	VecAXPY(TVD_FF21_2, 1.0 / 4.0 * dt, TVD_FF21_right_term1);
	VecAXPY(TVD_FF22_2, 1.0 / 4.0 * dt, TVD_FF22_right_term1);

	//TVD2
	get_FF_LSeq_right_term3D(
		TVD_FF00_2, TVD_FF01_2, TVD_FF02_2,
		TVD_FF10_2, TVD_FF11_2, TVD_FF12_2,
		TVD_FF20_2, TVD_FF21_2, TVD_FF22_2,
		TVD_FF00_right_term2, TVD_FF01_right_term2, TVD_FF02_right_term2,
		TVD_FF10_right_term2, TVD_FF11_right_term2, TVD_FF12_right_term2,
		TVD_FF20_right_term2, TVD_FF21_right_term2, TVD_FF22_right_term2);
	VecWAXPY(FF00, -1.0 / 12.0 * dt, TVD_FF00_right_term0, TVD_FF00_2);
	VecWAXPY(FF01, -1.0 / 12.0 * dt, TVD_FF01_right_term0, TVD_FF01_2);
	VecWAXPY(FF02, -1.0 / 12.0 * dt, TVD_FF02_right_term0, TVD_FF02_2);
	VecWAXPY(FF10, -1.0 / 12.0 * dt, TVD_FF10_right_term0, TVD_FF10_2);
	VecWAXPY(FF11, -1.0 / 12.0 * dt, TVD_FF11_right_term0, TVD_FF11_2);
	VecWAXPY(FF12, -1.0 / 12.0 * dt, TVD_FF12_right_term0, TVD_FF12_2);
	VecWAXPY(FF20, -1.0 / 12.0 * dt, TVD_FF20_right_term0, TVD_FF20_2);
	VecWAXPY(FF21, -1.0 / 12.0 * dt, TVD_FF21_right_term0, TVD_FF21_2);
	VecWAXPY(FF22, -1.0 / 12.0 * dt, TVD_FF22_right_term0, TVD_FF22_2);

	VecAXPY(FF00, -1.0 / 12.0 * dt, TVD_FF00_right_term1);
	VecAXPY(FF01, -1.0 / 12.0 * dt, TVD_FF01_right_term1);
	VecAXPY(FF02, -1.0 / 12.0 * dt, TVD_FF02_right_term1);
	VecAXPY(FF10, -1.0 / 12.0 * dt, TVD_FF10_right_term1);
	VecAXPY(FF11, -1.0 / 12.0 * dt, TVD_FF11_right_term1);
	VecAXPY(FF12, -1.0 / 12.0 * dt, TVD_FF12_right_term1);
	VecAXPY(FF20, -1.0 / 12.0 * dt, TVD_FF20_right_term1);
	VecAXPY(FF21, -1.0 / 12.0 * dt, TVD_FF21_right_term1);
	VecAXPY(FF22, -1.0 / 12.0 * dt, TVD_FF22_right_term1);

	VecAXPY(FF00, 8.0 / 12.0 * dt, TVD_FF00_right_term2);
	VecAXPY(FF01, 8.0 / 12.0 * dt, TVD_FF01_right_term2);
	VecAXPY(FF02, 8.0 / 12.0 * dt, TVD_FF02_right_term2);
	VecAXPY(FF10, 8.0 / 12.0 * dt, TVD_FF10_right_term2);
	VecAXPY(FF11, 8.0 / 12.0 * dt, TVD_FF11_right_term2);
	VecAXPY(FF12, 8.0 / 12.0 * dt, TVD_FF12_right_term2);
	VecAXPY(FF20, 8.0 / 12.0 * dt, TVD_FF20_right_term2);
	VecAXPY(FF21, 8.0 / 12.0 * dt, TVD_FF21_right_term2);
	VecAXPY(FF22, 8.0 / 12.0 * dt, TVD_FF22_right_term2);

	//修正形变梯度
	VecPointwiseMult(FF00, FF00, extract_interior_domain_vector);
	VecPointwiseMult(FF01, FF01, extract_interior_domain_vector);
	VecPointwiseMult(FF02, FF02, extract_interior_domain_vector);
	VecPointwiseMult(FF10, FF10, extract_interior_domain_vector);
	VecPointwiseMult(FF11, FF11, extract_interior_domain_vector);
	VecPointwiseMult(FF12, FF12, extract_interior_domain_vector);
	VecPointwiseMult(FF20, FF20, extract_interior_domain_vector);
	VecPointwiseMult(FF21, FF21, extract_interior_domain_vector);
	VecPointwiseMult(FF22, FF22, extract_interior_domain_vector);
}

//更新FF_inv-------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::update_FF_inv3D()
{
	//TVD0
	get_FF_inv_LSeq_right_term3D(
		FF_inv00, FF_inv01, FF_inv02,
		FF_inv10, FF_inv11, FF_inv12,
		FF_inv20, FF_inv21, FF_inv22,
		TVD_FF_inv00_right_term0, TVD_FF_inv01_right_term0, TVD_FF_inv02_right_term0,
		TVD_FF_inv10_right_term0, TVD_FF_inv11_right_term0, TVD_FF_inv12_right_term0,
		TVD_FF_inv20_right_term0, TVD_FF_inv21_right_term0, TVD_FF_inv22_right_term0);
	VecWAXPY(TVD_FF_inv00_1, dt, TVD_FF_inv00_right_term0, FF_inv00);
	VecWAXPY(TVD_FF_inv01_1, dt, TVD_FF_inv01_right_term0, FF_inv01);
	VecWAXPY(TVD_FF_inv02_1, dt, TVD_FF_inv02_right_term0, FF_inv02);
	VecWAXPY(TVD_FF_inv10_1, dt, TVD_FF_inv10_right_term0, FF_inv10);
	VecWAXPY(TVD_FF_inv11_1, dt, TVD_FF_inv11_right_term0, FF_inv11);
	VecWAXPY(TVD_FF_inv12_1, dt, TVD_FF_inv12_right_term0, FF_inv12);
	VecWAXPY(TVD_FF_inv20_1, dt, TVD_FF_inv20_right_term0, FF_inv20);
	VecWAXPY(TVD_FF_inv21_1, dt, TVD_FF_inv21_right_term0, FF_inv21);
	VecWAXPY(TVD_FF_inv22_1, dt, TVD_FF_inv22_right_term0, FF_inv22);

	//TVD1
	get_FF_inv_LSeq_right_term3D(
		TVD_FF_inv00_1, TVD_FF_inv01_1, TVD_FF_inv02_1,
		TVD_FF_inv10_1, TVD_FF_inv11_1, TVD_FF_inv12_1,
		TVD_FF_inv20_1, TVD_FF_inv21_1, TVD_FF_inv22_1,
		TVD_FF_inv00_right_term1, TVD_FF_inv01_right_term1, TVD_FF_inv02_right_term1,
		TVD_FF_inv10_right_term1, TVD_FF_inv11_right_term1, TVD_FF_inv12_right_term1,
		TVD_FF_inv20_right_term1, TVD_FF_inv21_right_term1, TVD_FF_inv22_right_term1);
	VecWAXPY(TVD_FF_inv00_2, -3.0 / 4.0 * dt, TVD_FF_inv00_right_term0, TVD_FF_inv00_1);
	VecWAXPY(TVD_FF_inv01_2, -3.0 / 4.0 * dt, TVD_FF_inv01_right_term0, TVD_FF_inv01_1);
	VecWAXPY(TVD_FF_inv02_2, -3.0 / 4.0 * dt, TVD_FF_inv02_right_term0, TVD_FF_inv02_1);
	VecWAXPY(TVD_FF_inv10_2, -3.0 / 4.0 * dt, TVD_FF_inv10_right_term0, TVD_FF_inv10_1);
	VecWAXPY(TVD_FF_inv11_2, -3.0 / 4.0 * dt, TVD_FF_inv11_right_term0, TVD_FF_inv11_1);
	VecWAXPY(TVD_FF_inv12_2, -3.0 / 4.0 * dt, TVD_FF_inv12_right_term0, TVD_FF_inv12_1);
	VecWAXPY(TVD_FF_inv20_2, -3.0 / 4.0 * dt, TVD_FF_inv20_right_term0, TVD_FF_inv20_1);
	VecWAXPY(TVD_FF_inv21_2, -3.0 / 4.0 * dt, TVD_FF_inv21_right_term0, TVD_FF_inv21_1);
	VecWAXPY(TVD_FF_inv22_2, -3.0 / 4.0 * dt, TVD_FF_inv22_right_term0, TVD_FF_inv22_1);

	VecAXPY(TVD_FF_inv00_2, 1.0 / 4.0 * dt, TVD_FF_inv00_right_term1);
	VecAXPY(TVD_FF_inv01_2, 1.0 / 4.0 * dt, TVD_FF_inv01_right_term1);
	VecAXPY(TVD_FF_inv02_2, 1.0 / 4.0 * dt, TVD_FF_inv02_right_term1);
	VecAXPY(TVD_FF_inv10_2, 1.0 / 4.0 * dt, TVD_FF_inv10_right_term1);
	VecAXPY(TVD_FF_inv11_2, 1.0 / 4.0 * dt, TVD_FF_inv11_right_term1);
	VecAXPY(TVD_FF_inv12_2, 1.0 / 4.0 * dt, TVD_FF_inv12_right_term1);
	VecAXPY(TVD_FF_inv20_2, 1.0 / 4.0 * dt, TVD_FF_inv20_right_term1);
	VecAXPY(TVD_FF_inv21_2, 1.0 / 4.0 * dt, TVD_FF_inv21_right_term1);
	VecAXPY(TVD_FF_inv22_2, 1.0 / 4.0 * dt, TVD_FF_inv22_right_term1);

	//TVD2
	get_FF_inv_LSeq_right_term3D(
		TVD_FF_inv00_2, TVD_FF_inv01_2, TVD_FF_inv02_2,
		TVD_FF_inv10_2, TVD_FF_inv11_2, TVD_FF_inv12_2,
		TVD_FF_inv20_2, TVD_FF_inv21_2, TVD_FF_inv22_2,
		TVD_FF_inv00_right_term2, TVD_FF_inv01_right_term2, TVD_FF_inv02_right_term2,
		TVD_FF_inv10_right_term2, TVD_FF_inv11_right_term2, TVD_FF_inv12_right_term2,
		TVD_FF_inv20_right_term2, TVD_FF_inv21_right_term2, TVD_FF_inv22_right_term2);
	VecWAXPY(FF_inv00, -1.0 / 12.0 * dt, TVD_FF_inv00_right_term0, TVD_FF_inv00_2);
	VecWAXPY(FF_inv01, -1.0 / 12.0 * dt, TVD_FF_inv01_right_term0, TVD_FF_inv01_2);
	VecWAXPY(FF_inv02, -1.0 / 12.0 * dt, TVD_FF_inv02_right_term0, TVD_FF_inv02_2);
	VecWAXPY(FF_inv10, -1.0 / 12.0 * dt, TVD_FF_inv10_right_term0, TVD_FF_inv10_2);
	VecWAXPY(FF_inv11, -1.0 / 12.0 * dt, TVD_FF_inv11_right_term0, TVD_FF_inv11_2);
	VecWAXPY(FF_inv12, -1.0 / 12.0 * dt, TVD_FF_inv12_right_term0, TVD_FF_inv12_2);
	VecWAXPY(FF_inv20, -1.0 / 12.0 * dt, TVD_FF_inv20_right_term0, TVD_FF_inv20_2);
	VecWAXPY(FF_inv21, -1.0 / 12.0 * dt, TVD_FF_inv21_right_term0, TVD_FF_inv21_2);
	VecWAXPY(FF_inv22, -1.0 / 12.0 * dt, TVD_FF_inv22_right_term0, TVD_FF_inv22_2);

	VecAXPY(FF_inv00, -1.0 / 12.0 * dt, TVD_FF_inv00_right_term1);
	VecAXPY(FF_inv01, -1.0 / 12.0 * dt, TVD_FF_inv01_right_term1);
	VecAXPY(FF_inv02, -1.0 / 12.0 * dt, TVD_FF_inv02_right_term1);
	VecAXPY(FF_inv10, -1.0 / 12.0 * dt, TVD_FF_inv10_right_term1);
	VecAXPY(FF_inv11, -1.0 / 12.0 * dt, TVD_FF_inv11_right_term1);
	VecAXPY(FF_inv12, -1.0 / 12.0 * dt, TVD_FF_inv12_right_term1);
	VecAXPY(FF_inv20, -1.0 / 12.0 * dt, TVD_FF_inv20_right_term1);
	VecAXPY(FF_inv21, -1.0 / 12.0 * dt, TVD_FF_inv21_right_term1);
	VecAXPY(FF_inv22, -1.0 / 12.0 * dt, TVD_FF_inv22_right_term1);

	VecAXPY(FF_inv00, 8.0 / 12.0 * dt, TVD_FF_inv00_right_term2);
	VecAXPY(FF_inv01, 8.0 / 12.0 * dt, TVD_FF_inv01_right_term2);
	VecAXPY(FF_inv02, 8.0 / 12.0 * dt, TVD_FF_inv02_right_term2);
	VecAXPY(FF_inv10, 8.0 / 12.0 * dt, TVD_FF_inv10_right_term2);
	VecAXPY(FF_inv11, 8.0 / 12.0 * dt, TVD_FF_inv11_right_term2);
	VecAXPY(FF_inv12, 8.0 / 12.0 * dt, TVD_FF_inv12_right_term2);
	VecAXPY(FF_inv20, 8.0 / 12.0 * dt, TVD_FF_inv20_right_term2);
	VecAXPY(FF_inv21, 8.0 / 12.0 * dt, TVD_FF_inv21_right_term2);
	VecAXPY(FF_inv22, 8.0 / 12.0 * dt, TVD_FF_inv22_right_term2);

	//修正形变梯度
	VecPointwiseMult(FF_inv00, FF_inv00, extract_interior_domain_vector);
	VecPointwiseMult(FF_inv01, FF_inv01, extract_interior_domain_vector);
	VecPointwiseMult(FF_inv02, FF_inv02, extract_interior_domain_vector);
	VecPointwiseMult(FF_inv10, FF_inv10, extract_interior_domain_vector);
	VecPointwiseMult(FF_inv11, FF_inv11, extract_interior_domain_vector);
	VecPointwiseMult(FF_inv12, FF_inv12, extract_interior_domain_vector);
	VecPointwiseMult(FF_inv20, FF_inv20, extract_interior_domain_vector);
	VecPointwiseMult(FF_inv21, FF_inv21, extract_interior_domain_vector);
	VecPointwiseMult(FF_inv22, FF_inv22, extract_interior_domain_vector);
}

//modified phi-----------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::modified_phi_based_on_intermediate_step()
{
	//定义变量----------------
	double account_num = 0;

	double disturbance = 1.e-4;
	double d_tau = 0.0005*dt;
	double hh = 1e-8;

	double * phi_temp;

	//赋值-------------------

	//全一向量
	VecSet(ones_vector, 1.0);

	//组装全局符号、边界标号、修正区域标号
	VecSet(phi_boundary_mark_vector, 0.0);

	VecGetArray(phi, &phi_temp);
	VecGetOwnershipRange(phi, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		double phi_value_temp = phi_temp[Ii - istart];

		double boundary_vector_value;
		
		//划定区域
		if (phi_value_temp > 0 + hh_phi && phi_value_temp < 1 - hh_phi){
			boundary_vector_value = 1.0;
			VecSetValues(phi_boundary_mark_vector, 1, &Ii, &boundary_vector_value, INSERT_VALUES);
		}
	}
	VecRestoreArray(phi, &phi_temp);

	VecAssemblyBegin(phi_boundary_mark_vector);
	VecAssemblyEnd(phi_boundary_mark_vector);
	
	//计算法线向量
	MatMult(extend_phi_matrix111111, phi, phi_at_p111111);
	MatMult(dx_matrix111111_to_000000, phi_at_p111111, phi_normal_x);
	MatMult(dy_matrix111111_to_000000, phi_at_p111111, phi_normal_y);
	MatMult(dz_matrix111111_to_000000, phi_at_p111111, phi_normal_z);

	VecPointwiseMult(square_phi_normal_x, phi_normal_x, phi_normal_x);
	VecPointwiseMult(square_phi_normal_y, phi_normal_y, phi_normal_y);
	VecPointwiseMult(normal_norm, phi_normal_z, phi_normal_z);

	VecAXPBYPCZ(normal_norm, 1.0, 1.0, 1.0, square_phi_normal_x, square_phi_normal_y);
	VecSqrtAbs(normal_norm);
	VecShift(normal_norm, 1e-10);
	VecReciprocal(normal_norm);

	VecPointwiseMult(normal_norm, normal_norm, extract_boundary_domain_vector);
	VecPointwiseMult(phi_normal_x, phi_normal_x, normal_norm);
	VecPointwiseMult(phi_normal_y, phi_normal_y, normal_norm);
	VecPointwiseMult(phi_normal_z, phi_normal_z, normal_norm);

	VecPointwiseMult(phi_normal_x, phi_normal_x, phi_boundary_mark_vector);
	VecPointwiseMult(phi_normal_y, phi_normal_y, phi_boundary_mark_vector);
	VecPointwiseMult(phi_normal_z, phi_normal_z, phi_boundary_mark_vector);
	
	for (int i = 0; i < 3; i++){
		//1
		MatMult(extend_phi_matrix111111, phi, phi_at_p111111);
		//VecAXPY(phi_at_p1111, 1.0, phi_auxiliary_vector1111);
		MatMult(dx_matrix111111_to_000000, phi_at_p111111, d_phi_x);
		MatMult(dy_matrix111111_to_000000, phi_at_p111111, d_phi_y);
		MatMult(dz_matrix111111_to_000000, phi_at_p111111, d_phi_z);
	
		//2
		VecWAXPY(one_redu_phi, -1.0, phi, ones_vector);
		VecPointwiseMult(modified_part1, phi, one_redu_phi);

		VecPointwiseMult(modified_part2, modified_part1, phi_normal_x);
		MatMult(extend_phi_matrix111111, modified_part2, phi_at_p111111);
		//VecAXPY(phi_at_p111111, 1.0, phi_auxiliary_vector111111);
		MatMult(dx_matrix111111_to_000000, phi_at_p111111, modified_part_x);
		VecAXPY(modified_part_x, -disturbance, d_phi_x);

		VecPointwiseMult(modified_part2, modified_part1, phi_normal_y);
		MatMult(extend_phi_matrix111111, modified_part2, phi_at_p111111);
		//VecAXPY(phi_at_p111111, 1.0, phi_auxiliary_vector111111);
		MatMult(dy_matrix111111_to_000000, phi_at_p111111, modified_part_y);
		VecAXPY(modified_part_y, -disturbance, d_phi_y);

		VecPointwiseMult(modified_part2, modified_part1, phi_normal_z);
		MatMult(extend_phi_matrix111111, modified_part2, phi_at_p111111);
		//VecAXPY(phi_at_p111111, 1.0, phi_auxiliary_vector111111);
		MatMult(dz_matrix111111_to_000000, phi_at_p111111, modified_part_z);
		VecAXPY(modified_part_z, -disturbance, d_phi_z);
	
		//右端项
		VecWAXPY(modified_right_term, 1.0, modified_part_x, modified_part_y);
		VecAXPY(modified_right_term, 1.0, modified_part_z);
		VecScale(modified_right_term, -1.0);
		VecPointwiseMult(modified_right_term, modified_right_term, phi_boundary_mark_vector);
	
		//更新
		VecAXPY(phi, d_tau, modified_right_term);
	}

	delete[] phi_temp;
}

//法线向量=========================================================================================
void wyh_solve_fluid_structure_interaction_system3D200630::solve_phi_normal()
{
	//边界范围
	assembly_extract_boundary_domain_vector3D(extract_boundary_domain_vector);

	//拓展phi
	MatMult(extend_phi_matrix111111, phi, phi_at_p111111);
	//VecAXPY(phi_at_p111111, 1.0, phi_auxiliary_vector111111);

	MatMult(dx_matrix111111_to_000000, phi_at_p111111, phi_normal_x);
	MatMult(dy_matrix111111_to_000000, phi_at_p111111, phi_normal_y);
	MatMult(dz_matrix111111_to_000000, phi_at_p111111, phi_normal_z);

	VecPointwiseMult(square_phi_normal_x, phi_normal_x, phi_normal_x);
	VecPointwiseMult(square_phi_normal_y, phi_normal_y, phi_normal_y);
	VecPointwiseMult(normal_norm, phi_normal_z, phi_normal_z);

	VecAXPBYPCZ(normal_norm, 1.0, 1.0, 1.0, square_phi_normal_x, square_phi_normal_y);
	VecSqrtAbs(normal_norm);
	VecShift(normal_norm, 1e-10);
	VecReciprocal(normal_norm);

	VecPointwiseMult(normal_norm, normal_norm, extract_boundary_domain_vector);
	VecPointwiseMult(phi_normal_x, phi_normal_x, normal_norm);
	VecPointwiseMult(phi_normal_y, phi_normal_y, normal_norm);
	VecPointwiseMult(phi_normal_z, phi_normal_z, normal_norm);
}

//***************************************************************************************************************************************************

//计算WENO===========================================================================================================================================

//***************************************************************************************************************************************************

//速度绝对值的最大值-----------------------------
double wyh_solve_fluid_structure_interaction_system3D200630::get_max_abs(Vec variable)
{
	double max_variable_abs;
	
	VecCopy(variable, variable_abs);
	VecAbs(variable_abs);
	VecMax(variable_abs, NULL, &max_variable_abs);
	
	return max_variable_abs;
}

//WENO at x====================================================================

//一阶导数---------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_WENO_d_variable_x_value(
	Vec variable_f_P, Vec variable_f_N, Vec d_variable_x)
{
	//权重
	get_WeightFP_at_x(variable_f_P);
	get_WeightFN_at_x(variable_f_N);

	//分量
	get_CellBoundaryFP_x(variable_f_P);
	get_CellBoundaryFN_x(variable_f_N);

	//
	VecPointwiseMult(FP_at_x0, OmegaP0_at_x, f_P0_at_x);	VecPointwiseMult(FN_at_x0, OmegaN0_at_x, f_N0_at_x);
	VecPointwiseMult(FP_at_x1, OmegaP1_at_x, f_P1_at_x);	VecPointwiseMult(FN_at_x1, OmegaN1_at_x, f_N1_at_x);
	VecPointwiseMult(FP_at_x2, OmegaP2_at_x, f_P2_at_x);	VecPointwiseMult(FN_at_x2, OmegaN2_at_x, f_N2_at_x);

	VecWAXPY(FP_at_x, 1.0, FP_at_x0, FP_at_x1);
	VecAXPY(FP_at_x, 1.0, FP_at_x2);
	VecWAXPY(FN_at_x, 1.0, FN_at_x0, FN_at_x1);
	VecAXPY(FN_at_x, 1.0, FN_at_x2);

	VecWAXPY(d_variable_x_temp, 1.0, FP_at_x, FN_at_x);

	MatMult(dx_matrix100000_or_010000_to_000000, d_variable_x_temp, d_variable_x);
}

//分量-------------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_CellBoundaryFP_x(Vec variable_f_P)
{
	MatMult(f_P0_x_matrix333333_to_100000, variable_f_P, f_P0_at_x);
	MatMult(f_P1_x_matrix333333_to_100000, variable_f_P, f_P1_at_x);
	MatMult(f_P2_x_matrix333333_to_100000, variable_f_P, f_P2_at_x);
}

void wyh_solve_fluid_structure_interaction_system3D200630::get_CellBoundaryFN_x(Vec variable_f_N)
{
	MatMult(f_N0_x_matrix333333_to_100000, variable_f_N, f_N0_at_x);
	MatMult(f_N1_x_matrix333333_to_100000, variable_f_N, f_N1_at_x);
	MatMult(f_N2_x_matrix333333_to_100000, variable_f_N, f_N2_at_x);
}

//权重-------------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_WeightFP_at_x(Vec variable_f_P)
{
	//part 1
	MatMult(ISP0_x_part1_matrix333333_to_100000, variable_f_P, ISP0_part1_value_at_x);
	MatMult(ISP1_x_part1_matrix333333_to_100000, variable_f_P, ISP1_part1_value_at_x);
	MatMult(ISP2_x_part1_matrix333333_to_100000, variable_f_P, ISP2_part1_value_at_x);

	VecPointwiseMult(ISP0_part1_value_at_x, ISP0_part1_value_at_x, ISP0_part1_value_at_x);
	VecPointwiseMult(ISP1_part1_value_at_x, ISP1_part1_value_at_x, ISP1_part1_value_at_x);
	VecPointwiseMult(ISP2_part1_value_at_x, ISP2_part1_value_at_x, ISP2_part1_value_at_x);

	VecScale(ISP0_part1_value_at_x, 13.0 / 12.0);
	VecScale(ISP1_part1_value_at_x, 13.0 / 12.0);
	VecScale(ISP2_part1_value_at_x, 13.0 / 12.0);

	//part 2
	MatMult(ISP0_x_part2_matrix333333_to_100000, variable_f_P, ISP0_part2_value_at_x);
	MatMult(ISP1_x_part2_matrix333333_to_100000, variable_f_P, ISP1_part2_value_at_x);
	MatMult(ISP2_x_part2_matrix333333_to_100000, variable_f_P, ISP2_part2_value_at_x);

	VecPointwiseMult(ISP0_part2_value_at_x, ISP0_part2_value_at_x, ISP0_part2_value_at_x);
	VecPointwiseMult(ISP1_part2_value_at_x, ISP1_part2_value_at_x, ISP1_part2_value_at_x);
	VecPointwiseMult(ISP2_part2_value_at_x, ISP2_part2_value_at_x, ISP2_part2_value_at_x);

	VecScale(ISP0_part2_value_at_x, 1.0 / 4.0);
	VecScale(ISP1_part2_value_at_x, 1.0 / 4.0);
	VecScale(ISP2_part2_value_at_x, 1.0 / 4.0);

	//求alpha
	double gammaP0 = 0.1;
	double gammaP1 = 0.6;
	double gammaP2 = 0.3;

	VecWAXPY(alphaP0_at_x, 1.0, ISP0_part1_value_at_x, ISP0_part2_value_at_x);
	VecWAXPY(alphaP1_at_x, 1.0, ISP1_part1_value_at_x, ISP1_part2_value_at_x);
	VecWAXPY(alphaP2_at_x, 1.0, ISP2_part1_value_at_x, ISP2_part2_value_at_x);

	VecShift(alphaP0_at_x, 1e-10);
	VecShift(alphaP1_at_x, 1e-10);
	VecShift(alphaP2_at_x, 1e-10);

	VecReciprocal(alphaP0_at_x);
	VecReciprocal(alphaP1_at_x);
	VecReciprocal(alphaP2_at_x);

	VecPointwiseMult(alphaP0_at_x, alphaP0_at_x, alphaP0_at_x);
	VecPointwiseMult(alphaP1_at_x, alphaP1_at_x, alphaP1_at_x);
	VecPointwiseMult(alphaP2_at_x, alphaP2_at_x, alphaP2_at_x);

	VecScale(alphaP0_at_x, gammaP0);
	VecScale(alphaP1_at_x, gammaP1);
	VecScale(alphaP2_at_x, gammaP2);

	//求Omega
	VecWAXPY(sum_alphaP_at_x, 1.0, alphaP0_at_x, alphaP1_at_x);
	VecAXPY(sum_alphaP_at_x, 1.0, alphaP2_at_x);

	VecShift(sum_alphaP_at_x, 1e-10);
	VecReciprocal(sum_alphaP_at_x);

	VecPointwiseMult(OmegaP0_at_x, alphaP0_at_x, sum_alphaP_at_x);
	VecPointwiseMult(OmegaP1_at_x, alphaP1_at_x, sum_alphaP_at_x);
	VecPointwiseMult(OmegaP2_at_x, alphaP2_at_x, sum_alphaP_at_x);
}

void wyh_solve_fluid_structure_interaction_system3D200630::get_WeightFN_at_x(Vec variable_f_N)
{
	//part 1
	MatMult(ISN0_x_part1_matrix333333_to_100000, variable_f_N, ISN0_part1_value_at_x);
	MatMult(ISN1_x_part1_matrix333333_to_100000, variable_f_N, ISN1_part1_value_at_x);
	MatMult(ISN2_x_part1_matrix333333_to_100000, variable_f_N, ISN2_part1_value_at_x);

	VecPointwiseMult(ISN0_part1_value_at_x, ISN0_part1_value_at_x, ISN0_part1_value_at_x);
	VecPointwiseMult(ISN1_part1_value_at_x, ISN1_part1_value_at_x, ISN1_part1_value_at_x);
	VecPointwiseMult(ISN2_part1_value_at_x, ISN2_part1_value_at_x, ISN2_part1_value_at_x);

	VecScale(ISN0_part1_value_at_x, 13.0 / 12.0);
	VecScale(ISN1_part1_value_at_x, 13.0 / 12.0);
	VecScale(ISN2_part1_value_at_x, 13.0 / 12.0);

	//part 2
	MatMult(ISN0_x_part2_matrix333333_to_100000, variable_f_N, ISN0_part2_value_at_x);
	MatMult(ISN1_x_part2_matrix333333_to_100000, variable_f_N, ISN1_part2_value_at_x);
	MatMult(ISN2_x_part2_matrix333333_to_100000, variable_f_N, ISN2_part2_value_at_x);

	VecPointwiseMult(ISN0_part2_value_at_x, ISN0_part2_value_at_x, ISN0_part2_value_at_x);
	VecPointwiseMult(ISN1_part2_value_at_x, ISN1_part2_value_at_x, ISN1_part2_value_at_x);
	VecPointwiseMult(ISN2_part2_value_at_x, ISN2_part2_value_at_x, ISN2_part2_value_at_x);

	VecScale(ISN0_part2_value_at_x, 1.0 / 4.0);
	VecScale(ISN1_part2_value_at_x, 1.0 / 4.0);
	VecScale(ISN2_part2_value_at_x, 1.0 / 4.0);

	//求alpha
	double gammaN0 = 0.3;
	double gammaN1 = 0.6;
	double gammaN2 = 0.1;

	VecWAXPY(alphaN0_at_x, 1.0, ISN0_part1_value_at_x, ISN0_part2_value_at_x);
	VecWAXPY(alphaN1_at_x, 1.0, ISN1_part1_value_at_x, ISN1_part2_value_at_x);
	VecWAXPY(alphaN2_at_x, 1.0, ISN2_part1_value_at_x, ISN2_part2_value_at_x);

	VecShift(alphaN0_at_x, 1e-10);
	VecShift(alphaN1_at_x, 1e-10);
	VecShift(alphaN2_at_x, 1e-10);

	VecReciprocal(alphaN0_at_x);
	VecReciprocal(alphaN1_at_x);
	VecReciprocal(alphaN2_at_x);

	VecPointwiseMult(alphaN0_at_x, alphaN0_at_x, alphaN0_at_x);
	VecPointwiseMult(alphaN1_at_x, alphaN1_at_x, alphaN1_at_x);
	VecPointwiseMult(alphaN2_at_x, alphaN2_at_x, alphaN2_at_x);

	VecScale(alphaN0_at_x, gammaN0);
	VecScale(alphaN1_at_x, gammaN1);
	VecScale(alphaN2_at_x, gammaN2);

	//求Omega
	VecWAXPY(sum_alphaN_at_x, 1.0, alphaN0_at_x, alphaN1_at_x);
	VecAXPY(sum_alphaN_at_x, 1.0, alphaN2_at_x);

	VecShift(sum_alphaN_at_x, 1e-10);
	VecReciprocal(sum_alphaN_at_x);

	VecPointwiseMult(OmegaN0_at_x, alphaN0_at_x, sum_alphaN_at_x);
	VecPointwiseMult(OmegaN1_at_x, alphaN1_at_x, sum_alphaN_at_x);
	VecPointwiseMult(OmegaN2_at_x, alphaN2_at_x, sum_alphaN_at_x);
}

//WENO at y====================================================================

//一阶导数---------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_WENO_d_variable_y_value(
	Vec variable_f_P, Vec variable_f_N, Vec d_variable_y)
{
	//权重
	get_WeightFP_at_y(variable_f_P);
	get_WeightFN_at_y(variable_f_N);

	//分量
	get_CellBoundaryFP_y(variable_f_P);
	get_CellBoundaryFN_y(variable_f_N);

	//
	VecPointwiseMult(FP_at_y0, OmegaP0_at_y, f_P0_at_y);	VecPointwiseMult(FN_at_y0, OmegaN0_at_y, f_N0_at_y);
	VecPointwiseMult(FP_at_y1, OmegaP1_at_y, f_P1_at_y);	VecPointwiseMult(FN_at_y1, OmegaN1_at_y, f_N1_at_y);
	VecPointwiseMult(FP_at_y2, OmegaP2_at_y, f_P2_at_y);	VecPointwiseMult(FN_at_y2, OmegaN2_at_y, f_N2_at_y);

	VecWAXPY(FP_at_y, 1.0, FP_at_y0, FP_at_y1);
	VecAXPY(FP_at_y, 1.0, FP_at_y2);
	VecWAXPY(FN_at_y, 1.0, FN_at_y0, FN_at_y1);
	VecAXPY(FN_at_y, 1.0, FN_at_y2);

	VecWAXPY(d_variable_y_temp, 1.0, FP_at_y, FN_at_y);

	MatMult(dy_matrix001000_or_000100_to_000000, d_variable_y_temp, d_variable_y);
}

//分量-------------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_CellBoundaryFP_y(Vec variable_f_P)
{
	MatMult(f_P0_y_matrix333333_to_001000, variable_f_P, f_P0_at_y);
	MatMult(f_P1_y_matrix333333_to_001000, variable_f_P, f_P1_at_y);
	MatMult(f_P2_y_matrix333333_to_001000, variable_f_P, f_P2_at_y);
}

void wyh_solve_fluid_structure_interaction_system3D200630::get_CellBoundaryFN_y(Vec variable_f_N)
{
	MatMult(f_N0_y_matrix333333_to_001000, variable_f_N, f_N0_at_y);
	MatMult(f_N1_y_matrix333333_to_001000, variable_f_N, f_N1_at_y);
	MatMult(f_N2_y_matrix333333_to_001000, variable_f_N, f_N2_at_y);
}

//权重-------------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_WeightFP_at_y(Vec variable_f_P)
{
	//part 1
	MatMult(ISP0_y_part1_matrix333333_to_001000, variable_f_P, ISP0_part1_value_at_y);
	MatMult(ISP1_y_part1_matrix333333_to_001000, variable_f_P, ISP1_part1_value_at_y);
	MatMult(ISP2_y_part1_matrix333333_to_001000, variable_f_P, ISP2_part1_value_at_y);

	VecPointwiseMult(ISP0_part1_value_at_y, ISP0_part1_value_at_y, ISP0_part1_value_at_y);
	VecPointwiseMult(ISP1_part1_value_at_y, ISP1_part1_value_at_y, ISP1_part1_value_at_y);
	VecPointwiseMult(ISP2_part1_value_at_y, ISP2_part1_value_at_y, ISP2_part1_value_at_y);

	VecScale(ISP0_part1_value_at_y, 13.0 / 12.0);
	VecScale(ISP1_part1_value_at_y, 13.0 / 12.0);
	VecScale(ISP2_part1_value_at_y, 13.0 / 12.0);

	//part 2
	MatMult(ISP0_y_part2_matrix333333_to_001000, variable_f_P, ISP0_part2_value_at_y);
	MatMult(ISP1_y_part2_matrix333333_to_001000, variable_f_P, ISP1_part2_value_at_y);
	MatMult(ISP2_y_part2_matrix333333_to_001000, variable_f_P, ISP2_part2_value_at_y);

	VecPointwiseMult(ISP0_part2_value_at_y, ISP0_part2_value_at_y, ISP0_part2_value_at_y);
	VecPointwiseMult(ISP1_part2_value_at_y, ISP1_part2_value_at_y, ISP1_part2_value_at_y);
	VecPointwiseMult(ISP2_part2_value_at_y, ISP2_part2_value_at_y, ISP2_part2_value_at_y);

	VecScale(ISP0_part2_value_at_y, 1.0 / 4.0);
	VecScale(ISP1_part2_value_at_y, 1.0 / 4.0);
	VecScale(ISP2_part2_value_at_y, 1.0 / 4.0);

	//求alpha
	double gammaP0 = 0.1;
	double gammaP1 = 0.6;
	double gammaP2 = 0.3;

	VecWAXPY(alphaP0_at_y, 1.0, ISP0_part1_value_at_y, ISP0_part2_value_at_y);
	VecWAXPY(alphaP1_at_y, 1.0, ISP1_part1_value_at_y, ISP1_part2_value_at_y);
	VecWAXPY(alphaP2_at_y, 1.0, ISP2_part1_value_at_y, ISP2_part2_value_at_y);

	VecShift(alphaP0_at_y, 1e-10);
	VecShift(alphaP1_at_y, 1e-10);
	VecShift(alphaP2_at_y, 1e-10);

	VecReciprocal(alphaP0_at_y);
	VecReciprocal(alphaP1_at_y);
	VecReciprocal(alphaP2_at_y);

	VecPointwiseMult(alphaP0_at_y, alphaP0_at_y, alphaP0_at_y);
	VecPointwiseMult(alphaP1_at_y, alphaP1_at_y, alphaP1_at_y);
	VecPointwiseMult(alphaP2_at_y, alphaP2_at_y, alphaP2_at_y);

	VecScale(alphaP0_at_y, gammaP0);
	VecScale(alphaP1_at_y, gammaP1);
	VecScale(alphaP2_at_y, gammaP2);

	//求Omega
	VecWAXPY(sum_alphaP_at_y, 1.0, alphaP0_at_y, alphaP1_at_y);
	VecAXPY(sum_alphaP_at_y, 1.0, alphaP2_at_y);

	VecShift(sum_alphaP_at_y, 1e-10);
	VecReciprocal(sum_alphaP_at_y);

	VecPointwiseMult(OmegaP0_at_y, alphaP0_at_y, sum_alphaP_at_y);
	VecPointwiseMult(OmegaP1_at_y, alphaP1_at_y, sum_alphaP_at_y);
	VecPointwiseMult(OmegaP2_at_y, alphaP2_at_y, sum_alphaP_at_y);
}

void wyh_solve_fluid_structure_interaction_system3D200630::get_WeightFN_at_y(Vec variable_f_N)
{
	//part 1
	MatMult(ISN0_y_part1_matrix333333_to_001000, variable_f_N, ISN0_part1_value_at_y);
	MatMult(ISN1_y_part1_matrix333333_to_001000, variable_f_N, ISN1_part1_value_at_y);
	MatMult(ISN2_y_part1_matrix333333_to_001000, variable_f_N, ISN2_part1_value_at_y);

	VecPointwiseMult(ISN0_part1_value_at_y, ISN0_part1_value_at_y, ISN0_part1_value_at_y);
	VecPointwiseMult(ISN1_part1_value_at_y, ISN1_part1_value_at_y, ISN1_part1_value_at_y);
	VecPointwiseMult(ISN2_part1_value_at_y, ISN2_part1_value_at_y, ISN2_part1_value_at_y);

	VecScale(ISN0_part1_value_at_y, 13.0 / 12.0);
	VecScale(ISN1_part1_value_at_y, 13.0 / 12.0);
	VecScale(ISN2_part1_value_at_y, 13.0 / 12.0);

	//part 2
	MatMult(ISN0_y_part2_matrix333333_to_001000, variable_f_N, ISN0_part2_value_at_y);
	MatMult(ISN1_y_part2_matrix333333_to_001000, variable_f_N, ISN1_part2_value_at_y);
	MatMult(ISN2_y_part2_matrix333333_to_001000, variable_f_N, ISN2_part2_value_at_y);

	VecPointwiseMult(ISN0_part2_value_at_y, ISN0_part2_value_at_y, ISN0_part2_value_at_y);
	VecPointwiseMult(ISN1_part2_value_at_y, ISN1_part2_value_at_y, ISN1_part2_value_at_y);
	VecPointwiseMult(ISN2_part2_value_at_y, ISN2_part2_value_at_y, ISN2_part2_value_at_y);

	VecScale(ISN0_part2_value_at_y, 1.0 / 4.0);
	VecScale(ISN1_part2_value_at_y, 1.0 / 4.0);
	VecScale(ISN2_part2_value_at_y, 1.0 / 4.0);

	//求alpha
	double gammaN0 = 0.3;
	double gammaN1 = 0.6;
	double gammaN2 = 0.1;

	VecWAXPY(alphaN0_at_y, 1.0, ISN0_part1_value_at_y, ISN0_part2_value_at_y);
	VecWAXPY(alphaN1_at_y, 1.0, ISN1_part1_value_at_y, ISN1_part2_value_at_y);
	VecWAXPY(alphaN2_at_y, 1.0, ISN2_part1_value_at_y, ISN2_part2_value_at_y);

	VecShift(alphaN0_at_y, 1e-10);
	VecShift(alphaN1_at_y, 1e-10);
	VecShift(alphaN2_at_y, 1e-10);

	VecReciprocal(alphaN0_at_y);
	VecReciprocal(alphaN1_at_y);
	VecReciprocal(alphaN2_at_y);

	VecPointwiseMult(alphaN0_at_y, alphaN0_at_y, alphaN0_at_y);
	VecPointwiseMult(alphaN1_at_y, alphaN1_at_y, alphaN1_at_y);
	VecPointwiseMult(alphaN2_at_y, alphaN2_at_y, alphaN2_at_y);

	VecScale(alphaN0_at_y, gammaN0);
	VecScale(alphaN1_at_y, gammaN1);
	VecScale(alphaN2_at_y, gammaN2);

	//求Omega
	VecWAXPY(sum_alphaN_at_y, 1.0, alphaN0_at_y, alphaN1_at_y);
	VecAXPY(sum_alphaN_at_y, 1.0, alphaN2_at_y);

	VecShift(sum_alphaN_at_y, 1e-10);
	VecReciprocal(sum_alphaN_at_y);

	VecPointwiseMult(OmegaN0_at_y, alphaN0_at_y, sum_alphaN_at_y);
	VecPointwiseMult(OmegaN1_at_y, alphaN1_at_y, sum_alphaN_at_y);
	VecPointwiseMult(OmegaN2_at_y, alphaN2_at_y, sum_alphaN_at_y);
}

//WENO at z====================================================================

//一阶导数---------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_WENO_d_variable_z_value(
	Vec variable_f_P, Vec variable_f_N, Vec d_variable_z)
{
	//权重
	get_WeightFP_at_z(variable_f_P);
	get_WeightFN_at_z(variable_f_N);

	//分量
	get_CellBoundaryFP_z(variable_f_P);
	get_CellBoundaryFN_z(variable_f_N);

	//
	VecPointwiseMult(FP_at_z0, OmegaP0_at_z, f_P0_at_z);	VecPointwiseMult(FN_at_z0, OmegaN0_at_z, f_N0_at_z);
	VecPointwiseMult(FP_at_z1, OmegaP1_at_z, f_P1_at_z);	VecPointwiseMult(FN_at_z1, OmegaN1_at_z, f_N1_at_z);
	VecPointwiseMult(FP_at_z2, OmegaP2_at_z, f_P2_at_z);	VecPointwiseMult(FN_at_z2, OmegaN2_at_z, f_N2_at_z);

	VecWAXPY(FP_at_z, 1.0, FP_at_z0, FP_at_z1);
	VecAXPY(FP_at_z, 1.0, FP_at_z2);
	VecWAXPY(FN_at_z, 1.0, FN_at_z0, FN_at_z1);
	VecAXPY(FN_at_z, 1.0, FN_at_z2);

	VecWAXPY(d_variable_z_temp, 1.0, FP_at_z, FN_at_z);

	MatMult(dz_matrix000010_or_000001_to_000000, d_variable_z_temp, d_variable_z);
}

//分量-------------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_CellBoundaryFP_z(Vec variable_f_P)
{
	MatMult(f_P0_z_matrix333333_to_000010, variable_f_P, f_P0_at_z);
	MatMult(f_P1_z_matrix333333_to_000010, variable_f_P, f_P1_at_z);
	MatMult(f_P2_z_matrix333333_to_000010, variable_f_P, f_P2_at_z);
}

void wyh_solve_fluid_structure_interaction_system3D200630::get_CellBoundaryFN_z(Vec variable_f_N)
{
	MatMult(f_N0_z_matrix333333_to_000010, variable_f_N, f_N0_at_z);
	MatMult(f_N1_z_matrix333333_to_000010, variable_f_N, f_N1_at_z);
	MatMult(f_N2_z_matrix333333_to_000010, variable_f_N, f_N2_at_z);
}

//权重-------------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::get_WeightFP_at_z(Vec variable_f_P)
{
	//part 1
	MatMult(ISP0_z_part1_matrix333333_to_000010, variable_f_P, ISP0_part1_value_at_z);
	MatMult(ISP1_z_part1_matrix333333_to_000010, variable_f_P, ISP1_part1_value_at_z);
	MatMult(ISP2_z_part1_matrix333333_to_000010, variable_f_P, ISP2_part1_value_at_z);

	VecPointwiseMult(ISP0_part1_value_at_z, ISP0_part1_value_at_z, ISP0_part1_value_at_z);
	VecPointwiseMult(ISP1_part1_value_at_z, ISP1_part1_value_at_z, ISP1_part1_value_at_z);
	VecPointwiseMult(ISP2_part1_value_at_z, ISP2_part1_value_at_z, ISP2_part1_value_at_z);

	VecScale(ISP0_part1_value_at_z, 13.0 / 12.0);
	VecScale(ISP1_part1_value_at_z, 13.0 / 12.0);
	VecScale(ISP2_part1_value_at_z, 13.0 / 12.0);

	//part 2
	MatMult(ISP0_z_part2_matrix333333_to_000010, variable_f_P, ISP0_part2_value_at_z);
	MatMult(ISP1_z_part2_matrix333333_to_000010, variable_f_P, ISP1_part2_value_at_z);
	MatMult(ISP2_z_part2_matrix333333_to_000010, variable_f_P, ISP2_part2_value_at_z);

	VecPointwiseMult(ISP0_part2_value_at_z, ISP0_part2_value_at_z, ISP0_part2_value_at_z);
	VecPointwiseMult(ISP1_part2_value_at_z, ISP1_part2_value_at_z, ISP1_part2_value_at_z);
	VecPointwiseMult(ISP2_part2_value_at_z, ISP2_part2_value_at_z, ISP2_part2_value_at_z);

	VecScale(ISP0_part2_value_at_z, 1.0 / 4.0);
	VecScale(ISP1_part2_value_at_z, 1.0 / 4.0);
	VecScale(ISP2_part2_value_at_z, 1.0 / 4.0);

	//求alpha
	double gammaP0 = 0.1;
	double gammaP1 = 0.6;
	double gammaP2 = 0.3;

	VecWAXPY(alphaP0_at_z, 1.0, ISP0_part1_value_at_z, ISP0_part2_value_at_z);
	VecWAXPY(alphaP1_at_z, 1.0, ISP1_part1_value_at_z, ISP1_part2_value_at_z);
	VecWAXPY(alphaP2_at_z, 1.0, ISP2_part1_value_at_z, ISP2_part2_value_at_z);

	VecShift(alphaP0_at_z, 1e-10);
	VecShift(alphaP1_at_z, 1e-10);
	VecShift(alphaP2_at_z, 1e-10);

	VecReciprocal(alphaP0_at_z);
	VecReciprocal(alphaP1_at_z);
	VecReciprocal(alphaP2_at_z);

	VecPointwiseMult(alphaP0_at_z, alphaP0_at_z, alphaP0_at_z);
	VecPointwiseMult(alphaP1_at_z, alphaP1_at_z, alphaP1_at_z);
	VecPointwiseMult(alphaP2_at_z, alphaP2_at_z, alphaP2_at_z);

	VecScale(alphaP0_at_z, gammaP0);
	VecScale(alphaP1_at_z, gammaP1);
	VecScale(alphaP2_at_z, gammaP2);

	//求Omega
	VecWAXPY(sum_alphaP_at_z, 1.0, alphaP0_at_z, alphaP1_at_z);
	VecAXPY(sum_alphaP_at_z, 1.0, alphaP2_at_z);

	VecShift(sum_alphaP_at_z, 1e-10);
	VecReciprocal(sum_alphaP_at_z);

	VecPointwiseMult(OmegaP0_at_z, alphaP0_at_z, sum_alphaP_at_z);
	VecPointwiseMult(OmegaP1_at_z, alphaP1_at_z, sum_alphaP_at_z);
	VecPointwiseMult(OmegaP2_at_z, alphaP2_at_z, sum_alphaP_at_z);
}

void wyh_solve_fluid_structure_interaction_system3D200630::get_WeightFN_at_z(Vec variable_f_N)
{
	//part 1
	MatMult(ISN0_z_part1_matrix333333_to_000010, variable_f_N, ISN0_part1_value_at_z);
	MatMult(ISN1_z_part1_matrix333333_to_000010, variable_f_N, ISN1_part1_value_at_z);
	MatMult(ISN2_z_part1_matrix333333_to_000010, variable_f_N, ISN2_part1_value_at_z);

	VecPointwiseMult(ISN0_part1_value_at_z, ISN0_part1_value_at_z, ISN0_part1_value_at_z);
	VecPointwiseMult(ISN1_part1_value_at_z, ISN1_part1_value_at_z, ISN1_part1_value_at_z);
	VecPointwiseMult(ISN2_part1_value_at_z, ISN2_part1_value_at_z, ISN2_part1_value_at_z);

	VecScale(ISN0_part1_value_at_z, 13.0 / 12.0);
	VecScale(ISN1_part1_value_at_z, 13.0 / 12.0);
	VecScale(ISN2_part1_value_at_z, 13.0 / 12.0);

	//part 2
	MatMult(ISN0_z_part2_matrix333333_to_000010, variable_f_N, ISN0_part2_value_at_z);
	MatMult(ISN1_z_part2_matrix333333_to_000010, variable_f_N, ISN1_part2_value_at_z);
	MatMult(ISN2_z_part2_matrix333333_to_000010, variable_f_N, ISN2_part2_value_at_z);

	VecPointwiseMult(ISN0_part2_value_at_z, ISN0_part2_value_at_z, ISN0_part2_value_at_z);
	VecPointwiseMult(ISN1_part2_value_at_z, ISN1_part2_value_at_z, ISN1_part2_value_at_z);
	VecPointwiseMult(ISN2_part2_value_at_z, ISN2_part2_value_at_z, ISN2_part2_value_at_z);

	VecScale(ISN0_part2_value_at_z, 1.0 / 4.0);
	VecScale(ISN1_part2_value_at_z, 1.0 / 4.0);
	VecScale(ISN2_part2_value_at_z, 1.0 / 4.0);

	//求alpha
	double gammaN0 = 0.3;
	double gammaN1 = 0.6;
	double gammaN2 = 0.1;

	VecWAXPY(alphaN0_at_z, 1.0, ISN0_part1_value_at_z, ISN0_part2_value_at_z);
	VecWAXPY(alphaN1_at_z, 1.0, ISN1_part1_value_at_z, ISN1_part2_value_at_z);
	VecWAXPY(alphaN2_at_z, 1.0, ISN2_part1_value_at_z, ISN2_part2_value_at_z);

	VecShift(alphaN0_at_z, 1e-10);
	VecShift(alphaN1_at_z, 1e-10);
	VecShift(alphaN2_at_z, 1e-10);

	VecReciprocal(alphaN0_at_z);
	VecReciprocal(alphaN1_at_z);
	VecReciprocal(alphaN2_at_z);

	VecPointwiseMult(alphaN0_at_z, alphaN0_at_z, alphaN0_at_z);
	VecPointwiseMult(alphaN1_at_z, alphaN1_at_z, alphaN1_at_z);
	VecPointwiseMult(alphaN2_at_z, alphaN2_at_z, alphaN2_at_z);

	VecScale(alphaN0_at_z, gammaN0);
	VecScale(alphaN1_at_z, gammaN1);
	VecScale(alphaN2_at_z, gammaN2);

	//求Omega
	VecWAXPY(sum_alphaN_at_z, 1.0, alphaN0_at_z, alphaN1_at_z);
	VecAXPY(sum_alphaN_at_z, 1.0, alphaN2_at_z);

	VecShift(sum_alphaN_at_z, 1e-10);
	VecReciprocal(sum_alphaN_at_z);

	VecPointwiseMult(OmegaN0_at_z, alphaN0_at_z, sum_alphaN_at_z);
	VecPointwiseMult(OmegaN1_at_z, alphaN1_at_z, sum_alphaN_at_z);
	VecPointwiseMult(OmegaN2_at_z, alphaN2_at_z, sum_alphaN_at_z);
}

//***************************************************************************************************************************************************

//WENO相关计算矩阵====================================================================================================================================

//***************************************************************************************************************************************************

//分量=============================================================================================

//分量-------------------------------------------

//Cell boundary FP
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_CellBoundaryFP_x_matrix333333_to_100000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(f_P0_x_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P0_x_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(f_P1_x_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P1_x_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(f_P2_x_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P2_x_matrix333333_to_100000, 3, NULL);
	MatZeroEntries(f_P0_x_matrix333333_to_100000);
	MatZeroEntries(f_P1_x_matrix333333_to_100000);
	MatZeroEntries(f_P2_x_matrix333333_to_100000);

	MatGetOwnershipRange(f_P0_x_matrix333333_to_100000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 5);
		PetscInt j = Ii02 - i * (ghost_m_num - 5);

		int Ij_0 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 2);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 4);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0 / 3.0, -7.0 / 6.0, 11.0 / 6.0 };
		double insert_value1_vector[3] = { -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0 };
		double insert_value2_vector[3] = { 1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0 };

		MatSetValues(f_P0_x_matrix333333_to_100000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(f_P1_x_matrix333333_to_100000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(f_P2_x_matrix333333_to_100000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(f_P0_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P0_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_P1_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P1_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_P2_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P2_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_CellBoundaryFP_y_matrix333333_to_001000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(f_P0_y_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P0_y_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(f_P1_y_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P1_y_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(f_P2_y_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P2_y_matrix333333_to_001000, 3, NULL);
	MatZeroEntries(f_P0_y_matrix333333_to_001000);
	MatZeroEntries(f_P1_y_matrix333333_to_001000);
	MatZeroEntries(f_P2_y_matrix333333_to_001000);

	MatGetOwnershipRange(f_P0_y_matrix333333_to_001000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_0 = (k + 3)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 3);
		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 2)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 4)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0 / 3.0, -7.0 / 6.0, 11.0 / 6.0 };
		double insert_value1_vector[3] = { -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0 };
		double insert_value2_vector[3] = { 1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0 };

		MatSetValues(f_P0_y_matrix333333_to_001000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(f_P1_y_matrix333333_to_001000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(f_P2_y_matrix333333_to_001000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(f_P0_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P0_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_P1_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P1_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_P2_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P2_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_CellBoundaryFP_z_matrix333333_to_000010()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(f_P0_z_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P0_z_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(f_P1_z_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P1_z_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(f_P2_z_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_P2_z_matrix333333_to_000010, 3, NULL);
	MatZeroEntries(f_P0_z_matrix333333_to_000010);
	MatZeroEntries(f_P1_z_matrix333333_to_000010);
	MatZeroEntries(f_P2_z_matrix333333_to_000010);

	MatGetOwnershipRange(f_P0_z_matrix333333_to_000010, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 2)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 4)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0 / 3.0, -7.0 / 6.0, 11.0 / 6.0 };
		double insert_value1_vector[3] = { -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0 };
		double insert_value2_vector[3] = { 1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0 };

		MatSetValues(f_P0_z_matrix333333_to_000010, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(f_P1_z_matrix333333_to_000010, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(f_P2_z_matrix333333_to_000010, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(f_P0_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P0_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_P1_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P1_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_P2_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_P2_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
}

//Cell boundary FN
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_CellBoundaryFN_x_matrix333333_to_100000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(f_N0_x_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N0_x_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(f_N1_x_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N1_x_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(f_N2_x_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N2_x_matrix333333_to_100000, 3, NULL);
	MatZeroEntries(f_N0_x_matrix333333_to_100000);
	MatZeroEntries(f_N1_x_matrix333333_to_100000);
	MatZeroEntries(f_N2_x_matrix333333_to_100000);

	MatGetOwnershipRange(f_N0_x_matrix333333_to_100000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 5);
		PetscInt j = Ii02 - i * (ghost_m_num - 5);

		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 2);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 4);
		int Ij_5 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 5);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0 };
		double insert_value1_vector[3] = { 1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0 };
		double insert_value2_vector[3] = { 11.0 / 6.0, -7.0 / 6.0, 1.0 / 3.0 };

		MatSetValues(f_N0_x_matrix333333_to_100000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(f_N1_x_matrix333333_to_100000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(f_N2_x_matrix333333_to_100000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(f_N0_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N0_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_N1_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N1_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_N2_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N2_x_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_CellBoundaryFN_y_matrix333333_to_001000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(f_N0_y_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N0_y_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(f_N1_y_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N1_y_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(f_N2_y_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N2_y_matrix333333_to_001000, 3, NULL);
	MatZeroEntries(f_N0_y_matrix333333_to_001000);
	MatZeroEntries(f_N1_y_matrix333333_to_001000);
	MatZeroEntries(f_N2_y_matrix333333_to_001000);

	MatGetOwnershipRange(f_N0_y_matrix333333_to_001000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 2)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 4)*ghost_m_num + (j + 3);
		int Ij_5 = (k + 3)*ghost_m_num*ghost_n_num + (i + 5)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0 };
		double insert_value1_vector[3] = { 1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0 };
		double insert_value2_vector[3] = { 11.0 / 6.0, -7.0 / 6.0, 1.0 / 3.0 };

		MatSetValues(f_N0_y_matrix333333_to_001000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(f_N1_y_matrix333333_to_001000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(f_N2_y_matrix333333_to_001000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(f_N0_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N0_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_N1_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N1_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_N2_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N2_y_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_CellBoundaryFN_z_matrix333333_to_000010()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(f_N0_z_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N0_z_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(f_N1_z_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N1_z_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(f_N2_z_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(f_N2_z_matrix333333_to_000010, 3, NULL);
	MatZeroEntries(f_N0_z_matrix333333_to_000010);
	MatZeroEntries(f_N1_z_matrix333333_to_000010);
	MatZeroEntries(f_N2_z_matrix333333_to_000010);

	MatGetOwnershipRange(f_N0_z_matrix333333_to_000010, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 2)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 4)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_5 = (k + 5)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0 };
		double insert_value1_vector[3] = { 1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0 };
		double insert_value2_vector[3] = { 11.0 / 6.0, -7.0 / 6.0, 1.0 / 3.0 };

		MatSetValues(f_N0_z_matrix333333_to_000010, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(f_N1_z_matrix333333_to_000010, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(f_N2_z_matrix333333_to_000010, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(f_N0_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N0_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_N1_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N1_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(f_N2_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(f_N2_z_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
}

//权重-------------------------------------------

//part 1 of the Weight FP
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part1_WeightFP_x_matrix333333_to_100000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISP0_x_part1_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP0_x_part1_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(ISP1_x_part1_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP1_x_part1_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(ISP2_x_part1_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP2_x_part1_matrix333333_to_100000, 3, NULL);
	MatZeroEntries(ISP0_x_part1_matrix333333_to_100000);
	MatZeroEntries(ISP1_x_part1_matrix333333_to_100000);
	MatZeroEntries(ISP2_x_part1_matrix333333_to_100000);

	MatGetOwnershipRange(ISP0_x_part1_matrix333333_to_100000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 5);
		PetscInt j = Ii02 - i * (ghost_m_num - 5);

		int Ij_0 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 2);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 4);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value1_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value2_vector[3] = { 1.0, -2.0, 1.0 };

		MatSetValues(ISP0_x_part1_matrix333333_to_100000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISP1_x_part1_matrix333333_to_100000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISP2_x_part1_matrix333333_to_100000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISP0_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP0_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP1_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP1_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP2_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP2_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part1_WeightFP_y_matrix333333_to_001000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISP0_y_part1_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP0_y_part1_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(ISP1_y_part1_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP1_y_part1_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(ISP2_y_part1_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP2_y_part1_matrix333333_to_001000, 3, NULL);
	MatZeroEntries(ISP0_y_part1_matrix333333_to_001000);
	MatZeroEntries(ISP1_y_part1_matrix333333_to_001000);
	MatZeroEntries(ISP2_y_part1_matrix333333_to_001000);

	MatGetOwnershipRange(ISP0_y_part1_matrix333333_to_001000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_0 = (k + 3)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 3);
		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 2)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 4)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value1_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value2_vector[3] = { 1.0, -2.0, 1.0 };

		MatSetValues(ISP0_y_part1_matrix333333_to_001000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISP1_y_part1_matrix333333_to_001000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISP2_y_part1_matrix333333_to_001000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISP0_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP0_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP1_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP1_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP2_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP2_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part1_WeightFP_z_matrix333333_to_000010()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISP0_z_part1_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP0_z_part1_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(ISP1_z_part1_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP1_z_part1_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(ISP2_z_part1_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP2_z_part1_matrix333333_to_000010, 3, NULL);
	MatZeroEntries(ISP0_z_part1_matrix333333_to_000010);
	MatZeroEntries(ISP1_z_part1_matrix333333_to_000010);
	MatZeroEntries(ISP2_z_part1_matrix333333_to_000010);

	MatGetOwnershipRange(ISP0_z_part1_matrix333333_to_000010, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 2)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 4)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value1_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value2_vector[3] = { 1.0, -2.0, 1.0 };

		MatSetValues(ISP0_z_part1_matrix333333_to_000010, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISP1_z_part1_matrix333333_to_000010, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISP2_z_part1_matrix333333_to_000010, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISP0_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP0_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP1_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP1_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP2_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP2_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
}

//part 2 of the Weight FP
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part2_WeightFP_x_matrix333333_to_100000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISP0_x_part2_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP0_x_part2_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(ISP1_x_part2_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP1_x_part2_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(ISP2_x_part2_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP2_x_part2_matrix333333_to_100000, 3, NULL);
	MatZeroEntries(ISP0_x_part2_matrix333333_to_100000);
	MatZeroEntries(ISP1_x_part2_matrix333333_to_100000);
	MatZeroEntries(ISP2_x_part2_matrix333333_to_100000);

	MatGetOwnershipRange(ISP0_x_part2_matrix333333_to_100000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 5);
		PetscInt j = Ii02 - i * (ghost_m_num - 5);

		int Ij_0 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 2);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 4);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0, -4.0, 3.0 };
		double insert_value1_vector[3] = { 1.0, 0.0, -1.0 };
		double insert_value2_vector[3] = { 3.0, -4.0, 1.0 };

		MatSetValues(ISP0_x_part2_matrix333333_to_100000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISP1_x_part2_matrix333333_to_100000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISP2_x_part2_matrix333333_to_100000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISP0_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP0_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP1_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP1_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP2_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP2_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part2_WeightFP_y_matrix333333_to_001000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISP0_y_part2_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP0_y_part2_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(ISP1_y_part2_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP1_y_part2_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(ISP2_y_part2_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP2_y_part2_matrix333333_to_001000, 3, NULL);
	MatZeroEntries(ISP0_y_part2_matrix333333_to_001000);
	MatZeroEntries(ISP1_y_part2_matrix333333_to_001000);
	MatZeroEntries(ISP2_y_part2_matrix333333_to_001000);

	MatGetOwnershipRange(ISP0_y_part2_matrix333333_to_001000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_0 = (k + 3)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 3);
		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 2)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 4)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0, -4.0, 3.0 };
		double insert_value1_vector[3] = { 1.0, 0.0, -1.0 };
		double insert_value2_vector[3] = { 3.0, -4.0, 1.0 };

		MatSetValues(ISP0_y_part2_matrix333333_to_001000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISP1_y_part2_matrix333333_to_001000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISP2_y_part2_matrix333333_to_001000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISP0_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP0_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP1_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP1_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP2_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP2_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part2_WeightFP_z_matrix333333_to_000010()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISP0_z_part2_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP0_z_part2_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(ISP1_z_part2_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP1_z_part2_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(ISP2_z_part2_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISP2_z_part2_matrix333333_to_000010, 3, NULL);
	MatZeroEntries(ISP0_z_part2_matrix333333_to_000010);
	MatZeroEntries(ISP1_z_part2_matrix333333_to_000010);
	MatZeroEntries(ISP2_z_part2_matrix333333_to_000010);

	MatGetOwnershipRange(ISP0_z_part2_matrix333333_to_000010, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 2)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 4)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_0, Ij_1, Ij_2 };
		int Ij_vector1[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector2[3] = { Ij_2, Ij_3, Ij_4 };

		double insert_value0_vector[3] = { 1.0, -4.0, 3.0 };
		double insert_value1_vector[3] = { 1.0, 0.0, -1.0 };
		double insert_value2_vector[3] = { 3.0, -4.0, 1.0 };

		MatSetValues(ISP0_z_part2_matrix333333_to_000010, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISP1_z_part2_matrix333333_to_000010, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISP2_z_part2_matrix333333_to_000010, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISP0_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP0_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP1_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP1_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISP2_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISP2_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
}

//part 1 of the Weight FN
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part1_WeightFN_x_matrix333333_to_100000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISN0_x_part1_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN0_x_part1_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(ISN1_x_part1_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN1_x_part1_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(ISN2_x_part1_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN2_x_part1_matrix333333_to_100000, 3, NULL);
	MatZeroEntries(ISN0_x_part1_matrix333333_to_100000);
	MatZeroEntries(ISN1_x_part1_matrix333333_to_100000);
	MatZeroEntries(ISN2_x_part1_matrix333333_to_100000);

	MatGetOwnershipRange(ISN0_x_part1_matrix333333_to_100000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 5);
		PetscInt j = Ii02 - i * (ghost_m_num - 5);

		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 2);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 4);
		int Ij_5 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 5);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value1_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value2_vector[3] = { 1.0, -2.0, 1.0 };

		MatSetValues(ISN0_x_part1_matrix333333_to_100000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISN1_x_part1_matrix333333_to_100000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISN2_x_part1_matrix333333_to_100000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISN0_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN0_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN1_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN1_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN2_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN2_x_part1_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part1_WeightFN_y_matrix333333_to_001000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISN0_y_part1_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN0_y_part1_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(ISN1_y_part1_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN1_y_part1_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(ISN2_y_part1_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN2_y_part1_matrix333333_to_001000, 3, NULL);
	MatZeroEntries(ISN0_y_part1_matrix333333_to_001000);
	MatZeroEntries(ISN1_y_part1_matrix333333_to_001000);
	MatZeroEntries(ISN2_y_part1_matrix333333_to_001000);

	MatGetOwnershipRange(ISN0_y_part1_matrix333333_to_001000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 2)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 4)*ghost_m_num + (j + 3);
		int Ij_5 = (k + 3)*ghost_m_num*ghost_n_num + (i + 5)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value1_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value2_vector[3] = { 1.0, -2.0, 1.0 };

		MatSetValues(ISN0_y_part1_matrix333333_to_001000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISN1_y_part1_matrix333333_to_001000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISN2_y_part1_matrix333333_to_001000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISN0_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN0_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN1_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN1_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN2_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN2_y_part1_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part1_WeightFN_z_matrix333333_to_000010()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISN0_z_part1_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN0_z_part1_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(ISN1_z_part1_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN1_z_part1_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(ISN2_z_part1_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN2_z_part1_matrix333333_to_000010, 3, NULL);
	MatZeroEntries(ISN0_z_part1_matrix333333_to_000010);
	MatZeroEntries(ISN1_z_part1_matrix333333_to_000010);
	MatZeroEntries(ISN2_z_part1_matrix333333_to_000010);

	MatGetOwnershipRange(ISN0_z_part1_matrix333333_to_000010, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 2)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 4)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_5 = (k + 5)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value1_vector[3] = { 1.0, -2.0, 1.0 };
		double insert_value2_vector[3] = { 1.0, -2.0, 1.0 };

		MatSetValues(ISN0_z_part1_matrix333333_to_000010, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISN1_z_part1_matrix333333_to_000010, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISN2_z_part1_matrix333333_to_000010, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISN0_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN0_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN1_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN1_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN2_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN2_z_part1_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
}

//part 2 of the Weight FN
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part2_WeightFN_x_matrix333333_to_100000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISN0_x_part2_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN0_x_part2_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(ISN1_x_part2_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN1_x_part2_matrix333333_to_100000, 3, NULL);
	MatMPIAIJSetPreallocation(ISN2_x_part2_matrix333333_to_100000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN2_x_part2_matrix333333_to_100000, 3, NULL);
	MatZeroEntries(ISN0_x_part2_matrix333333_to_100000);
	MatZeroEntries(ISN1_x_part2_matrix333333_to_100000);
	MatZeroEntries(ISN2_x_part2_matrix333333_to_100000);

	MatGetOwnershipRange(ISN0_x_part2_matrix333333_to_100000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 5)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 5);
		PetscInt j = Ii02 - i * (ghost_m_num - 5);

		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 2);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 4);
		int Ij_5 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 5);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { 1.0, -4.0, 3.0 };
		double insert_value1_vector[3] = { 1.0, 0.0, -1.0 };
		double insert_value2_vector[3] = { 3.0, -4.0, 1.0 };

		MatSetValues(ISN0_x_part2_matrix333333_to_100000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISN1_x_part2_matrix333333_to_100000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISN2_x_part2_matrix333333_to_100000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISN0_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN0_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN1_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN1_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN2_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN2_x_part2_matrix333333_to_100000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part2_WeightFN_y_matrix333333_to_001000()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISN0_y_part2_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN0_y_part2_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(ISN1_y_part2_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN1_y_part2_matrix333333_to_001000, 3, NULL);
	MatMPIAIJSetPreallocation(ISN2_y_part2_matrix333333_to_001000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN2_y_part2_matrix333333_to_001000, 3, NULL);
	MatZeroEntries(ISN0_y_part2_matrix333333_to_001000);
	MatZeroEntries(ISN1_y_part2_matrix333333_to_001000);
	MatZeroEntries(ISN2_y_part2_matrix333333_to_001000);

	MatGetOwnershipRange(ISN0_y_part2_matrix333333_to_001000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 5));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_1 = (k + 3)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 3)*ghost_m_num*ghost_n_num + (i + 2)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 3)*ghost_m_num*ghost_n_num + (i + 4)*ghost_m_num + (j + 3);
		int Ij_5 = (k + 3)*ghost_m_num*ghost_n_num + (i + 5)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { 1.0, -4.0, 3.0 };
		double insert_value1_vector[3] = { 1.0, 0.0, -1.0 };
		double insert_value2_vector[3] = { 3.0, -4.0, 1.0 };

		MatSetValues(ISN0_y_part2_matrix333333_to_001000, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISN1_y_part2_matrix333333_to_001000, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISN2_y_part2_matrix333333_to_001000, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISN0_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN0_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN1_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN1_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN2_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN2_y_part2_matrix333333_to_001000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_part2_WeightFN_z_matrix333333_to_000010()
{
	//边界
	int ex_j_L = 3;		int ex_j_H = 3;
	int ex_i_L = 3;		int ex_i_H = 3;
	int ex_k_L = 3;		int ex_k_H = 3;
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	//
	MatMPIAIJSetPreallocation(ISN0_z_part2_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN0_z_part2_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(ISN1_z_part2_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN1_z_part2_matrix333333_to_000010, 3, NULL);
	MatMPIAIJSetPreallocation(ISN2_z_part2_matrix333333_to_000010, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ISN2_z_part2_matrix333333_to_000010, 3, NULL);
	MatZeroEntries(ISN0_z_part2_matrix333333_to_000010);
	MatZeroEntries(ISN1_z_part2_matrix333333_to_000010);
	MatZeroEntries(ISN2_z_part2_matrix333333_to_000010);

	MatGetOwnershipRange(ISN0_z_part2_matrix333333_to_000010, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 6)*(ghost_n_num - 6));
		PetscInt i = Ii02 / (ghost_m_num - 6);
		PetscInt j = Ii02 - i * (ghost_m_num - 6);

		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_2 = (k + 2)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_3 = (k + 3)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_4 = (k + 4)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);
		int Ij_5 = (k + 5)*ghost_m_num*ghost_n_num + (i + 3)*ghost_m_num + (j + 3);

		int Ij_vector0[3] = { Ij_1, Ij_2, Ij_3 };
		int Ij_vector1[3] = { Ij_2, Ij_3, Ij_4 };
		int Ij_vector2[3] = { Ij_3, Ij_4, Ij_5 };

		double insert_value0_vector[3] = { 1.0, -4.0, 3.0 };
		double insert_value1_vector[3] = { 1.0, 0.0, -1.0 };
		double insert_value2_vector[3] = { 3.0, -4.0, 1.0 };

		MatSetValues(ISN0_z_part2_matrix333333_to_000010, 1, &Ii, 3, Ij_vector0, insert_value0_vector, INSERT_VALUES);
		MatSetValues(ISN1_z_part2_matrix333333_to_000010, 1, &Ii, 3, Ij_vector1, insert_value1_vector, INSERT_VALUES);
		MatSetValues(ISN2_z_part2_matrix333333_to_000010, 1, &Ii, 3, Ij_vector2, insert_value2_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ISN0_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN0_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN1_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN1_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(ISN2_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ISN2_z_part2_matrix333333_to_000010, MAT_FINAL_ASSEMBLY);
}

//***************************************************************************************************************************************************

//求导相关矩阵========================================================================================================================================

//***************************************************************************************************************************************************

//1阶差分==========================================================================================
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_dx_matrix100000_or_010000_to_000000()
{
	int ghost_m_num = m_num + 1;
	int ghost_n_num = n_num + 0;
	int ghost_o_num = o_num + 0;

	//
	MatMPIAIJSetPreallocation(dx_matrix100000_or_010000_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(dx_matrix100000_or_010000_to_000000, 2, NULL);
	MatZeroEntries(dx_matrix100000_or_010000_to_000000);

	MatGetOwnershipRange(dx_matrix100000_or_010000_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 1)*(ghost_n_num - 0));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 1)*(ghost_n_num - 0));
		PetscInt i = Ii02 / (ghost_m_num - 1);
		PetscInt j = Ii02 - i * (ghost_m_num - 1);

		int Ij_0 = k*ghost_m_num*ghost_n_num + i*ghost_m_num + (j + 0);
		int Ij_1 = k*ghost_m_num*ghost_n_num + i*ghost_m_num + (j + 1);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { -1.0 / hx, 1.0 / hx };

		MatSetValues(dx_matrix100000_or_010000_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(dx_matrix100000_or_010000_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(dx_matrix100000_or_010000_to_000000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_dy_matrix001000_or_000100_to_000000()
{
	int ghost_m_num = m_num + 0;
	int ghost_n_num = n_num + 1;
	int ghost_o_num = o_num + 0;

	//
	MatMPIAIJSetPreallocation(dy_matrix001000_or_000100_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(dy_matrix001000_or_000100_to_000000, 2, NULL);
	MatZeroEntries(dy_matrix001000_or_000100_to_000000);

	MatGetOwnershipRange(dy_matrix001000_or_000100_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 0)*(ghost_n_num - 1));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 0)*(ghost_n_num - 1));
		PetscInt i = Ii02 / (ghost_m_num - 0);
		PetscInt j = Ii02 - i * (ghost_m_num - 0);

		int Ij_0 = k*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_1 = k*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 0);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { -1.0 / hy, 1.0 / hy };

		MatSetValues(dy_matrix001000_or_000100_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(dy_matrix001000_or_000100_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(dy_matrix001000_or_000100_to_000000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_dz_matrix000010_or_000001_to_000000()
{
	int ghost_m_num = m_num + 0;
	int ghost_n_num = n_num + 0;
	int ghost_o_num = o_num + 1;

	//
	MatMPIAIJSetPreallocation(dz_matrix000010_or_000001_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(dz_matrix000010_or_000001_to_000000, 2, NULL);
	MatZeroEntries(dz_matrix000010_or_000001_to_000000);

	MatGetOwnershipRange(dz_matrix000010_or_000001_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 0)*(ghost_n_num - 0));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 0)*(ghost_n_num - 0));
		PetscInt i = Ii02 / (ghost_m_num - 0);
		PetscInt j = Ii02 - i * (ghost_m_num - 0);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { -1.0 / hz, 1.0 / hz };

		MatSetValues(dz_matrix000010_or_000001_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(dz_matrix000010_or_000001_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(dz_matrix000010_or_000001_to_000000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_dx_matrix111111_to_000000()
{
	int ghost_m_num = m_num + 2;
	int ghost_n_num = n_num + 2;
	int ghost_o_num = o_num + 2;

	//
	MatMPIAIJSetPreallocation(dx_matrix111111_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(dx_matrix111111_to_000000, 2, NULL);
	MatZeroEntries(dx_matrix111111_to_000000);

	MatGetOwnershipRange(dx_matrix111111_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt i = Ii02 / (ghost_m_num - 2);
		PetscInt j = Ii02 - i * (ghost_m_num - 2);

		int Ij_0 = (k + 1)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 2);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { -1.0 / (2 * hx), 1.0 / (2 * hx) };

		MatSetValues(dx_matrix111111_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(dx_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(dx_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_dy_matrix111111_to_000000()
{
	int ghost_m_num = m_num + 2;
	int ghost_n_num = n_num + 2;
	int ghost_o_num = o_num + 2;

	//
	MatMPIAIJSetPreallocation(dy_matrix111111_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(dy_matrix111111_to_000000, 2, NULL);
	MatZeroEntries(dy_matrix111111_to_000000);

	MatGetOwnershipRange(dy_matrix111111_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt i = Ii02 / (ghost_m_num - 2);
		PetscInt j = Ii02 - i * (ghost_m_num - 2);

		int Ij_0 = (k + 1)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 1);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 2)*ghost_m_num + (j + 1);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { -1.0 / (2 * hy), 1.0 / (2 * hy) };

		MatSetValues(dy_matrix111111_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(dy_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(dy_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_dz_matrix111111_to_000000()
{
	int ghost_m_num = m_num + 2;
	int ghost_n_num = n_num + 2;
	int ghost_o_num = o_num + 2;

	//
	MatMPIAIJSetPreallocation(dz_matrix111111_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(dz_matrix111111_to_000000, 2, NULL);
	MatZeroEntries(dz_matrix111111_to_000000);

	MatGetOwnershipRange(dz_matrix111111_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt i = Ii02 / (ghost_m_num - 2);
		PetscInt j = Ii02 - i * (ghost_m_num - 2);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 1);
		int Ij_1 = (k + 2)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 1);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { -1.0 / (2 * hz), 1.0 / (2 * hz) };

		MatSetValues(dz_matrix111111_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(dz_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(dz_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
}

//2阶中心差分------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_ddx_matrix111111_to_000000()
{
	int ghost_m_num = m_num + 2;
	int ghost_n_num = n_num + 2;
	int ghost_o_num = o_num + 2;

	//
	MatMPIAIJSetPreallocation(ddx_matrix111111_to_000000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ddx_matrix111111_to_000000, 3, NULL);
	MatZeroEntries(ddx_matrix111111_to_000000);

	MatGetOwnershipRange(ddx_matrix111111_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt i = Ii02 / (ghost_m_num - 2);
		PetscInt j = Ii02 - i * (ghost_m_num - 2);

		int Ij_0 = (k + 1)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 1)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 2);

		int Ij_vector[3] = { Ij_0, Ij_1, Ij_2 };
		double insert_value_vector[3] = { 1.0 / (hx * hx), -2.0 / (hx * hx), 1.0 / (hx * hx) };

		MatSetValues(ddx_matrix111111_to_000000, 1, &Ii, 3, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ddx_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ddx_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_ddy_matrix111111_to_000000()
{
	int ghost_m_num = m_num + 2;
	int ghost_n_num = n_num + 2;
	int ghost_o_num = o_num + 2;

	//
	MatMPIAIJSetPreallocation(ddy_matrix111111_to_000000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ddy_matrix111111_to_000000, 3, NULL);
	MatZeroEntries(ddy_matrix111111_to_000000);

	MatGetOwnershipRange(ddy_matrix111111_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt i = Ii02 / (ghost_m_num - 2);
		PetscInt j = Ii02 - i * (ghost_m_num - 2);

		int Ij_0 = (k + 1)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 1);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 1)*ghost_m_num*ghost_n_num + (i + 2)*ghost_m_num + (j + 1);

		int Ij_vector[3] = { Ij_0, Ij_1, Ij_2 };
		double insert_value_vector[3] = { 1.0 / (hy * hy), -2.0 / (hy * hy), 1.0 / (hy * hy) };

		MatSetValues(ddy_matrix111111_to_000000, 1, &Ii, 3, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ddy_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ddy_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_ddz_matrix111111_to_000000()
{
	int ghost_m_num = m_num + 2;
	int ghost_n_num = n_num + 2;
	int ghost_o_num = o_num + 2;

	//
	MatMPIAIJSetPreallocation(ddz_matrix111111_to_000000, 3, NULL, 3, NULL);
	MatSeqAIJSetPreallocation(ddz_matrix111111_to_000000, 3, NULL);
	MatZeroEntries(ddz_matrix111111_to_000000);

	MatGetOwnershipRange(ddz_matrix111111_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 2)*(ghost_n_num - 2));
		PetscInt i = Ii02 / (ghost_m_num - 2);
		PetscInt j = Ii02 - i * (ghost_m_num - 2);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 1);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 2)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 1);

		int Ij_vector[3] = { Ij_0, Ij_1, Ij_2 };
		double insert_value_vector[3] = { 1.0 / (hz * hz), -2.0 / (hz * hz), 1.0 / (hz * hz) };

		MatSetValues(ddz_matrix111111_to_000000, 1, &Ii, 3, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(ddz_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ddz_matrix111111_to_000000, MAT_FINAL_ASSEMBLY);
}

//***************************************************************************************************************************************************

//平均数据============================================================================================================================================

//***************************************************************************************************************************************************

//将(100000或010000类型数据求两点平均000000)=======================================================
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_average_matrix100000_or_010000_to_000000()
{
	int ghost_m_num = m_num + 1;
	int ghost_n_num = n_num + 0;
	int ghost_o_num = o_num + 0;

	//
	MatMPIAIJSetPreallocation(average_matrix100000_or_010000_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(average_matrix100000_or_010000_to_000000, 2, NULL);
	MatZeroEntries(average_matrix100000_or_010000_to_000000);

	MatGetOwnershipRange(average_matrix100000_or_010000_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 1)*(ghost_n_num - 0));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 1)*(ghost_n_num - 0));
		PetscInt i = Ii02 / (ghost_m_num - 1);
		PetscInt j = Ii02 - i * (ghost_m_num - 1);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 1);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { 1.0 / 2.0, 1.0 / 2.0 };

		MatSetValues(average_matrix100000_or_010000_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(average_matrix100000_or_010000_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(average_matrix100000_or_010000_to_000000, MAT_FINAL_ASSEMBLY);
}

//将(001000或000100类型数据求两点平均000000)=======================================================
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_average_matrix001000_or_000100_to_000000()
{
	int ghost_m_num = m_num + 0;
	int ghost_n_num = n_num + 1;
	int ghost_o_num = o_num + 0;

	//
	MatMPIAIJSetPreallocation(average_matrix001000_or_000100_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(average_matrix001000_or_000100_to_000000, 2, NULL);
	MatZeroEntries(average_matrix001000_or_000100_to_000000);

	MatGetOwnershipRange(average_matrix001000_or_000100_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 0)*(ghost_n_num - 1));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 0)*(ghost_n_num - 1));
		PetscInt i = Ii02 / (ghost_m_num - 0);
		PetscInt j = Ii02 - i * (ghost_m_num - 0);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 0)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 0);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { 1.0 / 2.0, 1.0 / 2.0 };

		MatSetValues(average_matrix001000_or_000100_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(average_matrix001000_or_000100_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(average_matrix001000_or_000100_to_000000, MAT_FINAL_ASSEMBLY);
}

//将(000010或000001类型数据求两点平均000000)=======================================================
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_average_matrix000010_or_000001_to_000000()
{
	int ghost_m_num = m_num + 0;
	int ghost_n_num = n_num + 0;
	int ghost_o_num = o_num + 1;

	//
	MatMPIAIJSetPreallocation(average_matrix000010_or_000001_to_000000, 2, NULL, 2, NULL);
	MatSeqAIJSetPreallocation(average_matrix000010_or_000001_to_000000, 2, NULL);
	MatZeroEntries(average_matrix000010_or_000001_to_000000);

	MatGetOwnershipRange(average_matrix000010_or_000001_to_000000, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 0)*(ghost_n_num - 0));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 0)*(ghost_n_num - 0));
		PetscInt i = Ii02 / (ghost_m_num - 0);
		PetscInt j = Ii02 - i * (ghost_m_num - 0);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 1)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);

		int Ij_vector[2] = { Ij_0, Ij_1 };
		double insert_value_vector[2] = { 1.0 / 2.0, 1.0 / 2.0 };

		MatSetValues(average_matrix000010_or_000001_to_000000, 1, &Ii, 2, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(average_matrix000010_or_000001_to_000000, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(average_matrix000010_or_000001_to_000000, MAT_FINAL_ASSEMBLY);
}

//将(344333或433433类型数据求四点平均333333)=======================================================
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_average_matrix344333_or_433433_to_333333()
{
	//边界
	int ghost_m_num = m_num + 7;
	int ghost_n_num = n_num + 7;
	int ghost_o_num = o_num + 6;

	//组装
	MatMPIAIJSetPreallocation(average_matrix344333_or_433433_to_333333, 4, NULL, 4, NULL);
	MatSeqAIJSetPreallocation(average_matrix344333_or_433433_to_333333, 4, NULL);
	MatZeroEntries(average_matrix344333_or_433433_to_333333);

	MatGetOwnershipRange(average_matrix344333_or_433433_to_333333, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 1)*(ghost_n_num - 1));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 1)*(ghost_n_num - 1));
		PetscInt i = Ii02 / (ghost_m_num - 1);
		PetscInt j = Ii02 - i * (ghost_m_num - 1);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 0)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 0);
		int Ij_3 = (k + 0)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 1);

		int Ij_vector[4] = { Ij_0, Ij_1, Ij_2, Ij_3 };
		double insert_value_vector[4] = { 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0 };

		MatSetValues(average_matrix344333_or_433433_to_333333, 1, &Ii, 4, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(average_matrix344333_or_433433_to_333333, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(average_matrix344333_or_433433_to_333333, MAT_FINAL_ASSEMBLY);
}

//将(343343或433334类型数据求四点平均333333)=======================================================
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_average_matrix343343_or_433334_to_333333()
{
	//边界
	int ghost_m_num = m_num + 7;
	int ghost_n_num = n_num + 6;
	int ghost_o_num = o_num + 7;

	//组装
	MatMPIAIJSetPreallocation(average_matrix343343_or_433334_to_333333, 4, NULL, 4, NULL);
	MatSeqAIJSetPreallocation(average_matrix343343_or_433334_to_333333, 4, NULL);
	MatZeroEntries(average_matrix343343_or_433334_to_333333);

	MatGetOwnershipRange(average_matrix343343_or_433334_to_333333, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 1)*(ghost_n_num - 0));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 1)*(ghost_n_num - 0));
		PetscInt i = Ii02 / (ghost_m_num - 1);
		PetscInt j = Ii02 - i * (ghost_m_num - 1);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 1);
		int Ij_2 = (k + 1)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_3 = (k + 1)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 1);

		int Ij_vector[4] = { Ij_0, Ij_1, Ij_2, Ij_3 };
		double insert_value_vector[4] = { 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0 };

		MatSetValues(average_matrix343343_or_433334_to_333333, 1, &Ii, 4, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(average_matrix343343_or_433334_to_333333, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(average_matrix343343_or_433334_to_333333, MAT_FINAL_ASSEMBLY);
}

//将(333443或334334类型数据求四点平均333333)=======================================================
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_average_matrix333443_or_334334_to_333333()
{
	//边界
	int ghost_m_num = m_num + 6;
	int ghost_n_num = n_num + 7;
	int ghost_o_num = o_num + 7;

	//组装
	MatMPIAIJSetPreallocation(average_matrix333443_or_334334_to_333333, 4, NULL, 4, NULL);
	MatSeqAIJSetPreallocation(average_matrix333443_or_334334_to_333333, 4, NULL);
	MatZeroEntries(average_matrix333443_or_334334_to_333333);

	MatGetOwnershipRange(average_matrix333443_or_334334_to_333333, &istart, &iend);
	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / ((ghost_m_num - 0)*(ghost_n_num - 1));
		PetscInt Ii02 = Ii - k*((ghost_m_num - 0)*(ghost_n_num - 1));
		PetscInt i = Ii02 / (ghost_m_num - 0);
		PetscInt j = Ii02 - i * (ghost_m_num - 0);

		int Ij_0 = (k + 0)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 0);
		int Ij_1 = (k + 0)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);
		int Ij_2 = (k + 1)*ghost_m_num*ghost_n_num + (i + 1)*ghost_m_num + (j + 0);
		int Ij_3 = (k + 1)*ghost_m_num*ghost_n_num + (i + 0)*ghost_m_num + (j + 0);

		int Ij_vector[4] = { Ij_0, Ij_1, Ij_2, Ij_3 };
		double insert_value_vector[4] = { 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0 };

		MatSetValues(average_matrix333443_or_334334_to_333333, 1, &Ii, 4, Ij_vector, insert_value_vector, INSERT_VALUES);
	}
	MatAssemblyBegin(average_matrix333443_or_334334_to_333333, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(average_matrix333443_or_334334_to_333333, MAT_FINAL_ASSEMBLY);
}

//***************************************************************************************************************************************************

//拓展数据点==========================================================================================================================================

//***************************************************************************************************************************************************

//拓展矩阵==========================================================================================

//非交错网格--------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_nonstagger_extend_velocity_u_matrix3D(Mat extend_u_matrix, int *extend_num, int velocity_boundary)
{
	//边界
	int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
	int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
	int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
	int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
	int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(extend_u_matrix, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(extend_u_matrix, 1, NULL);
	MatZeroEntries(extend_u_matrix);
	MatGetOwnershipRange(extend_u_matrix, &istart, &iend);

	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (ghost_m_num*ghost_n_num);
		PetscInt i = (Ii - k*ghost_m_num*ghost_n_num) / ghost_m_num;
		PetscInt j = (Ii - k*ghost_m_num*ghost_n_num) - i*ghost_m_num;
		PetscInt Ij;

		int j_mark;
		int i_mark;
		int k_mark;
		double insert_value;
		int j_judge = 0;
		int i_judge = 0;
		int k_judge = 0;

		//0边界条件-----
		if (velocity_boundary == 0){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0;
				MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//周期性边界条件
		if (velocity_boundary == 1){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_mark = m_num - (ex_j_L - j);
				}
				else{
					j_mark = j - (m_num + ex_j_L);
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_mark = n_num - (ex_i_L - i);
				}
				else{
					i_mark = i - (n_num + ex_i_L);
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_mark = o_num - (ex_k_L - k);
				}
				else{
					k_mark = k - (o_num + ex_k_L);
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//固壁边界条件--
		if (velocity_boundary == 2){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//方腔驱动流----
		if (velocity_boundary == 3){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//横向管道流----
		if (velocity_boundary == 4){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_judge = -1;
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_judge = 1;
					j_mark = m_num - 1;
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(extend_u_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(extend_u_matrix, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_nonstagger_extend_velocity_v_matrix3D(Mat extend_v_matrix, int *extend_num, int velocity_boundary)
{
	//边界
	int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
	int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
	int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
	int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
	int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(extend_v_matrix, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(extend_v_matrix, 1, NULL);
	MatZeroEntries(extend_v_matrix);
	MatGetOwnershipRange(extend_v_matrix, &istart, &iend);

	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (ghost_m_num*ghost_n_num);
		PetscInt i = (Ii - k*ghost_m_num*ghost_n_num) / ghost_m_num;
		PetscInt j = (Ii - k*ghost_m_num*ghost_n_num) - i*ghost_m_num;
		PetscInt Ij;

		int j_mark;
		int i_mark;
		int k_mark;
		double insert_value;
		int j_judge = 0;
		int i_judge = 0;
		int k_judge = 0;

		//0边界条件-----
		if (velocity_boundary == 0){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0;
				MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//周期性边界条件
		if (velocity_boundary == 1){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_mark = m_num - (ex_j_L - j);
				}
				else{
					j_mark = j - (m_num + ex_j_L);
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_mark = n_num - (ex_i_L - i);
				}
				else{
					i_mark = i - (n_num + ex_i_L);
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_mark = o_num - (ex_k_L - k);
				}
				else{
					k_mark = k - (o_num + ex_k_L);
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//固壁边界条件--
		if (velocity_boundary == 2){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//方腔驱动流----
		if (velocity_boundary == 3){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//横向管道流----
		if (velocity_boundary == 4){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_judge = -1;
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_judge = 1;
					j_mark = m_num - 1;
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(extend_v_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(extend_v_matrix, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_nonstagger_extend_velocity_w_matrix3D(Mat extend_w_matrix, int *extend_num, int velocity_boundary)
{
	//边界
	int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
	int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
	int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
	int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
	int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(extend_w_matrix, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(extend_w_matrix, 1, NULL);
	MatZeroEntries(extend_w_matrix);
	MatGetOwnershipRange(extend_w_matrix, &istart, &iend);

	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (ghost_m_num*ghost_n_num);
		PetscInt i = (Ii - k*ghost_m_num*ghost_n_num) / ghost_m_num;
		PetscInt j = (Ii - k*ghost_m_num*ghost_n_num) - i*ghost_m_num;
		PetscInt Ij;

		int j_mark;
		int i_mark;
		int k_mark;
		double insert_value;
		int j_judge = 0;
		int i_judge = 0;
		int k_judge = 0;

		//0边界条件-----
		if (velocity_boundary == 0){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0;
				MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//周期性边界条件
		if (velocity_boundary == 1){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_mark = m_num - (ex_j_L - j);
				}
				else{
					j_mark = j - (m_num + ex_j_L);
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_mark = n_num - (ex_i_L - i);
				}
				else{
					i_mark = i - (n_num + ex_i_L);
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_mark = o_num - (ex_k_L - k);
				}
				else{
					k_mark = k - (o_num + ex_k_L);
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//固壁边界条件--
		if (velocity_boundary == 2){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//方腔驱动流----
		if (velocity_boundary == 3){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//横向管道流----
		if (velocity_boundary == 4){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_judge = -1;
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_judge = 1;
					j_mark = m_num - 1;
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(extend_w_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(extend_w_matrix, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_nonstagger_extend_pressure_matrix3D(Mat extend_pressure_matrix, int *extend_num, int pressure_boundary)
{
	//边界
	int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
	int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
	int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
	int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
	int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(extend_pressure_matrix, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(extend_pressure_matrix, 1, NULL);
	MatZeroEntries(extend_pressure_matrix);
	MatGetOwnershipRange(extend_pressure_matrix, &istart, &iend);

	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (ghost_m_num*ghost_n_num);
		PetscInt i = (Ii - k*ghost_m_num*ghost_n_num) / ghost_m_num;
		PetscInt j = (Ii - k*ghost_m_num*ghost_n_num) - i*ghost_m_num;
		PetscInt Ij;

		int j_mark;
		int i_mark;
		int k_mark;
		double insert_value;
		int j_judge = 0;
		int i_judge = 0;
		int k_judge = 0;

		//0边界条件-----
		if (pressure_boundary == 0){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0;
				MatSetValues(extend_pressure_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//周期性边界条件
		if (pressure_boundary == 1){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_mark = m_num - (ex_j_L - j);
				}
				else{
					j_mark = j - (m_num + ex_j_L);
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_mark = n_num - (ex_i_L - i);
				}
				else{
					i_mark = i - (n_num + ex_i_L);
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_mark = o_num - (ex_k_L - k);
				}
				else{
					k_mark = k - (o_num + ex_k_L);
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_pressure_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//固壁边界条件--
		if (pressure_boundary == 2){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_pressure_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(extend_pressure_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(extend_pressure_matrix, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_nonstagger_extend_phi_matrix3D(Mat extend_phi_matrix, int *extend_num, int phi_boundary)
{
	//边界
	int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
	int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
	int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
	int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
	int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(extend_phi_matrix, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(extend_phi_matrix, 1, NULL);
	MatZeroEntries(extend_phi_matrix);
	MatGetOwnershipRange(extend_phi_matrix, &istart, &iend);

	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (ghost_m_num*ghost_n_num);
		PetscInt i = (Ii - k*ghost_m_num*ghost_n_num) / ghost_m_num;
		PetscInt j = (Ii - k*ghost_m_num*ghost_n_num) - i*ghost_m_num;
		PetscInt Ij;

		int j_mark;
		int i_mark;
		int k_mark;
		double insert_value;
		int j_judge = 0;
		int i_judge = 0;
		int k_judge = 0;

		//0边界条件-----
		if (phi_boundary == 0){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0;
				MatSetValues(extend_phi_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//周期性边界条件
		if (phi_boundary == 1){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_mark = m_num - (ex_j_L - j);
				}
				else{
					j_mark = j - (m_num + ex_j_L);
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_mark = n_num - (ex_i_L - i);
				}
				else{
					i_mark = i - (n_num + ex_i_L);
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_mark = o_num - (ex_k_L - k);
				}
				else{
					k_mark = k - (o_num + ex_k_L);
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_phi_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(extend_phi_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(extend_phi_matrix, MAT_FINAL_ASSEMBLY);
}

//交错网格----------------------------------------
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_stagger_extend_velocity_u_matrix3D(Mat extend_u_matrix, int *extend_num, int velocity_boundary)
{
	//边界
	int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
	int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
	int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
	int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
	int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(extend_u_matrix, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(extend_u_matrix, 1, NULL);
	MatZeroEntries(extend_u_matrix);
	MatGetOwnershipRange(extend_u_matrix, &istart, &iend);

	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (ghost_m_num*ghost_n_num);
		PetscInt i = (Ii - k * ghost_m_num*ghost_n_num) / ghost_m_num;
		PetscInt j = (Ii - k * ghost_m_num*ghost_n_num) - i*ghost_m_num;
		PetscInt Ij;

		int j_mark;
		int i_mark;
		int k_mark;
		double insert_value;
		int j_judge = 0;
		int i_judge = 0;
		int k_judge = 0;

		//0边界条件-----
		if (velocity_boundary == 0){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0;
				MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//周期性边界条件
		if (velocity_boundary == 1){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_mark = m_num - (ex_j_L - j);
				}
				else{
					j_mark = j - (m_num + ex_j_L);
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_mark = n_num - (ex_i_L - i);
				}
				else{
					i_mark = i - (n_num + ex_i_L);
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_mark = o_num - (ex_k_L - k);
				}
				else{
					k_mark = k - (o_num + ex_k_L);
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//固壁边界条件--
		if (velocity_boundary == 2){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_judge = -1;
					j_mark = ex_j_L - j;
				}
				else{
					if (j == m_num + ex_j_L){
						j_judge = 0;
					}
					else{
						j_judge = -1;
						j_mark = (m_num - 1) - (j - (m_num + ex_j_L)) + 1;
					}
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			if (j_judge != 0){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0 * j_judge * i_judge * k_judge;
				MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//方腔驱动流----
		if (velocity_boundary == 3){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_judge = -1;
					j_mark = ex_j_L - j;
				}
				else{
					if (j == m_num + ex_j_L){
						j_judge = 0;
					}
					else{
						j_judge = -1;
						j_mark = (m_num - 1) - (j - (m_num + ex_j_L)) + 1;
					}
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			if (j_judge != 0){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0 * j_judge * i_judge * k_judge;
				MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//横向管道流----
		if (velocity_boundary == 4){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_judge = -1;
					j_mark = ex_j_L - j;
				}
				else{
					j_judge = 1;
					j_mark = m_num - 1;
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0 * j_judge * i_judge * k_judge;
			MatSetValues(extend_u_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(extend_u_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(extend_u_matrix, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_stagger_extend_velocity_v_matrix3D(Mat extend_v_matrix, int *extend_num, int velocity_boundary)
{
	//边界
	int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
	int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
	int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
	int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
	int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(extend_v_matrix, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(extend_v_matrix, 1, NULL);
	MatZeroEntries(extend_v_matrix);
	MatGetOwnershipRange(extend_v_matrix, &istart, &iend);

	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (ghost_m_num*ghost_n_num);
		PetscInt i = (Ii - k * ghost_m_num*ghost_n_num) / ghost_m_num;
		PetscInt j = (Ii - k * ghost_m_num*ghost_n_num) - i*ghost_m_num;
		PetscInt Ij;

		int j_mark;
		int i_mark;
		int k_mark;
		double insert_value;
		int j_judge = 0;
		int i_judge = 0;
		int k_judge = 0;

		//0边界条件-----
		if (velocity_boundary == 0){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0;
				MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//周期性边界条件
		if (velocity_boundary == 1){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_mark = m_num - (ex_j_L - j);
				}
				else{
					j_mark = j - (m_num + ex_j_L);
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_mark = n_num - (ex_i_L - i);
				}
				else{
					i_mark = i - (n_num + ex_i_L);
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_mark = o_num - (ex_k_L - k);
				}
				else{
					k_mark = k - (o_num + ex_k_L);
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//固壁边界条件--
		if (velocity_boundary == 2){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_judge = -1;
					i_mark = ex_i_L - i;
				}
				else{
					if (i == n_num + ex_i_L){
						i_judge = 0;
					}
					else{
						i_judge = -1;
						i_mark = (n_num - 1) - (i - (n_num + ex_i_L)) + 1;
					}
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			if (i_judge != 0){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0 * j_judge * i_judge * k_judge;
				MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//方腔驱动流----
		if (velocity_boundary == 3){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_judge = -1;
					i_mark = ex_i_L - i;
				}
				else{
					if (i == n_num + ex_i_L){
						i_judge = 0;
					}
					else{
						i_judge = -1;
						i_mark = (n_num - 1) - (i - (n_num + ex_i_L)) + 1;
					}
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			if (i_judge != 0){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0 * j_judge * i_judge * k_judge;
				MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//横向管道流----
		if (velocity_boundary == 4){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_judge = -1;
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_judge = 1;
					j_mark = m_num - 1;
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_judge = -1;
					i_mark = ex_i_L - i;
				}
				else{
					if (i == n_num + ex_i_L){
						i_judge = 0;
					}
					else{
						i_judge = -1;
						i_mark = (n_num - 1) - (i - (n_num + ex_i_L)) + 1;
					}
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				k_judge = -1;
				if (k < ghost_kL_mark){
					k_mark = (ex_k_L - k) - 1;
				}
				else{
					k_mark = (o_num - 1) - (k - (o_num + ex_k_L));
				}
			}

			//
			if (i_judge != 0){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0 * j_judge * i_judge * k_judge;
				MatSetValues(extend_v_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}
	}

	MatAssemblyBegin(extend_v_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(extend_v_matrix, MAT_FINAL_ASSEMBLY);
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_stagger_extend_velocity_w_matrix3D(Mat extend_w_matrix, int *extend_num, int velocity_boundary)
{
	//边界
	int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
	int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
	int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
	int ghost_m_num = m_num + ex_j_L + ex_j_H;
	int ghost_n_num = n_num + ex_i_L + ex_i_H;
	int ghost_o_num = o_num + ex_k_L + ex_k_H;

	int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
	int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
	int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

	//组装
	MatMPIAIJSetPreallocation(extend_w_matrix, 1, NULL, 1, NULL);
	MatSeqAIJSetPreallocation(extend_w_matrix, 1, NULL);
	MatZeroEntries(extend_w_matrix);
	MatGetOwnershipRange(extend_w_matrix, &istart, &iend);

	for (PetscInt Ii = istart; Ii < iend; Ii++){
		PetscInt k = Ii / (ghost_m_num*ghost_n_num);
		PetscInt i = (Ii - k * ghost_m_num*ghost_n_num) / ghost_m_num;
		PetscInt j = (Ii - k * ghost_m_num*ghost_n_num) - i*ghost_m_num;
		PetscInt Ij;

		int j_mark;
		int i_mark;
		int k_mark;
		double insert_value;
		int j_judge = 0;
		int i_judge = 0;
		int k_judge = 0;

		//0边界条件-----
		if (velocity_boundary == 0){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0;
				MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//周期性边界条件
		if (velocity_boundary == 1){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_mark = m_num - (ex_j_L - j);
				}
				else{
					j_mark = j - (m_num + ex_j_L);
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_mark = i - ex_i_L;
			}
			else{
				if (i < ghost_iL_mark){
					i_mark = n_num - (ex_i_L - i);
				}
				else{
					i_mark = i - (n_num + ex_i_L);
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_mark = o_num - (ex_k_L - k);
				}
				else{
					k_mark = k - (o_num + ex_k_L);
				}
			}

			//
			Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
			insert_value = 1.0;
			MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
		}

		//固壁边界条件--
		if (velocity_boundary == 2){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_judge = -1;
					k_mark = ex_k_L - k;
				}
				else{
					if (k == o_num + ex_k_L){
						k_judge = 0;
					}
					else{
						k_judge = -1;
						k_mark = (o_num - 1) - (k - (o_num + ex_k_L)) + 1;
					}
				}
			}

			//
			if (k_judge != 0){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0 * j_judge * i_judge * k_judge;
				MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//方腔驱动流----
		if (velocity_boundary == 3){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				j_judge = -1;
				if (j < ghost_jL_mark){
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_mark = (m_num - 1) - (j - (m_num + ex_j_L));
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_judge = -1;
					k_mark = ex_k_L - k;
				}
				else{
					if (k == o_num + ex_k_L){
						k_judge = 0;
					}
					else{
						k_judge = -1;
						k_mark = (o_num - 1) - (k - (o_num + ex_k_L)) + 1;
					}
				}
			}

			//
			if (k_judge != 0){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0 * j_judge * i_judge * k_judge;
				MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}

		//横向管道流----
		if (velocity_boundary == 4){
			//x
			if (j >= ghost_jL_mark && j <= ghost_jH_mark){
				j_judge = 1;
				j_mark = j - ex_j_L;
			}
			else{
				if (j < ghost_jL_mark){
					j_judge = -1;
					j_mark = (ex_j_L - j) - 1;
				}
				else{
					j_judge = 1;
					j_mark = m_num - 1;
				}
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
				i_mark = i - ex_i_L;
			}
			else{
				i_judge = -1;
				if (i < ghost_iL_mark){
					i_mark = (ex_i_L - i) - 1;
				}
				else{
					i_mark = (n_num - 1) - (i - (n_num + ex_i_L));
				}
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
				k_mark = k - ex_k_L;
			}
			else{
				if (k < ghost_kL_mark){
					k_judge = -1;
					k_mark = ex_k_L - k;
				}
				else{
					if (k == o_num + ex_k_L){
						k_judge = 0;
					}
					else{
						k_judge = -1;
						k_mark = (o_num - 1) - (k - (o_num + ex_k_L)) + 1;
					}
				}
			}

			//
			if (k_judge != 0){
				Ij = k_mark*n_num*m_num + i_mark*m_num + j_mark;
				insert_value = 1.0 * j_judge * i_judge * k_judge;
				MatSetValues(extend_w_matrix, 1, &Ii, 1, &Ij, &insert_value, INSERT_VALUES);
			}
		}
	}

	MatAssemblyBegin(extend_w_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(extend_w_matrix, MAT_FINAL_ASSEMBLY);
}

//拓展向量==========================================================================================
void wyh_solve_fluid_structure_interaction_system3D200630::assembly_extend_velocity_u_auxiliary_vector3D(Vec u_auxiliary_vector, int *extend_num, int velocity_boundary)
{
	if (velocity_boundary == 0){
		VecSet(u_auxiliary_vector, 0);
	}

	if (velocity_boundary == 1){
		VecSet(u_auxiliary_vector, 0);
	}

	if (velocity_boundary == 2){
		VecSet(u_auxiliary_vector, 0);
	}

	if (velocity_boundary == 3){
		//边界
		int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
		int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
		int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
		int ghost_m_num = m_num + ex_j_L + ex_j_H;
		int ghost_n_num = n_num + ex_i_L + ex_i_H;
		int ghost_o_num = o_num + ex_k_L + ex_k_H;

		int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
		int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
		int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

		//组装
		VecGetOwnershipRange(u_auxiliary_vector, &istart, &iend);
		for (PetscInt Ii = istart; Ii < iend; Ii++) {
			PetscInt k = Ii / (ghost_m_num*ghost_n_num);
			PetscInt i = (Ii - k*ghost_m_num*ghost_n_num) / ghost_m_num;
			PetscInt j = (Ii - k*ghost_m_num*ghost_n_num) - i*ghost_m_num;
			PetscInt Ij;

			int j_mark;
			int i_mark;
			int k_mark;
			double insert_value;
			int j_judge = 0;
			int i_judge = 0;
			int k_judge = 0;

			//x
			//if (j >= ghost_jL_mark && j <= ghost_jH_mark){
			j_judge = 1;
			//}

			//y
			//if (i >= ghost_iL_mark && i <= ghost_iH_mark){
			i_judge = 1;
			//}

			//z
			if (k > ghost_kH_mark){
				k_judge = 1;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				insert_value = 2.0;
				VecSetValues(u_auxiliary_vector, 1, &Ii, &insert_value, INSERT_VALUES);
			}
		}
		VecAssemblyBegin(u_auxiliary_vector);
		VecAssemblyEnd(u_auxiliary_vector);
	}

	if (velocity_boundary == 4){
		//边界
		int ex_j_L = extend_num[0];		int ex_j_H = extend_num[1];
		int ex_i_L = extend_num[2];		int ex_i_H = extend_num[3];
		int ex_k_L = extend_num[4];		int ex_k_H = extend_num[5];
		int ghost_m_num = m_num + ex_j_L + ex_j_H;
		int ghost_n_num = n_num + ex_i_L + ex_i_H;
		int ghost_o_num = o_num + ex_k_L + ex_k_H;

		int ghost_jL_mark = ex_j_L;		int ghost_jH_mark = ex_j_L + m_num - 1;
		int ghost_iL_mark = ex_i_L;		int ghost_iH_mark = ex_i_L + n_num - 1;
		int ghost_kL_mark = ex_k_L;		int ghost_kH_mark = ex_k_L + o_num - 1;

		//组装
		VecGetOwnershipRange(u_auxiliary_vector, &istart, &iend);
		for (PetscInt Ii = istart; Ii < iend; Ii++) {
			PetscInt k = Ii / (ghost_m_num*ghost_n_num);
			PetscInt i = (Ii - k*ghost_m_num*ghost_n_num) / ghost_m_num;
			PetscInt j = (Ii - k*ghost_m_num*ghost_n_num) - i*ghost_m_num;
			PetscInt Ij;

			int j_mark;
			int i_mark;
			int k_mark;
			double insert_value;
			int j_judge = 0;
			int i_judge = 0;
			int k_judge = 0;

			//x
			if (j > ghost_jH_mark){
				j_judge = 1;
			}

			//y
			if (i >= ghost_iL_mark && i <= ghost_iH_mark){
				i_judge = 1;
			}

			//z
			if (k >= ghost_kL_mark && k <= ghost_kH_mark){
				k_judge = 1;
			}

			//
			if (j_judge == 1 && i_judge == 1 && k_judge == 1){
				insert_value = 2.0;
				VecSetValues(u_auxiliary_vector, 1, &Ii, &insert_value, INSERT_VALUES);
			}
		}
		VecAssemblyBegin(u_auxiliary_vector);
		VecAssemblyEnd(u_auxiliary_vector);
	}
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_extend_velocity_v_auxiliary_vector3D(Vec v_auxiliary_vector, int *extend_num, int velocity_boundary)
{
	if (velocity_boundary == 0){
		VecSet(v_auxiliary_vector, 0);
	}

	if (velocity_boundary == 1){
		VecSet(v_auxiliary_vector, 0);
	}

	if (velocity_boundary == 2){
		VecSet(v_auxiliary_vector, 0);
	}

	if (velocity_boundary == 3){
		VecSet(v_auxiliary_vector, 0);
	}

	if (velocity_boundary == 4){
		VecSet(v_auxiliary_vector, 0);
	}
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_extend_velocity_w_auxiliary_vector3D(Vec w_auxiliary_vector, int *extend_num, int velocity_boundary)
{
	if (velocity_boundary == 0){
		VecSet(w_auxiliary_vector, 0);
	}

	if (velocity_boundary == 1){
		VecSet(w_auxiliary_vector, 0);
	}

	if (velocity_boundary == 2){
		VecSet(w_auxiliary_vector, 0);
	}

	if (velocity_boundary == 3){
		VecSet(w_auxiliary_vector, 0);
	}

	if (velocity_boundary == 4){
		VecSet(w_auxiliary_vector, 0);
	}
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_extend_pressure_auxiliary_vector3D(Vec pressure_auxiliary_vector, int *extend_num, int pressure_boundary)
{
	if (pressure_boundary == 0){
		VecSet(pressure_auxiliary_vector, 0);
	}

	if (pressure_boundary == 1){
		VecSet(pressure_auxiliary_vector, 0);
	}

	if (pressure_boundary == 2){
		VecSet(pressure_auxiliary_vector, 0);
	}

	if (pressure_boundary == 3){
		VecSet(pressure_auxiliary_vector, 0);
	}

	if (pressure_boundary == 4){
		VecSet(pressure_auxiliary_vector, 0);
	}
}

void wyh_solve_fluid_structure_interaction_system3D200630::assembly_extend_phi_auxiliary_vector3D(Vec phi_auxiliary_vector, int *extend_num, int phi_boundary)
{
	if (phi_boundary == 0){
		VecSet(phi_auxiliary_vector, 0);
	}

	if (phi_boundary == 1){
		VecSet(phi_auxiliary_vector, 0);
	}

	if (phi_boundary == 2){
		VecSet(phi_auxiliary_vector, 0);
	}

	if (phi_boundary == 3){
		VecSet(phi_auxiliary_vector, 0);
	}

	if (phi_boundary == 4){
		VecSet(phi_auxiliary_vector, 0);
	}
}

//***************************************************************************************************************************************************

//Destroy============================================================================================================================================

//***************************************************************************************************************************************************
void wyh_solve_fluid_structure_interaction_system3D200630::destroy_variable()
{
	//KSP
	KSPDestroy(&KSP_NSeq_u);
	KSPDestroy(&KSP_NSeq_v);
	KSPDestroy(&KSP_NSeq_w);
	KSPDestroy(&KSP_NSeq_p);

	//速度
	VecDestroy(&velocity_u);			VecDestroy(&velocity_v);			VecDestroy(&velocity_w);
	VecDestroy(&velocity_u_at_p);		VecDestroy(&velocity_v_at_p);		VecDestroy(&velocity_w_at_p);

	//压力
	VecDestroy(&pressure);

	//phi
	VecDestroy(&phi);

	//FF
	VecDestroy(&FF00);					VecDestroy(&FF01);					VecDestroy(&FF02);
	VecDestroy(&FF10);					VecDestroy(&FF11);					VecDestroy(&FF12);
	VecDestroy(&FF20);					VecDestroy(&FF21);					VecDestroy(&FF22);

	// 拓展矩阵----------------------------------

	//velocity
	MatDestroy(&nonstagger_extend_u_matrix111111);		MatDestroy(&nonstagger_extend_v_matrix111111);		MatDestroy(&nonstagger_extend_w_matrix111111);
	MatDestroy(&nonstagger_extend_u_matrix333333);		MatDestroy(&nonstagger_extend_v_matrix333333);		MatDestroy(&nonstagger_extend_w_matrix333333);

	MatDestroy(&extend_u_matrix010000);			MatDestroy(&extend_v_matrix000100);			MatDestroy(&extend_w_matrix000001);
	VecDestroy(&u_auxiliary_vector010000);		//VecDestroy(&v_auxiliary_vector000100);		VecDestroy(&w_auxiliary_vector000001);

	MatDestroy(&extend_u_matrix111111);			MatDestroy(&extend_v_matrix111111);			MatDestroy(&extend_w_matrix111111);
	VecDestroy(&u_auxiliary_vector111111);		//VecDestroy(&v_auxiliary_vector111111);		VecDestroy(&w_auxiliary_vector111111);

	MatDestroy(&extend_u_matrix333333);			MatDestroy(&extend_v_matrix333333);			MatDestroy(&extend_w_matrix333333);
	VecDestroy(&u_auxiliary_vector333333);		//VecDestroy(&v_auxiliary_vector333333);		VecDestroy(&w_auxiliary_vector333333);

	MatDestroy(&extend_u_matrix344333);			MatDestroy(&extend_u_matrix343343);
	VecDestroy(&u_auxiliary_vector344333);		VecDestroy(&u_auxiliary_vector343343);

	MatDestroy(&extend_v_matrix433433);			MatDestroy(&extend_v_matrix333443);
	//VecDestroy(&v_auxiliary_vector433433);		VecDestroy(&v_auxiliary_vector333443);

	MatDestroy(&extend_w_matrix433334);			MatDestroy(&extend_w_matrix334334);
	//VecDestroy(&w_auxiliary_vector433334);		VecDestroy(&w_auxiliary_vector334334);
	
	//phi
	MatDestroy(&extend_phi_matrix100000);			MatDestroy(&extend_phi_matrix001000);		MatDestroy(&extend_phi_matrix000010);
	//VecDestroy(&phi_auxiliary_vector100000);		VecDestroy(&phi_auxiliary_vector001000);	VecDestroy(&phi_auxiliary_vector000010);

	MatDestroy(&extend_phi_matrix111111);
	//VecDestroy(&phi_auxiliary_vector111111);

	MatDestroy(&extend_phi_matrix333333);
	//VecDestroy(&phi_auxiliary_vector333333);

	//WENO 相关矩阵------------------------------
	MatDestroy(&f_P0_x_matrix333333_to_100000);	MatDestroy(&f_P0_y_matrix333333_to_001000);	MatDestroy(&f_P0_z_matrix333333_to_000010);
	MatDestroy(&f_P1_x_matrix333333_to_100000);	MatDestroy(&f_P1_y_matrix333333_to_001000);	MatDestroy(&f_P1_z_matrix333333_to_000010);
	MatDestroy(&f_P2_x_matrix333333_to_100000);	MatDestroy(&f_P2_y_matrix333333_to_001000);	MatDestroy(&f_P2_z_matrix333333_to_000010);

	MatDestroy(&f_N0_x_matrix333333_to_100000);	MatDestroy(&f_N0_y_matrix333333_to_001000);	MatDestroy(&f_N0_z_matrix333333_to_000010);
	MatDestroy(&f_N1_x_matrix333333_to_100000);	MatDestroy(&f_N1_y_matrix333333_to_001000);	MatDestroy(&f_N1_z_matrix333333_to_000010);
	MatDestroy(&f_N2_x_matrix333333_to_100000);	MatDestroy(&f_N2_y_matrix333333_to_001000);	MatDestroy(&f_N2_z_matrix333333_to_000010);

	MatDestroy(&ISP0_x_part1_matrix333333_to_100000);	MatDestroy(&ISP0_y_part1_matrix333333_to_001000);	MatDestroy(&ISP0_z_part1_matrix333333_to_000010);
	MatDestroy(&ISP1_x_part1_matrix333333_to_100000);	MatDestroy(&ISP1_y_part1_matrix333333_to_001000);	MatDestroy(&ISP1_z_part1_matrix333333_to_000010);
	MatDestroy(&ISP2_x_part1_matrix333333_to_100000);	MatDestroy(&ISP2_y_part1_matrix333333_to_001000);	MatDestroy(&ISP2_z_part1_matrix333333_to_000010);

	MatDestroy(&ISP0_x_part2_matrix333333_to_100000);	MatDestroy(&ISP0_y_part2_matrix333333_to_001000);	MatDestroy(&ISP0_z_part2_matrix333333_to_000010);
	MatDestroy(&ISP1_x_part2_matrix333333_to_100000);	MatDestroy(&ISP1_y_part2_matrix333333_to_001000);	MatDestroy(&ISP1_z_part2_matrix333333_to_000010);
	MatDestroy(&ISP2_x_part2_matrix333333_to_100000);	MatDestroy(&ISP2_y_part2_matrix333333_to_001000);	MatDestroy(&ISP2_z_part2_matrix333333_to_000010);

	MatDestroy(&ISN0_x_part1_matrix333333_to_100000);	MatDestroy(&ISN0_y_part1_matrix333333_to_001000);	MatDestroy(&ISN0_z_part1_matrix333333_to_000010);
	MatDestroy(&ISN1_x_part1_matrix333333_to_100000);	MatDestroy(&ISN1_y_part1_matrix333333_to_001000);	MatDestroy(&ISN1_z_part1_matrix333333_to_000010);
	MatDestroy(&ISN2_x_part1_matrix333333_to_100000);	MatDestroy(&ISN2_y_part1_matrix333333_to_001000);	MatDestroy(&ISN2_z_part1_matrix333333_to_000010);

	MatDestroy(&ISN0_x_part2_matrix333333_to_100000);	MatDestroy(&ISN0_y_part2_matrix333333_to_001000);	MatDestroy(&ISN0_z_part2_matrix333333_to_000010);
	MatDestroy(&ISN1_x_part2_matrix333333_to_100000);	MatDestroy(&ISN1_y_part2_matrix333333_to_001000);	MatDestroy(&ISN1_z_part2_matrix333333_to_000010);
	MatDestroy(&ISN2_x_part2_matrix333333_to_100000);	MatDestroy(&ISN2_y_part2_matrix333333_to_001000);	MatDestroy(&ISN2_z_part2_matrix333333_to_000010);

	//导数相关-----------------------------------
	MatDestroy(&dx_matrix100000_or_010000_to_000000);
	MatDestroy(&dy_matrix001000_or_000100_to_000000);
	MatDestroy(&dz_matrix000010_or_000001_to_000000);
	MatDestroy(&dx_matrix111111_to_000000);		MatDestroy(&dy_matrix111111_to_000000);		MatDestroy(&dz_matrix111111_to_000000);
	MatDestroy(&ddx_matrix111111_to_000000);	MatDestroy(&ddy_matrix111111_to_000000);	MatDestroy(&ddz_matrix111111_to_000000);

	//平均数据-----------------------------------
	MatDestroy(&average_matrix100000_or_010000_to_000000);
	MatDestroy(&average_matrix001000_or_000100_to_000000);
	MatDestroy(&average_matrix000010_or_000001_to_000000);

	MatDestroy(&average_matrix344333_or_433433_to_333333);
	MatDestroy(&average_matrix343343_or_433334_to_333333);
	MatDestroy(&average_matrix333443_or_334334_to_333333);
}
