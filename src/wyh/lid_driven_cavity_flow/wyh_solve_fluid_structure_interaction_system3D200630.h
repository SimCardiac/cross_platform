#ifndef wyh_solve_fluid_structure_interaction_system3D200630H
#define wyh_solve_fluid_structure_interaction_system3D200630H

#include <petsc.h>
#include <mpi.h>

class wyh_solve_fluid_structure_interaction_system3D200630
{
public:

	//***********************************************************************************************************************************************
	int root = 0;
	double pi = 3.141592653;
	int judge_num;

	PetscInt istart;
	PetscInt iend;

	//与实际问题相关变量*****************************************************************************************************************************
	
	//参数---------------------------------------
	double rho = 1.0;
	double mu_f = 1e-2;
	double mu_s = 0.5;

	double r = 0.2;
	double p_center[3] = { 0.6, 0.5, 0.5 };//+

	//边界条件-----------------------------------
	int velocity_boundary = 3;
	int pressure_boundary = 2;
	int phi_boundary = 0;
	double hh_phi = 0.0;

	//空间网格-----------------------------------
	int m_num = 128;
	int n_num = 128;
	int o_num = 128;

	//边界坐标-----------------------------------
	double x_L = 0.0;	double x_U = 1.0;
	double y_L = 0.0;	double y_U = 1.0;
	double z_L = 0.0;	double z_U = 1.0;

	//空间步长-----------------------------------
	double hx = (x_U - x_L) / m_num;
	double hy = (y_U - y_L) / n_num;
	double hz = (z_U - z_L) / o_num;

	double epsilon = 2 * hx;

	//时间步长-----------------------------------
	double dt = 0.125*hx;

	//系数---------------------------------------
	Mat NS_coeff_u3D;//+
	Mat NS_coeff_v3D;//+
	Mat NS_coeff_w3D;//+
	Mat NS_coeff_p3D;//+

	//修正变量-----------------------------------
	Vec u_NS_right_term_take_out_x_L;
	Vec v_NS_right_term_take_out_y_L;
	Vec w_NS_right_term_take_out_z_L;

	Vec velocity_u_star_take_out_x_L;
	Vec velocity_v_star_take_out_y_L;
	Vec velocity_w_star_take_out_z_L;

	//修改右端项
	Vec modified_u_NSeq_vector;//+

	//调整矩阵-----------------------------------

	//去掉边界的矩阵
	Mat take_out_x_L_boundary_matrix3D;//+
	Mat take_out_y_L_boundary_matrix3D;//+
	Mat take_out_z_L_boundary_matrix3D;//+

	//添加边界的矩阵
	Mat x_solution_out_x_L_extend_to_ori_size_matrix3D;//+
	Mat y_solution_out_y_L_extend_to_ori_size_matrix3D;//+
	Mat z_solution_out_z_L_extend_to_ori_size_matrix3D;//+

	//程序中全局变量*********************************************************************************************************************************

	//拓展规模=====================================================================================

	//空间拓展数组-------------------------------
	int extend_num100000[6] = { 1, 0, 0, 0, 0, 0 };
	int extend_num010000[6] = { 0, 1, 0, 0, 0, 0 };
	int extend_num001000[6] = { 0, 0, 1, 0, 0, 0 };
	int extend_num000100[6] = { 0, 0, 0, 1, 0, 0 };
	int extend_num000010[6] = { 0, 0, 0, 0, 1, 0 };
	int extend_num000001[6] = { 0, 0, 0, 0, 0, 1 };
	int extend_num111111[6] = { 1, 1, 1, 1, 1, 1 };
	int extend_num333333[6] = { 3, 3, 3, 3, 3, 3 };
	int extend_num344333[6] = { 3, 4, 4, 3, 3, 3 };
	int extend_num343343[6] = { 3, 4, 3, 3, 4, 3 };
	int extend_num433433[6] = { 4, 3, 3, 4, 3, 3 };
	int extend_num333443[6] = { 3, 3, 3, 4, 4, 3 };
	int extend_num433334[6] = { 4, 3, 3, 3, 3, 4 };
	int extend_num334334[6] = { 3, 3, 4, 3, 3, 4 };

	//空间拓展数值-------------------------------
	int all_mesh_num_take_out_x_L = (m_num - 1)*n_num*o_num;
	int all_mesh_num_take_out_y_L = m_num*(n_num - 1)*o_num;
	int all_mesh_num_take_out_z_L = m_num*n_num*(o_num - 1);

	int all_mesh_num = m_num*n_num*o_num;
	int all_mesh_num100 = (m_num + 1)*(n_num + 0)*(o_num + 0);
	int all_mesh_num010 = (m_num + 0)*(n_num + 1)*(o_num + 0);
	int all_mesh_num001 = (m_num + 0)*(n_num + 0)*(o_num + 1);
	int all_mesh_num222 = (m_num + 2)*(n_num + 2)*(o_num + 2);
	int all_mesh_num666 = (m_num + 6)*(n_num + 6)*(o_num + 6);
	int all_mesh_num776 = (m_num + 7)*(n_num + 7)*(o_num + 6);
	int all_mesh_num767 = (m_num + 7)*(n_num + 6)*(o_num + 7);
	int all_mesh_num677 = (m_num + 6)*(n_num + 7)*(o_num + 7);

	//向量变量=====================================================================================

	//基本变量-----------------------------------

	//速度...................
	Vec velocity_u;			Vec velocity_v;			Vec velocity_w;//+
	Vec velocity_u_at_p;	Vec velocity_v_at_p;	Vec velocity_w_at_p;//+

	//压力...................
	Vec pressure;//+

	//phi....................
	Vec phi;//+

	//FF.....................
	Vec FF00;	Vec FF01;	Vec FF02;//+
	Vec FF10;	Vec FF11;	Vec FF12;//+
	Vec FF20;	Vec FF21;	Vec FF22;//+

	//FF_inv.................
	Vec FF_inv00;	Vec FF_inv01;	Vec FF_inv02;//+
	Vec FF_inv10;	Vec FF_inv11;	Vec FF_inv12;//+
	Vec FF_inv20;	Vec FF_inv21;	Vec FF_inv22;//+

	//应力张量...............
	Vec tensor00, tensor01, tensor02;//+
	Vec tensor10, tensor11, tensor12;//+
	Vec tensor20, tensor21, tensor22;//+

	//应力...................
	Vec F_x_at_p;//+
	Vec F_y_at_p;//+
	Vec F_z_at_p;//+

	Vec tensor_body_x_at_p;//+
	Vec tensor_body_y_at_p;//+
	Vec tensor_body_z_at_p;//+

	Vec tensor_normal_x_at_p;//+
	Vec tensor_normal_y_at_p;//+
	Vec tensor_normal_z_at_p;//+

	//方程相关-----------------------------------

	//过渡项.................
	Vec velocity_u_star;//+
	Vec velocity_v_star;//+
	Vec velocity_w_star;//+
	Vec modified_term_at_p;

	//NS方程右端项...........
	Vec u_NS_right_term;//+
	Vec v_NS_right_term;//+
	Vec w_NS_right_term;//+
	Vec modified_equation_right_term;//+

	//应力...................
	Vec F_x_at_u;//+
	Vec F_y_at_v;//+
	Vec F_z_at_w;//+

	//压力...................
	Vec d_pressure_x_at_u;//+
	Vec d_pressure_y_at_v;//+
	Vec d_pressure_z_at_w;//+

	//对流项.................
	Vec u_NS_convection;//+
	Vec v_NS_convection;//+
	Vec w_NS_convection;//+

	Vec u_NS_convection_ori;
	Vec v_NS_convection_ori;
	Vec w_NS_convection_ori;

	//扩散项.................
	Vec u_NS_diffusion;//+
	Vec v_NS_diffusion;//+
	Vec w_NS_diffusion;//+

	//法线向量-----------------------------------
	Vec phi_normal_x;//+
	Vec phi_normal_y;//+
	Vec phi_normal_z;//+

	Vec normal_norm;//+
	Vec square_phi_normal_x;//+
	Vec square_phi_normal_y;//+

	//区域划分-----------------------------------

	//内部区域
	Vec extract_interior_domain_vector;//+

	//边界区域
	Vec extract_boundary_domain_vector;// 
	
	//修正phi------------------------------------
	Vec ones_vector;
	Vec one_redu_phi;

	Vec phi_boundary_mark_vector;
	Vec sign_boundary_mark;//+

	Vec modified_right_term;
	Vec modified_part1;
	Vec modified_part2;
	Vec modified_part_x;
	Vec modified_part_y;
	Vec modified_part_z;

	Vec d_phi_x;
	Vec d_phi_y;
	Vec d_phi_z;
	Vec d_phi_norm;

	//导数---------------------------------------

	//WENO相关
	Vec d_uphi_x, d_vphi_y;//+

	Vec d_uu_x, d_vu_y;//+
	Vec d_uv_x, d_vv_y;//+
	Vec d_uw_x, d_vw_y;//+

	//中间速度1阶导数
	Vec d_u_x_at_p; Vec d_u_y_at_p; Vec d_u_z_at_p;//+
	Vec d_v_x_at_p; Vec d_v_y_at_p; Vec d_v_z_at_p;//+
	Vec d_w_x_at_p; Vec d_w_y_at_p; Vec d_w_z_at_p;//+

	//过渡速度1阶导数
	Vec d_u_star_x;//+
	Vec d_v_star_y;//+

	//速度二阶导数
	Vec dd_u_xx, dd_u_yy;//+
	Vec dd_v_xx, dd_v_yy;//+
	Vec dd_w_xx, dd_w_yy;//+

	//应力...................
	Vec d_tensor00_x, d_tensor01_x, d_tensor02_x;//+
	Vec d_tensor10_y, d_tensor11_y, d_tensor12_y;//+
	Vec d_tensor20_z, d_tensor21_z, d_tensor22_z;//+

	Vec d_modified_term_x;//+
	Vec d_modified_term_y;//+
	Vec d_modified_term_z;//+

	//修正压力...............
	Vec dd_modified_term_xx;//+
	Vec dd_modified_term_yy;//+
	Vec lap_modified_term;//+

	//过渡项-------------------------------------
	Vec right_term_part0;//+
	Vec right_term_part1;//+
	Vec right_term_part1_0;//+
	Vec right_term_part1_1;//+

	Vec mult_part01, mult_part02;//+
	Vec tensor00_transition, tensor01_transition, tensor02_transition;//+
	Vec tensor10_transition, tensor11_transition, tensor12_transition;//+
	Vec tensor20_transition, tensor21_transition, tensor22_transition;//+

	//拓展变量-----------------------------------

	//速度
	Vec velocity_u_star_at_u010000;//+
	Vec velocity_v_star_at_v000100;//+
	Vec velocity_w_star_at_w000001;//+

	Vec velocity_u_at_u010000;//+
	Vec velocity_v_at_v000100;//+
	Vec velocity_w_at_w000001;//+

	Vec velocity_u_at_u111111;//+
	Vec velocity_v_at_v111111;//+
	Vec velocity_w_at_w111111;//+

	Vec velocity_u_at_u333333, velocity_u_at_v333333, velocity_u_at_w333333;//+
	Vec velocity_v_at_u333333, velocity_v_at_v333333, velocity_v_at_w333333;//+
	Vec velocity_w_at_u333333, velocity_w_at_v333333, velocity_w_at_w333333;//+

	Vec velocity_u_at_u344333, velocity_u_at_u343343;//+
	Vec velocity_v_at_v433433, velocity_v_at_v333443;//+
	Vec velocity_w_at_w433334, velocity_w_at_w334334;//+

	Vec velocity_u_at_p111111, velocity_v_at_p111111, velocity_w_at_p111111;//+
	Vec velocity_u_at_p333333, velocity_v_at_p333333, velocity_w_at_p333333;//+

	//压力
	Vec pressure_at_p100000;//+
	Vec pressure_at_p001000;//+
	Vec pressure_at_p000010;//+

	//phi
	Vec phi_at_p111111;//+

	//应力
	Vec F_x_at_p100000, F_y_at_p001000, F_z_at_p000010;//+

	Vec tensor00_111111, tensor01_111111, tensor02_111111;//+
	Vec tensor10_111111, tensor11_111111, tensor12_111111;//+
	Vec tensor20_111111, tensor21_111111, tensor22_111111;//+

	//修正项
	Vec modified_term_at_p100000;//+
	Vec modified_term_at_p001000;//+
	Vec modified_term_at_p000010;//+
	Vec modified_term_at_p111111;//+

	//WENO-------------------

	//phi
	Vec variable333333;//+
	Vec f_variable333333;//+
	Vec uphi_f_P333333, uphi_f_N333333;//+
	Vec vphi_f_P333333, vphi_f_N333333;//+
	Vec wphi_f_P333333, wphi_f_N333333;//+

	//速度
	Vec uu_f_P, uu_f_N;//+
	Vec vu_f_P, vu_f_N;//+
	Vec wu_f_P, wu_f_N;//+

	Vec uv_f_P, uv_f_N;//+
	Vec vv_f_P, vv_f_N;//+
	Vec wv_f_P, wv_f_N;//+

	Vec uw_f_P, uw_f_N;//+
	Vec vw_f_P, vw_f_N;//+
	Vec ww_f_P, ww_f_N;//+

	//WENO---------------------------------------

	//速度绝对值的最大值.....
	Vec variable_abs;//+

	//WENO_at_x..............

	//过渡量
	Vec ISP0_part1_value_at_x, ISP1_part1_value_at_x, ISP2_part1_value_at_x;//+
	Vec ISP0_part2_value_at_x, ISP1_part2_value_at_x, ISP2_part2_value_at_x;//+
	Vec ISN0_part1_value_at_x, ISN1_part1_value_at_x, ISN2_part1_value_at_x;//+
	Vec ISN0_part2_value_at_x, ISN1_part2_value_at_x, ISN2_part2_value_at_x;//+
	Vec sum_alphaP_at_x;//+
	Vec sum_alphaN_at_x;//+
	Vec alphaP0_at_x, alphaP1_at_x, alphaP2_at_x;//+
	Vec alphaN0_at_x, alphaN1_at_x, alphaN2_at_x;//+

	//分量、权重和导数
	Vec d_variable_x_temp;//+
	Vec OmegaP0_at_x, OmegaP1_at_x, OmegaP2_at_x;//+
	Vec OmegaN0_at_x, OmegaN1_at_x, OmegaN2_at_x;//+
	Vec f_P0_at_x, f_P1_at_x, f_P2_at_x;//+
	Vec f_N0_at_x, f_N1_at_x, f_N2_at_x;//+
	Vec FP_at_x, FP_at_x0, FP_at_x1, FP_at_x2;//+
	Vec FN_at_x, FN_at_x0, FN_at_x1, FN_at_x2;//+

	//WENO_at_y..............

	//过渡量
	Vec ISP0_part1_value_at_y, ISP1_part1_value_at_y, ISP2_part1_value_at_y;//+
	Vec ISP0_part2_value_at_y, ISP1_part2_value_at_y, ISP2_part2_value_at_y;//+
	Vec ISN0_part1_value_at_y, ISN1_part1_value_at_y, ISN2_part1_value_at_y;//+
	Vec ISN0_part2_value_at_y, ISN1_part2_value_at_y, ISN2_part2_value_at_y;//+
	Vec sum_alphaP_at_y;//+
	Vec sum_alphaN_at_y;//+
	Vec alphaP0_at_y, alphaP1_at_y, alphaP2_at_y;//+
	Vec alphaN0_at_y, alphaN1_at_y, alphaN2_at_y;//+

	//分量、权重和导数
	Vec d_variable_y_temp;//+
	Vec OmegaP0_at_y, OmegaP1_at_y, OmegaP2_at_y;//+
	Vec OmegaN0_at_y, OmegaN1_at_y, OmegaN2_at_y;//+
	Vec f_P0_at_y, f_P1_at_y, f_P2_at_y;//+
	Vec f_N0_at_y, f_N1_at_y, f_N2_at_y;//+
	Vec FP_at_y, FP_at_y0, FP_at_y1, FP_at_y2;//+
	Vec FN_at_y, FN_at_y0, FN_at_y1, FN_at_y2;//+

	//WENO_at_z..............

	//过渡量
	Vec ISP0_part1_value_at_z, ISP1_part1_value_at_z, ISP2_part1_value_at_z;//+
	Vec ISP0_part2_value_at_z, ISP1_part2_value_at_z, ISP2_part2_value_at_z;//+
	Vec ISN0_part1_value_at_z, ISN1_part1_value_at_z, ISN2_part1_value_at_z;//+
	Vec ISN0_part2_value_at_z, ISN1_part2_value_at_z, ISN2_part2_value_at_z;//+
	Vec sum_alphaP_at_z;//+
	Vec sum_alphaN_at_z;//+
	Vec alphaP0_at_z, alphaP1_at_z, alphaP2_at_z;//+
	Vec alphaN0_at_z, alphaN1_at_z, alphaN2_at_z;//+

	//分量、权重和导数
	Vec d_variable_z_temp;//+
	Vec OmegaP0_at_z, OmegaP1_at_z, OmegaP2_at_z;//+
	Vec OmegaN0_at_z, OmegaN1_at_z, OmegaN2_at_z;//+
	Vec f_P0_at_z, f_P1_at_z, f_P2_at_z;//+
	Vec f_N0_at_z, f_N1_at_z, f_N2_at_z;//+
	Vec FP_at_z, FP_at_z0, FP_at_z1, FP_at_z2;//+
	Vec FN_at_z, FN_at_z0, FN_at_z1, FN_at_z2;//+

	//TVD变量------------------------------------

	//phi....................
	Vec TVD_phi1;//+
	Vec TVD_phi2;//+
	Vec TVD_phi_LSeq_right_term0;//+
	Vec TVD_phi_LSeq_right_term1;//+
	Vec TVD_phi_LSeq_right_term2;//+

	//FF.....................
	Vec TVD_FF00_1;	Vec TVD_FF00_2;//+
	Vec TVD_FF01_1;	Vec TVD_FF01_2;//+
	Vec TVD_FF02_1;	Vec TVD_FF02_2;//+
	Vec TVD_FF10_1;	Vec TVD_FF10_2;//+
	Vec TVD_FF11_1;	Vec TVD_FF11_2;//+
	Vec TVD_FF12_1;	Vec TVD_FF12_2;//+
	Vec TVD_FF20_1;	Vec TVD_FF20_2;//+
	Vec TVD_FF21_1;	Vec TVD_FF21_2;//+
	Vec TVD_FF22_1;	Vec TVD_FF22_2;//+

	Vec TVD_FF00_right_term0;	Vec TVD_FF00_right_term1;	Vec TVD_FF00_right_term2;//+
	Vec TVD_FF01_right_term0;	Vec TVD_FF01_right_term1;	Vec TVD_FF01_right_term2;//+
	Vec TVD_FF02_right_term0;	Vec TVD_FF02_right_term1;	Vec TVD_FF02_right_term2;//+
	Vec TVD_FF10_right_term0;	Vec TVD_FF10_right_term1;	Vec TVD_FF10_right_term2;//+
	Vec TVD_FF11_right_term0;	Vec TVD_FF11_right_term1;	Vec TVD_FF11_right_term2;//+
	Vec TVD_FF12_right_term0;	Vec TVD_FF12_right_term1;	Vec TVD_FF12_right_term2;//+
	Vec TVD_FF20_right_term0;	Vec TVD_FF20_right_term1;	Vec TVD_FF20_right_term2;//+
	Vec TVD_FF21_right_term0;	Vec TVD_FF21_right_term1;	Vec TVD_FF21_right_term2;//+
	Vec TVD_FF22_right_term0;	Vec TVD_FF22_right_term1;	Vec TVD_FF22_right_term2;//+

	//FF_inv.....................
	Vec TVD_FF_inv00_1;	Vec TVD_FF_inv00_2;//+
	Vec TVD_FF_inv01_1;	Vec TVD_FF_inv01_2;//+
	Vec TVD_FF_inv02_1;	Vec TVD_FF_inv02_2;//+
	Vec TVD_FF_inv10_1;	Vec TVD_FF_inv10_2;//+
	Vec TVD_FF_inv11_1;	Vec TVD_FF_inv11_2;//+
	Vec TVD_FF_inv12_1;	Vec TVD_FF_inv12_2;//+
	Vec TVD_FF_inv20_1;	Vec TVD_FF_inv20_2;//+
	Vec TVD_FF_inv21_1;	Vec TVD_FF_inv21_2;//+
	Vec TVD_FF_inv22_1;	Vec TVD_FF_inv22_2;//+

	Vec TVD_FF_inv00_right_term0;	Vec TVD_FF_inv00_right_term1;	Vec TVD_FF_inv00_right_term2;//+
	Vec TVD_FF_inv01_right_term0;	Vec TVD_FF_inv01_right_term1;	Vec TVD_FF_inv01_right_term2;//+
	Vec TVD_FF_inv02_right_term0;	Vec TVD_FF_inv02_right_term1;	Vec TVD_FF_inv02_right_term2;//+
	Vec TVD_FF_inv10_right_term0;	Vec TVD_FF_inv10_right_term1;	Vec TVD_FF_inv10_right_term2;//+
	Vec TVD_FF_inv11_right_term0;	Vec TVD_FF_inv11_right_term1;	Vec TVD_FF_inv11_right_term2;//+
	Vec TVD_FF_inv12_right_term0;	Vec TVD_FF_inv12_right_term1;	Vec TVD_FF_inv12_right_term2;//+
	Vec TVD_FF_inv20_right_term0;	Vec TVD_FF_inv20_right_term1;	Vec TVD_FF_inv20_right_term2;//+
	Vec TVD_FF_inv21_right_term0;	Vec TVD_FF_inv21_right_term1;	Vec TVD_FF_inv21_right_term2;//+
	Vec TVD_FF_inv22_right_term0;	Vec TVD_FF_inv22_right_term1;	Vec TVD_FF_inv22_right_term2;//+

	//KSP===========================================================================================
	KSP KSP_NSeq_u;//+
	KSP KSP_NSeq_v;//+
	KSP KSP_NSeq_w;//+
	KSP KSP_NSeq_p;//+

	//矩阵变量=====================================================================================

	//导数相关矩阵-------------------------------
	Mat dx_matrix100000_or_010000_to_000000;//+
	Mat dy_matrix001000_or_000100_to_000000;//+
	Mat dz_matrix000010_or_000001_to_000000;//+

	Mat dx_matrix111111_to_000000;//+
	Mat dy_matrix111111_to_000000;//+
	Mat dz_matrix111111_to_000000;//+ 

	Mat ddx_matrix111111_to_000000;//+
	Mat ddy_matrix111111_to_000000;//+
	Mat ddz_matrix111111_to_000000;//+

	//平均数据矩阵-------------------------------
	Mat average_matrix100000_or_010000_to_000000;//+
	Mat average_matrix001000_or_000100_to_000000;//+
	Mat average_matrix000010_or_000001_to_000000;//+

	Mat average_matrix344333_or_433433_to_333333;//+
	Mat average_matrix343343_or_433334_to_333333;//+
	Mat average_matrix333443_or_334334_to_333333;//+
	
	//拓展矩阵-----------------------------------

	//速度
	Mat nonstagger_extend_u_matrix111111;//+
	Mat nonstagger_extend_v_matrix111111;//+
	Mat nonstagger_extend_w_matrix111111;//+

	Mat nonstagger_extend_u_matrix333333;//+
	Mat nonstagger_extend_v_matrix333333;//+
	Mat nonstagger_extend_w_matrix333333;//+

	Mat extend_u_matrix010000;		Vec u_auxiliary_vector010000;	//+
	Mat extend_v_matrix000100;		//Vec v_auxiliary_vector000100; //+
	Mat extend_w_matrix000001;		//Vec w_auxiliary_vector000001; //+

	Mat extend_u_matrix111111;		Vec u_auxiliary_vector111111;	//+
	Mat extend_v_matrix111111;		//Vec v_auxiliary_vector111111; //+
	Mat extend_w_matrix111111;		//Vec w_auxiliary_vector111111;	//+

	Mat extend_u_matrix333333;		Vec u_auxiliary_vector333333;	//+
	Mat extend_v_matrix333333;		//Vec v_auxiliary_vector333333; //+
	Mat extend_w_matrix333333;		//Vec w_auxiliary_vector333333; //+

	Mat extend_u_matrix344333;		Vec u_auxiliary_vector344333;	//+
	Mat extend_u_matrix343343;		Vec u_auxiliary_vector343343;	//+

	Mat extend_v_matrix433433;		//Vec v_auxiliary_vector433433;	//+
	Mat extend_v_matrix333443;		//Vec v_auxiliary_vector333443; //+

	Mat extend_w_matrix433334;		//Vec w_auxiliary_vector433334;	//+
	Mat extend_w_matrix334334;		//Vec w_auxiliary_vector334334; //+

	//压力
	Mat extend_pressure_matrix100000;		//Vec pressure_auxiliary_vector100000; //+
	Mat extend_pressure_matrix001000;		//Vec pressure_auxiliary_vector001000; //+
	Mat extend_pressure_matrix000010;		//Vec pressure_auxiliary_vector000010; //+

	Mat extend_pressure_matrix111111;		//Vec pressure_auxiliary_vector111111; //+

	//phi
	Mat extend_phi_matrix100000;		//Vec phi_auxiliary_vector100000; //+
	Mat extend_phi_matrix001000;		//Vec phi_auxiliary_vector001000; //+
	Mat extend_phi_matrix000010;		//Vec phi_auxiliary_vector000010; //+

	Mat extend_phi_matrix111111;		//Vec phi_auxiliary_vector111111; //+
	Mat extend_phi_matrix333333;		//Vec phi_auxiliary_vector333333; //+

	//WENO相关矩阵-------------------------------

	//x
	Mat f_P0_x_matrix333333_to_100000;	Mat f_P1_x_matrix333333_to_100000;	Mat f_P2_x_matrix333333_to_100000;//+
	Mat f_N0_x_matrix333333_to_100000;	Mat f_N1_x_matrix333333_to_100000;	Mat f_N2_x_matrix333333_to_100000;//+
	Mat ISP0_x_part1_matrix333333_to_100000;	Mat ISP1_x_part1_matrix333333_to_100000;	Mat ISP2_x_part1_matrix333333_to_100000;//+
	Mat ISP0_x_part2_matrix333333_to_100000;	Mat ISP1_x_part2_matrix333333_to_100000;	Mat ISP2_x_part2_matrix333333_to_100000;//+
	Mat ISN0_x_part1_matrix333333_to_100000;	Mat ISN1_x_part1_matrix333333_to_100000;	Mat ISN2_x_part1_matrix333333_to_100000;//+
	Mat ISN0_x_part2_matrix333333_to_100000;	Mat ISN1_x_part2_matrix333333_to_100000;	Mat ISN2_x_part2_matrix333333_to_100000;//+

	//y
	Mat f_P0_y_matrix333333_to_001000;	Mat f_P1_y_matrix333333_to_001000;	Mat f_P2_y_matrix333333_to_001000;//+
	Mat f_N0_y_matrix333333_to_001000;	Mat f_N1_y_matrix333333_to_001000;	Mat f_N2_y_matrix333333_to_001000;//+
	Mat ISP0_y_part1_matrix333333_to_001000; 	Mat ISP1_y_part1_matrix333333_to_001000;	Mat ISP2_y_part1_matrix333333_to_001000;//+
	Mat ISP0_y_part2_matrix333333_to_001000;	Mat ISP1_y_part2_matrix333333_to_001000;	Mat ISP2_y_part2_matrix333333_to_001000;//+
	Mat ISN0_y_part1_matrix333333_to_001000;	Mat ISN1_y_part1_matrix333333_to_001000;	Mat ISN2_y_part1_matrix333333_to_001000;//+
	Mat ISN0_y_part2_matrix333333_to_001000;	Mat ISN1_y_part2_matrix333333_to_001000;	Mat ISN2_y_part2_matrix333333_to_001000;//+

	//z
	Mat f_P0_z_matrix333333_to_000010;	Mat f_P1_z_matrix333333_to_000010;	Mat f_P2_z_matrix333333_to_000010;//+
	Mat f_N0_z_matrix333333_to_000010;	Mat f_N1_z_matrix333333_to_000010;	Mat f_N2_z_matrix333333_to_000010;//+
	Mat ISP0_z_part1_matrix333333_to_000010;	Mat ISP1_z_part1_matrix333333_to_000010;	Mat ISP2_z_part1_matrix333333_to_000010;//+
	Mat ISP0_z_part2_matrix333333_to_000010;	Mat ISP1_z_part2_matrix333333_to_000010;	Mat ISP2_z_part2_matrix333333_to_000010;//+
	Mat ISN0_z_part1_matrix333333_to_000010;	Mat ISN1_z_part1_matrix333333_to_000010;	Mat ISN2_z_part1_matrix333333_to_000010;//+
	Mat ISN0_z_part2_matrix333333_to_000010;	Mat ISN1_z_part2_matrix333333_to_000010;	Mat ISN2_z_part2_matrix333333_to_000010;//+
		
	//函数*****************************************************************************************

	//设定变量大小并赋值------------------------------

	//变量大小--------------------
	void set_variable_size(MPI_Comm PETSC_COMM_WORLD);

	//变量赋值--------------------
	void prepare_assembly_matrix_and_vector();

	//实际问题的计算----------------------------------

	//计算-----------------------

	//更新速度与压力
	void solve_fluid_structure_interaction_system3D();

	//本构模型计算应力
	void get_NS_elastic_force3D();

	//初始值
	void get_initial_condition();

	//方程组系数------------------

	//流体速度系数
	void assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_u3D(Mat NS_coeff_u3D);

	void assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_v3D(Mat NS_coeff_v3D);

	void assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_w3D(Mat NS_coeff_w3D);

	//流体压力系数
	void assembly_lid_driven_cavity_flow_boundary_condition_NS_coeff_p3D(Mat NS_coeff_p3D);

	//实际问题的向量与矩阵--------

	//处理NS方程右端项
	void assembly_modified_u_NSeq_vector_for_lid_driven_cavity_flow_boundary_condition3D();

	//去掉边界数据
	void assembly_take_out_x_L_boundary_matrix3D();

	void assembly_take_out_y_L_boundary_matrix3D();

	void assembly_take_out_z_L_boundary_matrix3D();

	//添加边界数据
	void assembly_x_solution_out_x_L_extend_to_ori_size_matrix3D();

	void assembly_y_solution_out_y_L_extend_to_ori_size_matrix3D();

	void assembly_z_solution_out_z_L_extend_to_ori_size_matrix3D();

	//提取边界区域的矩阵
	void assembly_extract_boundary_domain_vector3D(Vec extract_boundary_domain_vector);

	//提取内部区域的矩阵
	void assembly_extract_interior_domain_vector3D(Vec extract_interior_domain_vector);

	//求解N-S方程-------------------------------------

	//修正速度及压力项------------

	//修正速度
	void modified_velocity3D();

	//修正压力
	void modified_pressure3D();

	//求解NS方程所涉及的右端项-----

	//修正方程的右端项
	void get_modified_equation_right_term3D();

	//NS右端项
	void get_fluid_solid_NSeq_right_term3D();

	//压力的梯度
	void get_d_pressure_at_stagger3D();

	//扩散项
	void get_NS_diffusion3D();

	//对流项
	void get_NS_convection3D();

	//更新FF与phi------------------------------------

	//转换速度--------------------

	//获得中心速度
	void get_velocity_at_p3D();

	//速度导数
	void get_d_velocity_at_p3D();

	//求解level输运方程-----------

	//phi输运方程右端项
	void get_variable_LSeq_right_term3D(Vec variable, Vec variable_LSeq_right_term);

	//FF输运方程右端项
	void get_FF_LSeq_right_term3D(
		Vec variable00, Vec variable01, Vec variable02,
		Vec variable10, Vec variable11, Vec variable12,
		Vec variable20, Vec variable21, Vec variable22,
		Vec variable_right_term00, Vec variable_right_term01, Vec variable_right_term02,
		Vec variable_right_term10, Vec variable_right_term11, Vec variable_right_term12,
		Vec variable_right_term20, Vec variable_right_term21, Vec variable_right_term22);

	void get_FF_inv_LSeq_right_term3D(
		Vec variable00, Vec variable01, Vec variable02,
		Vec variable10, Vec variable11, Vec variable12,
		Vec variable20, Vec variable21, Vec variable22,
		Vec variable_right_term00, Vec variable_right_term01, Vec variable_right_term02,
		Vec variable_right_term10, Vec variable_right_term11, Vec variable_right_term12,
		Vec variable_right_term20, Vec variable_right_term21, Vec variable_right_term22);

	//更新phi
	void update_phi3D();

	//更新FF
	void update_FF3D();

	//更新FF
	void update_FF_inv3D();

	//修正phi
	void modified_phi_based_on_intermediate_step();

	//内法线
	void solve_phi_normal();

	//计算WENO---------------------------------------

	//速度绝对值的最大值
	double get_max_abs(Vec variable);

	//WENO at x

	//一阶导数
	void get_WENO_d_variable_x_value(Vec variable_f_P, Vec variable_f_N, Vec d_variable_x);

	//分量
	void get_CellBoundaryFP_x(Vec variable_f_P);
	void get_CellBoundaryFN_x(Vec variable_f_N);

	//权重
	void get_WeightFP_at_x(Vec variable_f_P);
	void get_WeightFN_at_x(Vec variable_f_N);

	//WENO at y

	//一阶导数
	void get_WENO_d_variable_y_value(Vec variable_f_P, Vec variable_f_N, Vec d_variable_y);

	//分量
	void get_CellBoundaryFP_y(Vec variable_f_P);
	void get_CellBoundaryFN_y(Vec variable_f_N);

	//权重
	void get_WeightFP_at_y(Vec variable_f_P);
	void get_WeightFN_at_y(Vec variable_f_N);

	//WENO at z

	//一阶导数
	void get_WENO_d_variable_z_value(Vec variable_f_P, Vec variable_f_N, Vec d_variable_z);

	//分量
	void get_CellBoundaryFP_z(Vec variable_f_P);
	void get_CellBoundaryFN_z(Vec variable_f_N);

	//权重
	void get_WeightFP_at_z(Vec variable_f_P);
	void get_WeightFN_at_z(Vec variable_f_N);

	//WENO相关计算矩阵--------------------------------

	//分量-----------------------

	//Cell boundary FP
	void assembly_CellBoundaryFP_x_matrix333333_to_100000();
	void assembly_CellBoundaryFP_y_matrix333333_to_001000();
	void assembly_CellBoundaryFP_z_matrix333333_to_000010();

	//Cell boundary FN
	void assembly_CellBoundaryFN_x_matrix333333_to_100000();
	void assembly_CellBoundaryFN_y_matrix333333_to_001000();
	void assembly_CellBoundaryFN_z_matrix333333_to_000010();

	//权重-----------------------

	//part 1 of the Weight FP
	void assembly_part1_WeightFP_x_matrix333333_to_100000();
	void assembly_part1_WeightFP_y_matrix333333_to_001000();
	void assembly_part1_WeightFP_z_matrix333333_to_000010();

	//part 2 of the Weight FP
	void assembly_part2_WeightFP_x_matrix333333_to_100000();
	void assembly_part2_WeightFP_y_matrix333333_to_001000();
	void assembly_part2_WeightFP_z_matrix333333_to_000010();

	//part 1 of the Weight FN
	void assembly_part1_WeightFN_x_matrix333333_to_100000();
	void assembly_part1_WeightFN_y_matrix333333_to_001000();
	void assembly_part1_WeightFN_z_matrix333333_to_000010();

	//part 2 of the Weight FN
	void assembly_part2_WeightFN_x_matrix333333_to_100000();
	void assembly_part2_WeightFN_y_matrix333333_to_001000();
	void assembly_part2_WeightFN_z_matrix333333_to_000010();

	//求导相关矩阵------------------------------------

	//1阶差分--------------------

	//向前或向后差分
	void assembly_dx_matrix100000_or_010000_to_000000();
	void assembly_dy_matrix001000_or_000100_to_000000();
	void assembly_dz_matrix000010_or_000001_to_000000();

	//中心差分
	void assembly_dx_matrix111111_to_000000();
	void assembly_dy_matrix111111_to_000000();
	void assembly_dz_matrix111111_to_000000();

	//2阶差分--------------------
	void assembly_ddx_matrix111111_to_000000();
	void assembly_ddy_matrix111111_to_000000();
	void assembly_ddz_matrix111111_to_000000();

	//平均数据----------------------------------------

	//将(100000或010000类型数据求两点平均000000
	void assembly_average_matrix100000_or_010000_to_000000();

	//将(001000或000100类型数据求两点平均000000
	void assembly_average_matrix001000_or_000100_to_000000();

	//将(000010或000001类型数据求两点平均000000
	void assembly_average_matrix000010_or_000001_to_000000();

	//将(344333或433433类型数据求四点平均333333)
	void assembly_average_matrix344333_or_433433_to_333333();

	//将(343343或433334类型数据求四点平均333333)
	void assembly_average_matrix343343_or_433334_to_333333();

	//将(333443或334334类型数据求四点平均333333)
	void assembly_average_matrix333443_or_334334_to_333333();

	//拓展数据点--------------------------------------

	//拓展矩阵--------------------

	//非交错网格
	void assembly_nonstagger_extend_velocity_u_matrix3D(Mat extend_u_matrix, int *extend_num, int velocity_boundary);
	void assembly_nonstagger_extend_velocity_v_matrix3D(Mat extend_v_matrix, int *extend_num, int velocity_boundary);
	void assembly_nonstagger_extend_velocity_w_matrix3D(Mat extend_w_matrix, int *extend_num, int velocity_boundary);
	void assembly_nonstagger_extend_pressure_matrix3D(Mat extend_pressure_matrix, int *extend_num, int pressure_boundary);
	void assembly_nonstagger_extend_phi_matrix3D(Mat extend_phi_matrix, int *extend_num, int phi_boundary);

	//交错网格
	void assembly_stagger_extend_velocity_u_matrix3D(Mat extend_u_matrix, int *extend_num, int velocity_boundary);
	void assembly_stagger_extend_velocity_v_matrix3D(Mat extend_v_matrix, int *extend_num, int velocity_boundary);
	void assembly_stagger_extend_velocity_w_matrix3D(Mat extend_w_matrix, int *extend_num, int velocity_boundary);

	//拓展向量
	void assembly_extend_velocity_u_auxiliary_vector3D(Vec u_auxiliary_vector, int *extend_num, int velocity_boundary);
	void assembly_extend_velocity_v_auxiliary_vector3D(Vec v_auxiliary_vector, int *extend_num, int velocity_boundary);
	void assembly_extend_velocity_w_auxiliary_vector3D(Vec w_auxiliary_vector, int *extend_num, int velocity_boundary);
	void assembly_extend_pressure_auxiliary_vector3D(Vec pressure_auxiliary_vector, int *extend_num, int pressure_boundary);
	void assembly_extend_phi_auxiliary_vector3D(Vec phi_auxiliary_vector, int *extend_num, int phi_boundary);

	//Destroy----------------------------------------
	void destroy_variable();
};
#endif
