
#ifndef MODELCREATOR_FLEX_H
#define MODELCREATOR_FLEX_H

#include <BTCSet.h>
#include <Vector.h>

class System;

struct model_parameters_flex
{

};

class ModelCreator_Flex
{
public:
    ModelCreator_Flex();
    bool Create_Flex(System *system);

    // ASM1 Model Properties
    const double v_mu_H=1.742;
    const double v_K_S=8.065;
    const double v_K_MH=0.254;
    const double v_K_OH=0.030;
    const double v_K_NOH=0.067;
    const double v_eta_g=0.582;
    const double v_b_H=0.641;
    const double v_K_NH=0.014;
    const double v_mu_M=0.720;
    const double v_K_MM=0.119;
    const double v_K_OM=0.033;
    const double v_K_NOM=0.629;
    const double v_b_M=0.098;
    const double v_mu_A=1.059;
    const double v_K_NHA=0.952;
    const double v_K_NOA=0.020;
    const double v_K_OA=0.238;
    const double v_b_A=0.185;
    const double v_eta_h=0.736;
    const double v_K_h=1.088;
    const double v_K_X=0.095;
    const double v_K_a=0.057;
    const double v_Y_H=0.611;
    const double v_Y_HM=0.400;
    const double v_Y_A=0.240;
    const double v_Y_M=0.462;
    const double v_f_p=0.080;
    const double v_i_XB=0.058;
    const double v_i_VSSB=1.420;
    const double v_i_VSSi=2.000;
    const double v_i_VSSs=1.800;
    const double v_i_VSSP=1.420;
    const double v_i_MeOH=1.500;
    const double v_k_LO2=190.446; //Aearation coeff

    // Constituents Concentration
    const double v_S_i_concentration=14.14;
    const double v_S_S_concentration=0.5;
    const double v_S_O_concentration=0.1;
    const double v_S_NH_concentration=1;
    const double v_S_M_concentration=0;
    const double v_S_NO_concentration=1;
    const double v_S_ND_concentration=1;

    const double v_X_BM_concentration=1450; // ?
    const double v_X_BH_concentration=460; // ?
    const double v_X_BA_concentration=180; // ?
    const double v_X_S_concentration=9;
    const double v_X_p_concentration=350; // ?
    const double v_X_ND_concentration=1;

    const double v_a_rate_coefficient=v_k_LO2; // 5 ~ 10 1/hr : 120 ~ 240 1/day
    const double v_a_saturation=8.55;

    const double p_31=0.57;
    const double p_32=1.42;
    const double p_2=0.08;
    const double p_36=0.05;
    const double p_3=0.075;

    // OHQ Model properties
    const double v_settling_vel=10000; // X_b : Unit: m/day
    const double v_r_storage=17500; // Reactor: Intial Storage : Unit: m3

    const double v_r_volume=v_r_storage; // Reactor: Intial Volume : Unit: m3

    const double v_s_t_storage=200; // Settling element top: initial storage
    const double v_s_t_bottom_elevation=1; // Settling element top: bottom elevation
    const double v_s_b_storage=200; // Settling element bottom: initial storage
    const double v_s_b_bottom_elevation=0; // Settling element bottom: bottom elevation
    const double v_st_sb_area=1e6; // Link: Settling element top to Settling element bottom: area

    const double v_c_bottom_elevation=1; // Clarifer: bottom elevation (head)
    const double v_was_bottom_elevation=0; // WAS: bottom elevation(head)

    const double v_flow_factor_i=2e10; // Flow factor for storage-base flow -- inside
    const double v_flow_factor_o=2e9; // Flow factor for storage-base flow -- out of


    // Water Balance (Constant)
    const double v_r_constant_flow=800; // Reactor: Constant flow : Unit: m3/day
    const double v_r_st_constant_flow=1700; // Link: Reactor to Settling element top: flow
    const double v_st_c_constant_flow=750; // Link: Settling element top to Clarifier: flow
    const double v_st_sb_constant_flow=950; // Link: Settling element top to Settling element bottom: flow
    const double v_sb_r_constant_flow=900; // Link: Settling element bottom to Reactor: flow
    const double v_sb_was_constant_flow=50; // Link: Settling element bottom to WAS: flow

    int n_tanks = 8;
    vector<bool> aeration_v = {true,true,true,false,false,false,true,false}; // (n_tanks)

    const double v_t_storage=v_r_storage/n_tanks; // Tank: Intial Storage : Unit: m3
    const double v_t_constant_flow=v_r_constant_flow; // Tank: Constant flow : Unit: m3/day
    const double v_t_t_st_constant_flow=v_r_st_constant_flow; // Link: Tank to Tank or Tank to Settling element top: flow
    const double v_sb_t_constant_flow=v_sb_r_constant_flow; // Link: Settling element bottom to Tank: flow

    const double v_t_volume=v_t_storage; // Tank: Intial Volume : Unit: m3
    const double v_s_t_volume=v_s_t_storage; // Settling element top: initial volume : Unit: m3
    const double v_s_b_volume=v_s_b_storage; // Settling element bottom: initial volume : Unit: m3

    // Temperature effect
    const double v_rt=0; // Ref Temp (+20 as a room temp included in the Temp file)
    const double v_rtt=0; // Ref Temp (+20 as a room temp included in the Temp file)

    const double v_mu_H_af=1.072;
    const double v_K_S_af=1.03;
    const double v_mu_M_af=1.09;
    const double v_b_M_af=1.03;
    const double v_mu_A_af=1.072;
    const double v_b_A_af=1.03;
    const double v_K_h_af=1.03;
    const double v_K_a_af=1.03;


    // Calibration
    const double v_k_LO2_high=250; //Aearation coeff UL
    const double v_k_LO2_low=50; //Aearation coeff LL

private:
    const double pi = 3.141521;
};

#endif // MODELCREATOR_FLEX_H
