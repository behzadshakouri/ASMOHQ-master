#include "modelcreator_flex.h"
#include "System.h"
#include "QString"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

ModelCreator_Flex::ModelCreator_Flex()
{
}


bool ModelCreator_Flex::Create_Flex(System *system)
{
    system->GetQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/main_components.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/unsaturated_soil.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Well.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/wastewater.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/mass_transfer.json");
    system->ReadSystemSettingsTemplate("/home/behzad/Projects/OpenHydroQual/resources/settings.json");

    bool OUP=true; // True for using OUProcess, False for using DeNite Data

    bool St=false; // True for using Simulation Time is Days, False for using Start and End Date

    const double Simulation_time=1; // Simulation Time in Days
    const double dt = 0.01; // Time Step

    const double Simulation_start_time=40210; // Simulation Start Date
    const double Simulation_end_time=40359; // Simulation End Date

    const double Simulation_time_Calc = Simulation_end_time - Simulation_start_time;

    // Model Configuration

    // Consistuents
    Constituent S_S;
    S_S.SetQuantities(system, "Constituent");
    S_S.SetName("S_S");
    S_S.SetType("Constituent");
    system->AddConstituent(S_S,false);

    Constituent S_O;
    S_O.SetQuantities(system, "Constituent");
    S_O.SetName("S_O");
    S_O.SetType("Constituent");
    system->AddConstituent(S_O,false);

    Constituent S_NH;
    S_NH.SetQuantities(system, "Constituent");
    S_NH.SetName("S_NH");
    S_NH.SetType("Constituent");
    system->AddConstituent(S_NH,false);

    Constituent S_M;
    S_M.SetQuantities(system, "Constituent");
    S_M.SetName("S_M");
    S_M.SetType("Constituent");
    system->AddConstituent(S_M,false);

    Constituent S_NO;
    S_NO.SetQuantities(system, "Constituent");
    S_NO.SetName("S_NO");
    S_NO.SetType("Constituent");
    system->AddConstituent(S_NO,false);

    Constituent S_ND;
    S_ND.SetQuantities(system, "Constituent");
    S_ND.SetName("S_ND");
    S_ND.SetType("Constituent");
    system->AddConstituent(S_ND,false);

    Constituent X_BM;
    X_BM.SetQuantities(system, "Particle");
    X_BM.SetName("X_BM");
    X_BM.SetType("Particle");
    X_BM.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_BM,false);

    Constituent X_BH;
    X_BH.SetQuantities(system, "Particle");
    X_BH.SetName("X_BH");
    X_BH.SetType("Particle");
    X_BH.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_BH,false);

    Constituent X_BA;
    X_BA.SetQuantities(system, "Particle");
    X_BA.SetName("X_BA");
    X_BA.SetType("Particle");
    X_BA.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_BA,false);

    Constituent X_S;
    X_S.SetQuantities(system, "Particle");
    X_S.SetName("X_S");
    X_S.SetType("Particle");
    X_S.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_S,false);

    Constituent X_p;
    X_p.SetQuantities(system, "Particle");
    X_p.SetName("X_p");
    X_p.SetType("Particle");
    X_p.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_p,false);

    Constituent X_ND;
    X_ND.SetQuantities(system, "Particle");
    X_ND.SetName("X_ND");
    X_ND.SetType("Particle");
    X_ND.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_ND,false);


    // Reaction Parameters
    RxnParameter mu_H;
    mu_H.SetQuantities(system,"ReactionParameter");
    mu_H.SetName("mu_H");
    mu_H.SetVal("base_value",v_mu_H);
    system->AddReactionParameter(mu_H, false);

    RxnParameter K_S;
    K_S.SetQuantities(system,"ReactionParameter");
    K_S.SetName("K_S");
    K_S.SetVal("base_value",v_K_S);
    system->AddReactionParameter(K_S, false);

    RxnParameter K_MH;
    K_MH.SetQuantities(system,"ReactionParameter");
    K_MH.SetName("K_MH");
    K_MH.SetVal("base_value",v_K_MH);
    system->AddReactionParameter(K_MH, false);

    RxnParameter K_OH;
    K_OH.SetQuantities(system,"ReactionParameter");
    K_OH.SetName("K_OH");
    K_OH.SetVal("base_value",v_K_OH);
    system->AddReactionParameter(K_OH, false);

    RxnParameter K_NOH;
    K_NOH.SetQuantities(system,"ReactionParameter");
    K_NOH.SetName("K_NOH");
    K_NOH.SetVal("base_value",v_K_NOH);
    system->AddReactionParameter(K_NOH, false);

    RxnParameter eta_g;
    eta_g.SetQuantities(system,"ReactionParameter");
    eta_g.SetName("eta_g");
    eta_g.SetVal("base_value",v_eta_g);
    system->AddReactionParameter(eta_g, false);

    RxnParameter b_H;
    b_H.SetQuantities(system,"ReactionParameter");
    b_H.SetName("b_H");
    b_H.SetVal("base_value",v_b_H);
    system->AddReactionParameter(b_H, false);

    RxnParameter K_NH;
    K_NH.SetQuantities(system,"ReactionParameter");
    K_NH.SetName("K_NH");
    K_NH.SetVal("base_value",v_K_NH);
    system->AddReactionParameter(K_NH, false);

    RxnParameter mu_M;
    mu_M.SetQuantities(system,"ReactionParameter");
    mu_M.SetName("mu_M");
    mu_M.SetVal("base_value",v_mu_M);
    system->AddReactionParameter(mu_M, false);

    RxnParameter K_MM;
    K_MM.SetQuantities(system,"ReactionParameter");
    K_MM.SetName("K_MM");
    K_MM.SetVal("base_value",v_K_MM);
    system->AddReactionParameter(K_MM, false);

    RxnParameter K_OM;
    K_OM.SetQuantities(system,"ReactionParameter");
    K_OM.SetName("K_OM");
    K_OM.SetVal("base_value",v_K_OM);
    system->AddReactionParameter(K_OM, false);

    RxnParameter K_NOM;
    K_NOM.SetQuantities(system,"ReactionParameter");
    K_NOM.SetName("K_NOM");
    K_NOM.SetVal("base_value",v_K_NOM);
    system->AddReactionParameter(K_NOM, false);

    RxnParameter b_M;
    b_M.SetQuantities(system,"ReactionParameter");
    b_M.SetName("b_M");
    b_M.SetVal("base_value",v_b_M);
    system->AddReactionParameter(b_M, false);

    RxnParameter mu_A;
    mu_A.SetQuantities(system,"ReactionParameter");
    mu_A.SetName("mu_A");
    mu_A.SetVal("base_value",v_mu_A);
    system->AddReactionParameter(mu_A, false);

    RxnParameter K_NHA;
    K_NHA.SetQuantities(system,"ReactionParameter");
    K_NHA.SetName("K_NHA");
    K_NHA.SetVal("base_value",v_K_NHA);
    system->AddReactionParameter(K_NHA, false);

    RxnParameter K_NOA;
    K_NOA.SetQuantities(system,"ReactionParameter");
    K_NOA.SetName("K_NOA");
    K_NOA.SetVal("base_value",v_K_NOA);
    system->AddReactionParameter(K_NOA, false);

    RxnParameter K_OA;
    K_OA.SetQuantities(system,"ReactionParameter");
    K_OA.SetName("K_OA");
    K_OA.SetVal("base_value",v_K_OA);
    system->AddReactionParameter(K_OA, false);

    RxnParameter b_A;
    b_A.SetQuantities(system,"ReactionParameter");
    b_A.SetName("b_A");
    b_A.SetVal("base_value",v_b_A);
    system->AddReactionParameter(b_A, false);

    RxnParameter eta_h;
    eta_h.SetQuantities(system,"ReactionParameter");
    eta_h.SetName("eta_h");
    eta_h.SetVal("base_value",v_eta_h);
    system->AddReactionParameter(eta_h, false);

    RxnParameter K_h;
    K_h.SetQuantities(system,"ReactionParameter");
    K_h.SetName("K_h");
    K_h.SetVal("base_value",v_K_h);
    system->AddReactionParameter(K_h, false);

    RxnParameter K_X;
    K_X.SetQuantities(system,"ReactionParameter");
    K_X.SetName("K_X");
    K_X.SetVal("base_value",v_K_X);
    system->AddReactionParameter(K_X, false);

    RxnParameter K_a;
    K_a.SetQuantities(system,"ReactionParameter");
    K_a.SetName("K_a");
    K_a.SetVal("base_value",v_K_a);
    system->AddReactionParameter(K_a, false);

    RxnParameter Y_H;
    Y_H.SetQuantities(system,"ReactionParameter");
    Y_H.SetName("Y_H");
    Y_H.SetVal("base_value",v_Y_H);
    system->AddReactionParameter(Y_H, false);

    RxnParameter Y_HM;
    Y_HM.SetQuantities(system,"ReactionParameter");
    Y_HM.SetName("Y_HM");
    Y_HM.SetVal("base_value",v_Y_HM);
    system->AddReactionParameter(Y_HM, false);

    RxnParameter Y_A;
    Y_A.SetQuantities(system,"ReactionParameter");
    Y_A.SetName("Y_A");
    Y_A.SetVal("base_value",v_Y_A);
    system->AddReactionParameter(Y_A, false);

    RxnParameter Y_M;
    Y_M.SetQuantities(system,"ReactionParameter");
    Y_M.SetName("Y_M");
    Y_M.SetVal("base_value",v_Y_M);
    system->AddReactionParameter(Y_M, false);

    RxnParameter f_p;
    f_p.SetQuantities(system,"ReactionParameter");
    f_p.SetName("f_p");
    f_p.SetVal("base_value",v_f_p);
    system->AddReactionParameter(f_p, false);

    RxnParameter i_XB;
    i_XB.SetQuantities(system,"ReactionParameter");
    i_XB.SetName("i_XB");
    i_XB.SetVal("base_value",v_i_XB);
    system->AddReactionParameter(i_XB, false);

    RxnParameter i_VSSB;
    i_VSSB.SetQuantities(system,"ReactionParameter");
    i_VSSB.SetName("i_VSSB");
    i_VSSB.SetVal("base_value",v_i_VSSB);
    system->AddReactionParameter(i_VSSB, false);

    RxnParameter i_VSSi;
    i_VSSi.SetQuantities(system,"ReactionParameter");
    i_VSSi.SetName("i_VSSi");
    i_VSSi.SetVal("base_value",v_i_VSSi);
    system->AddReactionParameter(i_VSSi, false);

    RxnParameter i_VSSs;
    i_VSSs.SetQuantities(system,"ReactionParameter");
    i_VSSs.SetName("i_VSSs");
    i_VSSs.SetVal("base_value",v_i_VSSs);
    system->AddReactionParameter(i_VSSs, false);

    RxnParameter i_VSSP;
    i_VSSP.SetQuantities(system,"ReactionParameter");
    i_VSSP.SetName("i_VSSP");
    i_VSSP.SetVal("base_value",v_i_VSSP);
    system->AddReactionParameter(i_VSSP, false);

    RxnParameter i_MeOH;
    i_MeOH.SetQuantities(system,"ReactionParameter");
    i_MeOH.SetName("i_MeOH");
    i_MeOH.SetVal("base_value",v_i_MeOH);
    system->AddReactionParameter(i_MeOH, false);

    RxnParameter k_LO2;
    k_LO2.SetQuantities(system,"ReactionParameter");
    k_LO2.SetName("k_LO2");
    k_LO2.SetVal("base_value",v_k_LO2);
    system->AddReactionParameter(k_LO2, false);


    // Sources
    Source Aeration;
    Aeration.SetQuantities(system, "atmospheric exchange");
    Aeration.SetType("atmospheric exchange");
    Aeration.SetName("Aeration");
    Aeration.SetVal("rate_coefficient",v_a_rate_coefficient);
    Aeration.SetVal("saturation",v_a_saturation);
    system->AddSource(Aeration, false);

    // Reactions
    Reaction AerobicGH; // Reaction 1
    AerobicGH.SetQuantities(system,"Reaction");
    AerobicGH.SetName("AerobicGH");
    AerobicGH.SetProperty("S_S:stoichiometric_constant","(0-1/Y_H)");
    AerobicGH.SetProperty("S_O:stoichiometric_constant","(0-(1-Y_H)/Y_H)");
    AerobicGH.SetProperty("S_NH:stoichiometric_constant","(0-i_XB)");
    AerobicGH.SetProperty("X_BH:stoichiometric_constant","(1)");
    AerobicGH.SetProperty("rate_expression","(mu_H*S_S/(K_S+S_S)*S_O/(K_OH+S_O)*S_NH/(K_NH+S_NH)*X_BH)");
    system->AddReaction(AerobicGH,false);

    Reaction AerobicGHM; // Reaction 2
    AerobicGHM.SetQuantities(system,"Reaction");
    AerobicGHM.SetName("AerobicGHM");
    AerobicGHM.SetProperty("S_M:stoichiometric_constant","(0-1/Y_HM)");
    AerobicGHM.SetProperty("S_O:stoichiometric_constant","(0-(1-Y_H)/Y_H)");
    AerobicGHM.SetProperty("S_NH:stoichiometric_constant","(0-i_XB)");
    AerobicGHM.SetProperty("X_BH:stoichiometric_constant","(1)");
    AerobicGHM.SetProperty("rate_expression","(mu_H*S_M/(K_MH+S_M)*S_O/(K_OH+S_O)*(S_NH/K_NH+S_NH)*X_BH)");
    system->AddReaction(AerobicGHM,false);

    Reaction AnoxicGH; // Reaction 3
    AnoxicGH.SetQuantities(system,"Reaction");
    AnoxicGH.SetName("AnoxicGH");
    AnoxicGH.SetProperty("S_S:stoichiometric_constant","(0-1/Y_H)");
    AnoxicGH.SetProperty("S_NO:stoichiometric_constant","(0-(1-Y_H)/(2.86*Y_H))");
    AnoxicGH.SetProperty("S_NH:stoichiometric_constant","(0-i_XB)");
    AnoxicGH.SetProperty("X_BH:stoichiometric_constant","(1)");
    AnoxicGH.SetProperty("rate_expression","(mu_H*S_S/(K_S+S_S)*K_OH/(K_OH+S_O)*S_NO/(K_NOH+S_NO)*S_NH/(K_NH+S_NH)*eta_g*X_BH)");
    system->AddReaction(AnoxicGH,false);

    Reaction AnoxicGM; // Reaction 4
    AnoxicGM.SetQuantities(system,"Reaction");
    AnoxicGM.SetName("AnoxicGM");
    AnoxicGM.SetProperty("S_M:stoichiometric_constant","(0-1/Y_M)");
    AnoxicGM.SetProperty("S_NO:stoichiometric_constant","(0-(1-Y_M)/(2.86*Y_M))");
    AnoxicGM.SetProperty("S_NH:stoichiometric_constant","(0-i_XB)");
    AnoxicGM.SetProperty("X_BM:stoichiometric_constant","(1)");
    AnoxicGM.SetProperty("rate_expression","(mu_M*S_M/(K_MM+S_M)*K_OM/(K_OM+S_O)*S_NO/(K_NOM+S_NO)*S_NH/(K_NH+S_NH)*X_BM)");
    system->AddReaction(AnoxicGM,false);

    Reaction AerobicGA; // Reaction 5
    AerobicGA.SetQuantities(system,"Reaction");
    AerobicGA.SetName("AerobicGA");
    AerobicGA.SetProperty("S_O:stoichiometric_constant","(0-(4.57-Y_A)/Y_A)");
    AerobicGA.SetProperty("S_NH:stoichiometric_constant","(0-(i_XB+1/Y_A))");
    AerobicGA.SetProperty("X_BA:stoichiometric_constant","(1)");
    AerobicGA.SetProperty("S_NO:stoichiometric_constant","(1/Y_A)");
    AerobicGA.SetProperty("rate_expression","(mu_A*S_O/(K_OA+S_O)*S_NH/(K_NHA+S_NH)*X_BA)");
    system->AddReaction(AerobicGA,false);

    Reaction DecayH; // Reaction 6
    DecayH.SetQuantities(system,"Reaction");
    DecayH.SetName("DecayH");
    DecayH.SetProperty("X_BH:stoichiometric_constant","(-1)");
    DecayH.SetProperty("X_S:stoichiometric_constant","(1-f_p)");
    DecayH.SetProperty("X_p:stoichiometric_constant","(f_p)");
    DecayH.SetProperty("X_ND:stoichiometric_constant","((1-f_p)*i_XB)");
    DecayH.SetProperty("rate_expression","(b_H*(S_O/(K_OH+S_O)+eta_h*S_NO/(K_NOH+S_NO)*K_OH/(K_OH+S_O))*X_BH)");
    system->AddReaction(DecayH,false);

    Reaction DecayM; // Reaction 7
    DecayM.SetQuantities(system,"Reaction");
    DecayM.SetName("DecayM");
    DecayM.SetProperty("X_BM:stoichiometric_constant","(-1)");
    DecayM.SetProperty("X_S:stoichiometric_constant","(1-f_p)");
    DecayM.SetProperty("X_p:stoichiometric_constant","(f_p)");
    DecayM.SetProperty("X_ND:stoichiometric_constant","((1-f_p)*i_XB)");
    DecayM.SetProperty("rate_expression","(b_M*(S_O/(K_OM+S_O)+eta_h*S_NO/(K_NOM+S_NO)*K_OM/(K_OM+S_O))*X_BM)");
    system->AddReaction(DecayM,false);

    Reaction DecayA; // Reaction 8
    DecayA.SetQuantities(system,"Reaction");
    DecayA.SetName("DecayA");
    DecayA.SetProperty("X_BA:stoichiometric_constant","(-1)");
    DecayA.SetProperty("X_S:stoichiometric_constant","(1-f_p)");
    DecayA.SetProperty("X_p:stoichiometric_constant","(f_p)");
    DecayA.SetProperty("X_ND:stoichiometric_constant","((1-f_p)*i_XB)");
    DecayA.SetProperty("rate_expression","(b_A*(S_O/(K_OA+S_O)+eta_h*S_NO/(K_NOA+S_NO)*K_OA/(K_OA+S_O))*X_BA)");
    system->AddReaction(DecayA,false);

    Reaction AmmonificationSON; // Reaction 9
    AmmonificationSON.SetQuantities(system,"Reaction");
    AmmonificationSON.SetName("AmmonificationSON");
    AmmonificationSON.SetProperty("S_ND:stoichiometric_constant","(-1)");
    AmmonificationSON.SetProperty("S_NH:stoichiometric_constant","(1)");
    AmmonificationSON.SetProperty("rate_expression","(K_a*S_ND*(X_BH+X_BM+X_BA))");
    system->AddReaction(AmmonificationSON,false);

    Reaction HydrolysisEO; // Reaction 10
    HydrolysisEO.SetQuantities(system,"Reaction");
    HydrolysisEO.SetName("HydrolysisEO");
    HydrolysisEO.SetProperty("X_S:stoichiometric_constant","(-1)");
    HydrolysisEO.SetProperty("S_S:stoichiometric_constant","(1)");
    HydrolysisEO.SetProperty("rate_expression","(K_h*(X_S/(X_BH+X_BM+X_BA))/(K_X+(X_S/(X_BH+X_BM+X_BA)))*((S_O/(K_OH+S_O))+(eta_h*(S_NO/(K_NOH+S_NO))*(K_OH/(K_OH+S_O))))*(X_BH+X_BM+X_BA))");
    system->AddReaction(HydrolysisEO,false);


    Reaction HydrolysisEON; // Reaction 11
    HydrolysisEON.SetQuantities(system,"Reaction");
    HydrolysisEON.SetName("HydrolysisEON");
    HydrolysisEON.SetProperty("X_ND:stoichiometric_constant","(-1)");
    HydrolysisEON.SetProperty("S_ND:stoichiometric_constant","(1)");
    HydrolysisEON.SetProperty("rate_expression","((X_ND/X_S)*K_h*(X_S/(X_BH+X_BM+X_BA))/(K_X+(X_S/(X_BH+X_BM+X_BA)))*((S_O/(K_OH+S_O))+(eta_h*(S_NO/(K_NOH+S_NO))*(K_OH/(K_OH+S_O))))*(X_BH+X_BM+X_BA))");
    system->AddReaction(HydrolysisEON,false);


    // Settling Elements
    Block Stl_element_top;
    Stl_element_top.SetQuantities(system, "Settling element");
    Stl_element_top.SetName("Settling element top");
    Stl_element_top.SetType("Settling element");
    Stl_element_top.SetVal("Coagulant:concentration",0);
    Stl_element_top.SetVal("Settled_Particles:concentration",0);
    Stl_element_top.SetVal("Solids:concentration",0);
    Stl_element_top.SetVal("Storage",v_s_t_storage);
    Stl_element_top.SetVal("bottom_elevation",v_s_t_bottom_elevation);
    Stl_element_top.SetVal("Volume",v_s_t_volume);
    Stl_element_top.SetVal("x",800);
    Stl_element_top.SetVal("y",600);
    system->AddBlock(Stl_element_top,false);

    Block Stl_element_bottom;
    Stl_element_bottom.SetQuantities(system, "Settling element");
    Stl_element_bottom.SetName("Settling element bottom");
    Stl_element_bottom.SetType("Settling element");
    Stl_element_bottom.SetVal("Coagulant:concentration",0);
    Stl_element_bottom.SetVal("Settled_Particles:concentration",0);
    Stl_element_bottom.SetVal("Solids:concentration",0);
    Stl_element_bottom.SetVal("Storage",v_s_b_storage);
    Stl_element_bottom.SetVal("bottom_elevation",v_s_b_bottom_elevation);
    Stl_element_bottom.SetVal("Volume",v_s_b_volume);
    Stl_element_bottom.SetVal("x",800);
    Stl_element_bottom.SetVal("y",1000);
    system->AddBlock(Stl_element_bottom,false);

/*
    // Fixed Head Blocks
    Block fh_clarifier;
    fh_clarifier.SetQuantities(system, "fixed_head");
    fh_clarifier.SetName("Clarifier");
    fh_clarifier.SetType("fixed_head");
    fh_clarifier.SetVal("head",v_c_bottom_elevation);
    fh_clarifier.SetVal("S_S:concentration",0);
    fh_clarifier.SetVal("X_b:concentration",0);
    fh_clarifier.SetVal("Storage",100000);
    fh_clarifier.SetVal("x",1200);
    fh_clarifier.SetVal("y",600);
    system->AddBlock(fh_clarifier,false);

    Block fh_WAS;
    fh_WAS.SetQuantities(system, "fixed_head");
    fh_WAS.SetName("WAS");
    fh_WAS.SetType("fixed_head");
    fh_WAS.SetVal("head",v_was_bottom_elevation);
    fh_WAS.SetVal("S_S:concentration",0);
    fh_WAS.SetVal("X_b:concentration",0);
    fh_WAS.SetVal("Storage",100000);
    fh_WAS.SetVal("x",1200);
    fh_WAS.SetVal("y",1000);
    system->AddBlock(fh_WAS,false);
*/

    // Time-variable Fixed Head Blocks
    Block fh_clarifier;
    fh_clarifier.SetQuantities(system, "time_variable_fixed_head");
    fh_clarifier.SetName("Clarifier");
    fh_clarifier.SetType("time_variable_fixed_head");
    fh_clarifier.SetVal("head",v_c_bottom_elevation);
    fh_clarifier.SetVal("S_S:concentration",0);
    fh_clarifier.SetVal("X_b:concentration",0);
    fh_clarifier.SetVal("Storage",100000);
    fh_clarifier.SetVal("x",1200);
    fh_clarifier.SetVal("y",600);
    system->AddBlock(fh_clarifier,false);

    Block fh_WAS;
    fh_WAS.SetQuantities(system, "time_variable_fixed_head");
    fh_WAS.SetName("WAS");
    fh_WAS.SetType("time_variable_fixed_head");
    fh_WAS.SetVal("head",v_was_bottom_elevation);
    fh_WAS.SetVal("S_S:concentration",0);
    fh_WAS.SetVal("X_b:concentration",0);
    fh_WAS.SetVal("Storage",100000);
    fh_WAS.SetVal("x",1200);
    fh_WAS.SetVal("y",1000);
    system->AddBlock(fh_WAS,false);


    // Reactor_Flex Block
    for (int i=0; i<n_tanks; i++)
    {
        Block Reactor_Flex;
        Reactor_Flex.SetQuantities(system, "Reactor_Flex");
        Reactor_Flex.SetName("Reactor_Flex(" + aquiutils::numbertostring(i+1)+")");
        Reactor_Flex.SetType("Reactor_Flex");

        Reactor_Flex.SetVal("S_S:concentration",v_S_S_concentration);
        Reactor_Flex.SetVal("S_O:concentration",v_S_O_concentration);
        Reactor_Flex.SetVal("S_M:concentration",v_S_M_concentration);
        Reactor_Flex.SetVal("S_NO:concentration",v_S_NO_concentration);
        Reactor_Flex.SetVal("S_NH:concentration",v_S_NH_concentration);
        Reactor_Flex.SetVal("S_ND:concentration",v_S_ND_concentration);

        Reactor_Flex.SetVal("X_S:concentration",v_X_S_concentration);
        Reactor_Flex.SetVal("X_p:concentration",v_X_p_concentration);
        Reactor_Flex.SetVal("X_BH:concentration",v_X_BH_concentration);
        Reactor_Flex.SetVal("X_BM:concentration",v_X_BM_concentration);
        Reactor_Flex.SetVal("X_BA:concentration",v_X_BA_concentration);
        Reactor_Flex.SetVal("X_ND:concentration",v_X_ND_concentration);

        Reactor_Flex.SetVal("Storage",v_t_storage);
        //Reactor_Flex.SetVal("Storage",v_r_storage);
        //if (i==0) Reactor_Flex.SetVal("constant_inflow",v_r_constant_flow);
        //Reactor_Flex.SetVal("constant_inflow",v_r_constant_flow);
        //if (i==4) Reactor_Flex.SetVal("S_M:External mass flow time-series","/home/behzad/Projects/ASM_Models/S_M_mfr.csv");
        Reactor_Flex.SetVal("Volume",v_t_volume);
        Reactor_Flex.SetVal("x",400*(i-n_tanks));
        if (i==0) Reactor_Flex.SetVal("y",1200);
        else if (i>0) Reactor_Flex.SetVal("y",800);
        system->AddBlock(Reactor_Flex,false);

        // Aeration to all tanks
        //system->block("Reactor_Flex(" + aquiutils::numbertostring(i+1)+")")->SetProperty("S_O:external_source","Aeration"); // Reactor_Flex does not have source "Aeration", so we have to call it!

        // Aeration to selected tanks by aeration_v
        if (aeration_v[i]) system->block("Reactor_Flex(" + aquiutils::numbertostring(i+1)+")")->SetProperty("S_O:external_source","Aeration"); // Reactor_Flex does not have source "Aeration", so we have to call it!
    }


    // Producing Constant Inflows

    CTimeSeries<double> S_S_inflow_concentration;
    S_S_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_S_concentration);
    S_S_inflow_concentration.writefile("/home/behzad/Projects/ASM_Models/S_S_constant_inflow_concentration.txt");

    CTimeSeries<double> X_S_inflow_concentration;
    X_S_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_X_S_concentration);
    X_S_inflow_concentration.writefile("/home/behzad/Projects/ASM_Models/X_S_constant_inflow_concentration.txt");

    CTimeSeries<double> X_p_inflow_concentration;
    X_p_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_X_p_concentration);
    X_p_inflow_concentration.writefile("/home/behzad/Projects/ASM_Models/X_p_constant_inflow_concentration.txt");

    CTimeSeries<double> S_NO_inflow_concentration;
    S_NO_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_NO_concentration);
    S_NO_inflow_concentration.writefile("/home/behzad/Projects/ASM_Models/S_NO_constant_inflow_concentration.txt");

    CTimeSeries<double> S_NH_inflow_concentration;
    S_NH_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_NH_concentration);
    S_NH_inflow_concentration.writefile("/home/behzad/Projects/ASM_Models/S_NH_constant_inflow_concentration.txt");

    CTimeSeries<double> S_ND_inflow_concentration;
    S_ND_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_ND_concentration);
    S_ND_inflow_concentration.writefile("/home/behzad/Projects/ASM_Models/S_ND_constant_inflow_concentration.txt");

    CTimeSeries<double> X_ND_inflow_concentration;
    X_ND_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_X_ND_concentration);
    X_ND_inflow_concentration.writefile("/home/behzad/Projects/ASM_Models/X_ND_constant_inflow_concentration.txt");


    // Constituents Inflow Calculations from Data (Time Variable)

#ifdef Behzad
    CTimeSeriesSet<double> Inflow_DeNit("/home/behzad/Projects/ASM_Models/Data/DeNit_Influent_Lump.txt",true); // Inflow (Q)
    CTimeSeriesSet<double> Inflow_DeNit_MeOH("/home/behzad/Projects/ASM_Models/Data/DeNit_MeOH.txt",true); // Methanol (Externtal S_M)
    CTimeSeriesSet<double> Inflow_DeNit_wasteflow("/home/behzad/Projects/ASM_Models/Data/DeNit_wasteflow.txt",true); // Wasteflow (WAS)
    CTimeSeriesSet<double> Inflow_DeNit_returnflow("/home/behzad/Projects/ASM_Models/Data/DeNit_returnflow.txt",true); // Returnflow (RAS)


#else
    CTimeSeriesSet<double> Inflow_DeNit("/home/arash/Projects/ASM_Models/Data/DeNit_Influent_Lump.txt",true); // Inflow (Q)
    CTimeSeriesSet<double> Inflow_DeNit_MeOH("/home/arash/Projects/ASM_Models/Data/DeNit_MeOH.txt",true); // Methanol (Externtal S_M)
    CTimeSeriesSet<double> Inflow_DeNit_wasteflow("/home/arash/Projects/ASM_Models/Data/DeNit_wasteflow.txt",true); // Wasteflow (WAS)
    CTimeSeriesSet<double> Inflow_DeNit_returnflow("/home/arash/Projects/ASM_Models/Data/DeNit_returnflow.txt",true); // Returnflow (RAS)

#endif


    if (!OUP)

    {
    // Determining Coefficients

    const double c_S_i=p_31; // S_i, 8 //--
    const double c_S_S=1-p_31; // S_S, 8
    const double c_X_S=0.71*p_32*(1-p_2); // X_S, 2
    const double c_X_p=0.71*p_32*p_2; // X_p, 2
    const double c_S_NO=1; // S_NO, 10
    const double c_S_NH=1; // S_NH, 9
    const double c_S_ND=p_36*(1-p_31); // S_ND, 8
    const double c_X_ND=p_3*p_32; // X_ND, 2

    // Inflows and Flows

    CTimeSeries<double> Inflow_Q; // Discharge (m3/day)

    CTimeSeries<double> Flow_WAS; // Wasteflow (WAS) (m3/day)
    CTimeSeries<double> Flow_RAS; // Returnflow (RAS) (m3/day)

    CTimeSeries<double> Flow_r_r_st; // Reactor_Flex to Reactor_Flex + Reactor_Flex to Settling element top flow (m3/day)
    CTimeSeries<double> Flow_st_sb; // Settling element top to Settling element bottom flow (m3/day)
    CTimeSeries<double> Flow_st_c; // Settling element top to Clarifier flow (m3/day)

    // Defining Constituents Inflows

    CTimeSeries<double> Inflow_S_i; //--
    CTimeSeries<double> Inflow_S_S;
    CTimeSeries<double> Inflow_X_S;
    CTimeSeries<double> Inflow_X_p;
    CTimeSeries<double> Inflow_S_NO;
    CTimeSeries<double> Inflow_S_NH;
    CTimeSeries<double> Inflow_S_ND;
    CTimeSeries<double> Inflow_X_ND;

    CTimeSeries<double> Inflow_MeOH; // Methanol

    // Calculating Inflows according DeNite data (Time Variable)

    Inflow_Q=Inflow_DeNit.BTC[0]; // Discharge (m3/day)

    Flow_WAS=Inflow_DeNit_wasteflow.BTC[0]; // Wasteflow (WAS) (m3/day)
    Flow_RAS=Inflow_DeNit_returnflow.BTC[0]; // Returnflow (RAS) (m3/day)

    Flow_r_r_st=Inflow_Q+Flow_RAS;
    Flow_st_sb=Flow_RAS+Flow_WAS;
    Flow_st_c=Inflow_Q-Flow_WAS;

    Inflow_S_i=c_S_i*Inflow_DeNit.BTC[8]; //--
    Inflow_S_S=c_S_S*Inflow_DeNit.BTC[8];
    Inflow_X_S=c_X_S*Inflow_DeNit.BTC[2];
    Inflow_X_p=c_X_p*Inflow_DeNit.BTC[2];
    Inflow_S_NO=c_S_NO*Inflow_DeNit.BTC[10];
    Inflow_S_NH=c_S_NH*Inflow_DeNit.BTC[9];
    Inflow_S_ND=c_S_ND*Inflow_DeNit.BTC[8];
    Inflow_X_ND=c_X_ND*Inflow_DeNit.BTC[2];

    Inflow_MeOH=Inflow_DeNit_MeOH.BTC[0];; // Methanol

    // Writing to File

    Inflow_Q.writefile("/home/behzad/Projects/ASM_Models/Q_tvif.txt"); //

    Flow_WAS.writefile("/home/behzad/Projects/ASM_Models/WAS_tvf.txt"); //
    Flow_RAS.writefile("/home/behzad/Projects/ASM_Models/RAS_tvf.txt"); //

    Flow_r_r_st.writefile("/home/behzad/Projects/ASM_Models/r_r_st_tvf.txt"); //
    Flow_st_sb.writefile("/home/behzad/Projects/ASM_Models/st_sb_tvf.txt"); //
    Flow_st_c.writefile("/home/behzad/Projects/ASM_Models/st_c_tvf.txt"); //

    Inflow_S_i.writefile("/home/behzad/Projects/ASM_Models/S_i_Inflow_concentration.txt");
    Inflow_S_S.writefile("/home/behzad/Projects/ASM_Models/S_S_Inflow_concentration.txt");
    Inflow_X_S.writefile("/home/behzad/Projects/ASM_Models/X_S_Inflow_concentration.txt");
    Inflow_X_p.writefile("/home/behzad/Projects/ASM_Models/X_p_Inflow_concentration.txt");
    Inflow_S_NO.writefile("/home/behzad/Projects/ASM_Models/S_NO_Inflow_concentration.txt");
    Inflow_S_NH.writefile("/home/behzad/Projects/ASM_Models/S_NH_Inflow_concentration.txt");
    Inflow_S_ND.writefile("/home/behzad/Projects/ASM_Models/S_ND_Inflow_concentration.txt");
    Inflow_X_ND.writefile("/home/behzad/Projects/ASM_Models/X_ND_Inflow_concentration.txt");

    Inflow_MeOH.writefile("/home/behzad/Projects/ASM_Models/S_M_mfr.txt"); // Methanol

    // Assigning Constituents Inflow concentrations to the system (Reactor_Flex)

    // Constant inflows by given concentration
/*
    system->block("Reactor_Flex")->SetProperty("S_S:constant_inflow_concentration","/home/behzad/Projects/ASM_Models/S_S_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("X_S:constant_inflow_concentration","/home/behzad/Projects/ASM_Models/X_S_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("X_p:constant_inflow_concentration","/home/behzad/Projects/ASM_Models/X_p_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("S_NO:constant_inflow_concentration","/home/behzad/Projects/ASM_Models/S_NO_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("S_NH:constant_inflow_concentration","/home/behzad/Projects/ASM_Models/S_NH_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("S_ND:constant_inflow_concentration","/home/behzad/Projects/ASM_Models/S_ND_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("X_ND:constant_inflow_concentration","/home/behzad/Projects/ASM_Models/X_ND_constant_inflow_concentration.txt");
*/


    // Time Variable inflows according DeNite data

    system->block("Reactor_Flex(1)")->SetProperty("time_variable_inflow","/home/behzad/Projects/ASM_Models/Q_tvif.txt"); // Discharge (m3/day)

    system->block("Reactor_Flex(1)")->SetProperty("S_S:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/S_S_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("X_S:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/X_S_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("X_p:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/X_p_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("S_NO:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/S_NO_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("S_NH:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/S_NH_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("S_ND:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/S_ND_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("X_ND:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/X_ND_Inflow_concentration.txt");

    system->block("Reactor_Flex(5)")->SetProperty("S_M:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/S_M_mfr.txt");

    // Flex_flow Links for Reactor_Flex
    for (int i=0; i<n_tanks-1; i++)
    {   Link l_r_st;
        l_r_st.SetQuantities(system, "Flex_flow");
        l_r_st.SetName("Reactor_Flex(" + aquiutils::numbertostring(i+1)+"_" + aquiutils::numbertostring(i+2) + ")");
        l_r_st.SetType("Flex_flow");
        //l_r_st.SetProperty("flow", "/home/behzad/Projects/ASM_Models/r_r_st_tvf.txt");
        l_r_st.SetVal("flow_factor",v_flow_factor);
        system->AddLink(l_r_st, "Reactor_Flex(" + aquiutils::numbertostring(i+1)+")", "Reactor_Flex(" + aquiutils::numbertostring(i+2)+")", false);
    }


    // Flex_flow Links
    Link l_r_st;
    l_r_st.SetQuantities(system, "Flex_flow");
    l_r_st.SetName("Reactor_Flex(" + aquiutils::numbertostring(n_tanks) + ") - Settling element top");
    l_r_st.SetType("Flex_flow");
    //l_r_st.SetProperty("flow", "/home/behzad/Projects/ASM_Models/r_r_st_tvf.txt");
    l_r_st.SetVal("flow_factor",v_flow_factor);
    system->AddLink(l_r_st, "Reactor_Flex(" + aquiutils::numbertostring(n_tanks) + ")", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Flex_flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Flex_flow");
    //l_st_c.SetProperty("flow", "/home/behzad/Projects/ASM_Models/st_c_tvf.txt");
    l_st_c.SetVal("flow_factor",v_flow_factor);
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_sb_was;
    l_sb_was.SetQuantities(system, "Flex_flow");
    l_sb_was.SetName("Settling element bottom - WAS");
    l_sb_was.SetType("Flex_flow");
    //l_sb_was.SetProperty("flow", "/home/behzad/Projects/ASM_Models/WAS_tvf.txt");
    l_sb_was.SetVal("flow_factor",v_flow_factor);
    system->AddLink(l_sb_was, "Settling element bottom", "WAS", false);

    // Time-Dependent interface Link
    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface (time variable)");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface (time variable)");
    l_st_sb.SetProperty("flow", "/home/behzad/Projects/ASM_Models/st_sb_tvf.txt");
    l_st_sb.SetVal("area", v_st_sb_area);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    // Time-Dependent flow Links
    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Time-Dependent flow");
    l_sb_r.SetName("Settling element bottom - Reactor_Flex(1)");
    l_sb_r.SetType("Time-Dependent flow");
    l_sb_r.SetProperty("flow", "/home/behzad/Projects/ASM_Models/RAS_tvf.txt");
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor_Flex(1)", false);

    }

    if (OUP)

    {
    // Determining Inflows by OUProcess

        // Inflows

        CTimeSeries<double> OUP_Inflow_Q; // Discharge (m3/day)

        CTimeSeries<double> OUP_Flow_WAS; // Wasteflow (WAS) (m3/day)
        CTimeSeries<double> OUP_Flow_RAS; // Returnflow (RAS) (m3/day)

        CTimeSeries<double> OUP_Flow_r_r_st; // Reactor_Flex to Reactor_Flex + Reactor_Flex to Settling element top flow (m3/day)
        CTimeSeries<double> OUP_Flow_st_sb; // Settling element top to Settling element bottom flow (m3/day)
        CTimeSeries<double> OUP_Flow_st_c; // Settling element top to Clarifier flow (m3/day)

        // Defining Constituents Inflows

        CTimeSeries<double> OUP_Inflow_S_i; //--
        CTimeSeries<double> OUP_Inflow_S_S;
        CTimeSeries<double> OUP_Inflow_X_S;
        CTimeSeries<double> OUP_Inflow_X_p;
        CTimeSeries<double> OUP_Inflow_S_NO;
        CTimeSeries<double> OUP_Inflow_S_NH;
        CTimeSeries<double> OUP_Inflow_S_ND;
        CTimeSeries<double> OUP_Inflow_X_ND;

        CTimeSeries<double> OUP_Inflow_MeOH; // Methanol

    CTimeSeries<double> OUP_Inflow_Q_NS; // Discharge (m3/day)
    OUP_Inflow_Q_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_Q_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_Q_NS_tvif.csv");
    vector<double> Q_params; Q_params.push_back(10); Q_params.push_back(1);
    OUP_Inflow_Q = OUP_Inflow_Q_NS.MapfromNormalScoreToDistribution("lognormal", Q_params);
    OUP_Inflow_Q.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_Q_tvif.csv");

    CTimeSeries<double> OUP_Flow_WAS_NS; // Wasteflow (WAS) (m3/day)
    OUP_Flow_WAS_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Flow_WAS_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Flow_WAS_NS_tvf.csv");
    vector<double> WAS_params; WAS_params.push_back(3); WAS_params.push_back(0.8);
    OUP_Flow_WAS = OUP_Flow_WAS_NS.MapfromNormalScoreToDistribution("lognormal", WAS_params);
    OUP_Flow_WAS.writefile("/home/behzad/Projects/ASM_Models/OUP_Flow_WAS_tvf.csv");

    CTimeSeries<double> OUP_Flow_RAS_NS; // Returnflow (RAS) (m3/day)
    OUP_Flow_RAS_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Flow_RAS_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Flow_RAS_NS_tvf.csv");
    vector<double> RAS_params; RAS_params.push_back(9.2); RAS_params.push_back(1);
    OUP_Flow_RAS = OUP_Flow_RAS_NS.MapfromNormalScoreToDistribution("lognormal", RAS_params);
    OUP_Flow_RAS.writefile("/home/behzad/Projects/ASM_Models/OUP_Flow_RAS_tvf.csv");


    OUP_Flow_r_r_st=OUP_Inflow_Q+OUP_Flow_RAS;
    OUP_Flow_st_sb=OUP_Flow_RAS+OUP_Flow_WAS;
    OUP_Flow_st_c=OUP_Inflow_Q-OUP_Flow_WAS;


    OUP_Flow_r_r_st.writefile("/home/behzad/Projects/ASM_Models/OUP_r_r_st_tvf.csv"); //
    OUP_Flow_st_sb.writefile("/home/behzad/Projects/ASM_Models/OUP_st_sb_tvf.csv"); //
    OUP_Flow_st_c.writefile("/home/behzad/Projects/ASM_Models/OUP_st_c_tvf.csv"); //


    CTimeSeries<double> OUP_Inflow_S_i_NS;  //--
    OUP_Inflow_S_i_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_S_i_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_i_NS.csv");
    vector<double> S_i_params; S_i_params.push_back(0.2); S_i_params.push_back(0.5);
    OUP_Inflow_S_i = OUP_Inflow_S_i_NS.MapfromNormalScoreToDistribution("lognormal", S_i_params);
    OUP_Inflow_S_i.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_i.csv");

    CTimeSeries<double> OUP_Inflow_S_S_NS;
    OUP_Inflow_S_S_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_S_S_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_S_NS.csv");
    vector<double> S_S_params; S_S_params.push_back(1); S_S_params.push_back(0.4);
    OUP_Inflow_S_S = OUP_Inflow_S_S_NS.MapfromNormalScoreToDistribution("lognormal", S_S_params);
    OUP_Inflow_S_S.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_S.csv");

    CTimeSeries<double> OUP_Inflow_X_S_NS;
    OUP_Inflow_X_S_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_X_S_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_X_S_NS.csv");
    vector<double> X_S_params; X_S_params.push_back(1.5); X_S_params.push_back(0.5);
    OUP_Inflow_X_S = OUP_Inflow_X_S_NS.MapfromNormalScoreToDistribution("lognormal", X_S_params);
    OUP_Inflow_X_S.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_X_S.csv");

    CTimeSeries<double> OUP_Inflow_X_p_NS;
    OUP_Inflow_X_p_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_X_p_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_X_p_NS.csv");
    vector<double> X_p_params; X_p_params.push_back(1.2); X_p_params.push_back(0.4);
    OUP_Inflow_X_p = OUP_Inflow_X_p_NS.MapfromNormalScoreToDistribution("lognormal", X_p_params);
    OUP_Inflow_X_p.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_X_p.csv");

    CTimeSeries<double> OUP_Inflow_S_NO_NS;
    OUP_Inflow_S_NO_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_S_NO_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_NO_NS.csv");
    vector<double> S_NO_params; S_NO_params.push_back(0.8); S_NO_params.push_back(0.4);
    OUP_Inflow_S_NO = OUP_Inflow_S_NO_NS.MapfromNormalScoreToDistribution("lognormal", S_NO_params);
    OUP_Inflow_S_NO.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_NO.csv");

    CTimeSeries<double> OUP_Inflow_S_NH_NS;
    OUP_Inflow_S_NH_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_S_NH_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_NH_NS.csv");
    vector<double> S_NH_params; S_NH_params.push_back(1.3); S_NH_params.push_back(0.5);
    OUP_Inflow_S_NH = OUP_Inflow_S_NH_NS.MapfromNormalScoreToDistribution("lognormal", S_NH_params);
    OUP_Inflow_S_NH.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_NH.csv");

    CTimeSeries<double> OUP_Inflow_S_ND_NS;
    OUP_Inflow_S_ND_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_S_ND_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_ND_NS.csv");
    vector<double> S_ND_params; S_ND_params.push_back(0.7); S_ND_params.push_back(0.3);
    OUP_Inflow_S_ND = OUP_Inflow_S_ND_NS.MapfromNormalScoreToDistribution("lognormal", S_ND_params);
    OUP_Inflow_S_ND.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_S_ND.csv");

    CTimeSeries<double> OUP_Inflow_X_ND_NS;
    OUP_Inflow_X_ND_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_X_ND_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_X_ND_NS.csv");
    vector<double> X_ND_params; X_ND_params.push_back(1.4); X_ND_params.push_back(0.3);
    OUP_Inflow_X_ND = OUP_Inflow_X_ND_NS.MapfromNormalScoreToDistribution("lognormal", X_ND_params);
    OUP_Inflow_X_ND.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_X_ND.csv");

    CTimeSeries<double> OUP_Inflow_MeOH_NS; // Methanol
    OUP_Inflow_MeOH_NS.CreateOUProcess(0,Simulation_time_Calc,dt,1);
    OUP_Inflow_MeOH_NS.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_MeOH_NS.csv");
    vector<double> MeOH_params; MeOH_params.push_back(3); MeOH_params.push_back(0.6);
    OUP_Inflow_MeOH = OUP_Inflow_MeOH_NS.MapfromNormalScoreToDistribution("lognormal", MeOH_params);
    OUP_Inflow_MeOH.writefile("/home/behzad/Projects/ASM_Models/OUP_Inflow_MeOH.csv");

    // Time Variable inflows by OUProcess

    system->block("Reactor_Flex(1)")->SetProperty("time_variable_inflow","/home/behzad/Projects/ASM_Models/OUP_Inflow_Q_tvif.csv"); // Discharge (m3/day)

    system->block("Reactor_Flex(1)")->SetProperty("S_S:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/OUP_Inflow_S_S.csv");
    system->block("Reactor_Flex(1)")->SetProperty("X_S:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/OUP_Inflow_X_S.csv");
    system->block("Reactor_Flex(1)")->SetProperty("X_p:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/OUP_Inflow_X_p.csv");
    system->block("Reactor_Flex(1)")->SetProperty("S_NO:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/OUP_Inflow_S_NO.csv");
    system->block("Reactor_Flex(1)")->SetProperty("S_NH:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/OUP_Inflow_S_NH.csv");
    system->block("Reactor_Flex(1)")->SetProperty("S_ND:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/OUP_Inflow_S_ND.csv");
    system->block("Reactor_Flex(1)")->SetProperty("X_ND:time_variable_inflow_concentration","/home/behzad/Projects/ASM_Models/OUP_Inflow_X_ND.csv");

    system->block("Reactor_Flex(5)")->SetProperty("S_M:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/OUP_Inflow_MeOH.csv");


    // Flex_flow Links for Reactor_Flex
    for (int i=0; i<n_tanks-1; i++)
    {   Link l_r_st;
        l_r_st.SetQuantities(system, "Flex_flow");
        l_r_st.SetName("Reactor_Flex(" + aquiutils::numbertostring(i+1)+"_" + aquiutils::numbertostring(i+2) + ")");
        l_r_st.SetType("Flex_flow");
        //l_r_st.SetProperty("flow", "/home/behzad/Projects/ASM_Models/OUP_r_r_st_tvf.csv");
        l_r_st.SetVal("flow_factor",v_flow_factor);
        system->AddLink(l_r_st, "Reactor_Flex(" + aquiutils::numbertostring(i+1)+")", "Reactor_Flex(" + aquiutils::numbertostring(i+2)+")", false);
    }


    // Flex_flow Links
    Link l_r_st;
    l_r_st.SetQuantities(system, "Flex_flow");
    l_r_st.SetName("Reactor_Flex(" + aquiutils::numbertostring(n_tanks) + ") - Settling element top");
    l_r_st.SetType("Flex_flow");
    //l_r_st.SetProperty("flow", "/home/behzad/Projects/ASM_Models/OUP_r_r_st_tvf.csv");
    l_r_st.SetVal("flow_factor",v_flow_factor);
    system->AddLink(l_r_st, "Reactor_Flex(" + aquiutils::numbertostring(n_tanks) + ")", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Flex_flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Flex_flow");
    //l_st_c.SetProperty("flow", "/home/behzad/Projects/ASM_Models/OUP_st_c_tvf.csv");
    l_st_c.SetVal("flow_factor",v_flow_factor);
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_sb_was;
    l_sb_was.SetQuantities(system, "Flex_flow");
    l_sb_was.SetName("Settling element bottom - WAS");
    l_sb_was.SetType("Flex_flow");
    //l_sb_was.SetProperty("flow", "/home/behzad/Projects/ASM_Models/OUP_Flow_WAS_NS_tvf.csv");
    l_sb_was.SetVal("flow_factor",v_flow_factor);
    system->AddLink(l_sb_was, "Settling element bottom", "WAS", false);

    // Time-Dependent interface Link
    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface (time variable)");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface (time variable)");
    l_st_sb.SetProperty("flow", "/home/behzad/Projects/ASM_Models/OUP_st_sb_tvf.csv");
    l_st_sb.SetVal("area", v_st_sb_area);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    // Time-Dependent flow Links
    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Time-Dependent flow");
    l_sb_r.SetName("Settling element bottom - Reactor_Flex(1)");
    l_sb_r.SetType("Time-Dependent flow");
    l_sb_r.SetProperty("flow", "/home/behzad/Projects/ASM_Models/OUP_Flow_RAS_tvf.csv");
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor_Flex(1)", false);

}

    // Observations
    Observation total_inflow;

    total_inflow.SetQuantities(system, "Observation");
    total_inflow.SetProperty("expression","inflow");
    total_inflow.SetProperty("object","Reactor_Flex(1)");
    total_inflow.SetName("Inflow");
    total_inflow.SetType("Observation");
    system->AddObservation(total_inflow,false);

    Observation WAS_flow;

    WAS_flow.SetQuantities(system, "Observation");
    WAS_flow.SetProperty("expression","flow");
    WAS_flow.SetProperty("object","Settling element bottom - WAS");
    WAS_flow.SetName("WAS_Flow");
    WAS_flow.SetType("Observation");
    system->AddObservation(WAS_flow,false);

    Observation RAS_flow;

    RAS_flow.SetQuantities(system, "Observation");
    RAS_flow.SetProperty("expression","flow");
    RAS_flow.SetProperty("object","Settling element bottom - Reactor_Flex(1)");
    RAS_flow.SetName("RAS_Flow");
    RAS_flow.SetType("Observation");
    system->AddObservation(RAS_flow,false);

    /*
    Observation S_i_inflow_cn; //--

    S_i_inflow_cn.SetQuantities(system, "Observation");
    S_i_inflow_cn.SetProperty("expression","S_i:time_variable_inflow_concentration");
    S_i_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    S_i_inflow_cn.SetName("S_i_InflowConcentration");
    S_i_inflow_cn.SetType("Observation");
    system->AddObservation(S_i_inflow_cn,false);
    */

    Observation S_S_inflow_cn;

    S_S_inflow_cn.SetQuantities(system, "Observation");
    S_S_inflow_cn.SetProperty("expression","S_S:time_variable_inflow_concentration");
    S_S_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    S_S_inflow_cn.SetName("S_S_InflowConcentration");
    S_S_inflow_cn.SetType("Observation");
    system->AddObservation(S_S_inflow_cn,false);

    Observation X_S_inflow_cn; //?

    X_S_inflow_cn.SetQuantities(system, "Observation");
    X_S_inflow_cn.SetProperty("expression","X_S:time_variable_inflow_concentration");
    X_S_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    X_S_inflow_cn.SetName("X_S_InflowConcentration");
    X_S_inflow_cn.SetType("Observation");
    system->AddObservation(X_S_inflow_cn,false);

    Observation X_p_inflow_cn; //?

    X_p_inflow_cn.SetQuantities(system, "Observation");
    X_p_inflow_cn.SetProperty("expression","X_p:time_variable_inflow_concentration");
    X_p_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    X_p_inflow_cn.SetName("X_p_InflowConcentration");
    X_p_inflow_cn.SetType("Observation");
    system->AddObservation(X_p_inflow_cn,false);

    Observation S_NO_inflow_cn;

    S_NO_inflow_cn.SetQuantities(system, "Observation");
    S_NO_inflow_cn.SetProperty("expression","S_NO:time_variable_inflow_concentration");
    S_NO_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    S_NO_inflow_cn.SetName("S_NO_InflowConcentration");
    S_NO_inflow_cn.SetType("Observation");
    system->AddObservation(S_NO_inflow_cn,false);

    Observation S_NH_inflow_cn;

    S_NH_inflow_cn.SetQuantities(system, "Observation");
    S_NH_inflow_cn.SetProperty("expression","S_NH:time_variable_inflow_concentration");
    S_NH_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    S_NH_inflow_cn.SetName("S_NH_InflowConcentration");
    S_NH_inflow_cn.SetType("Observation");
    system->AddObservation(S_NH_inflow_cn,false);

    Observation S_ND_inflow_cn;

    S_ND_inflow_cn.SetQuantities(system, "Observation");
    S_ND_inflow_cn.SetProperty("expression","S_ND:time_variable_inflow_concentration");
    S_ND_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    S_ND_inflow_cn.SetName("S_ND_InflowConcentration");
    S_ND_inflow_cn.SetType("Observation");
    system->AddObservation(S_ND_inflow_cn,false);

    Observation X_ND_inflow_cn; //?

    X_ND_inflow_cn.SetQuantities(system, "Observation");
    X_ND_inflow_cn.SetProperty("expression","X_ND:time_variable_inflow_concentration");
    X_ND_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    X_ND_inflow_cn.SetName("X_ND_InflowConcentration");
    X_ND_inflow_cn.SetType("Observation");
    system->AddObservation(X_ND_inflow_cn,false);

    Observation S_M_inflow_cn; // Methonal : Reactor 5

    S_M_inflow_cn.SetQuantities(system, "Observation");
    S_M_inflow_cn.SetProperty("expression","S_M:time_variable_inflow_concentration");
    S_M_inflow_cn.SetProperty("object","Reactor_Flex(5)");
    S_M_inflow_cn.SetName("S_M_InflowConcentration");
    S_M_inflow_cn.SetType("Observation");
    system->AddObservation(S_M_inflow_cn,false);

    Observation X_BH_cn;

    X_BH_cn.SetQuantities(system, "Observation");
    X_BH_cn.SetProperty("expression","X_BH:concentration");
    X_BH_cn.SetProperty("object","Reactor_Flex(1)");
    X_BH_cn.SetName("X_BH_Concentration");
    X_BH_cn.SetType("Observation");
    system->AddObservation(X_BH_cn,false);

    /*
    Observation coagulant_concentration;

    coagulant_concentration.SetQuantities(system, "Observation");
    coagulant_concentration.SetProperty("expression","Coagulant:external_mass_flow_timeseries");
    coagulant_concentration.SetProperty("object","Settling element (1)");
    coagulant_concentration.SetName("Coagulant_Concentration");
    coagulant_concentration.SetType("Observation");
    system->AddObservation(coagulant_concentration,false);

    Observation solids_concentration;

    solids_concentration.SetQuantities(system, "Observation");
    solids_concentration.SetProperty("expression","Solids:concentration");
    solids_concentration.SetProperty("object","Settling element (1)");
    solids_concentration.SetName("Solids_Concentration");
    solids_concentration.SetType("Observation");
    system->AddObservation(solids_concentration,false);
    */


    if (St)
    {
        system->SetSettingsParameter("simulation_end_time",Simulation_time);
    }

    else if (!St)
    {
    system->SetSettingsParameter("simulation_start_time",Simulation_start_time);
    system->SetSettingsParameter("simulation_end_time",Simulation_end_time);
    }

    system->SetSystemSettings();
    cout<<"Populate functions"<<endl;
    system->PopulateOperatorsFunctions();
    cout<<"Variable parents"<<endl;
    system->SetVariableParents();
    return true;
}







