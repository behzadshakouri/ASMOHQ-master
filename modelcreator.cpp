#include "modelcreator.h"
#include "System.h"
#include "QString"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

ModelCreator::ModelCreator()
{
}


bool ModelCreator::Create(System *system)
{
#ifdef Behzad 
    system->GetQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/main_components.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/unsaturated_soil.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Well.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/wastewater.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/mass_transfer.json");
    system->ReadSystemSettingsTemplate("/home/behzad/Projects/OpenHydroQual/resources/settings.json");
#endif // Behzad 
#ifdef Arash_Windows
    system->GetQuanTemplate("../../OpenHydroQual/resources/main_components.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/unsaturated_soil.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Well.json");
    system->AppendQuanTemplate("../../OpenHydroQual/resources/wastewater.json");
    system->AppendQuanTemplate("../../OpenHydroQual/resources/mass_transfer.json");
    system->ReadSystemSettingsTemplate("../../OpenHydroQual/resources/settings.json");
#endif // Window_Arash

    
    bool St=false; // True for using Simulation Time is Days, False for using Start and End Date

    const double Simulation_time=1; // Simulation Time in Days

    const double Simulation_start_time=40210; // Simulation Start Date
    const double Simulation_end_time=40359; // Simulation End Date

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
    AerobicGHM.SetProperty("rate_expression","(mu_H*S_M/(K_MH+S_M)*S_O/(K_OH+S_O)*S_NH/(K_NH+S_NH)*X_BH)");
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

    /*
    CTimeSeries<double> CoagNS;
    CoagNS.CreateOUProcess(0,Simulation_time,0.05,1);
    CoagNS.writefile(Workingfolder + "coagulant_mfr_NS.csv");
    vector<double> c_params; c_params.push_back(2.5); c_params.push_back(0.6);
    CTimeSeries<double> Coag = CoagNS.MapfromNormalScoreToDistribution("lognormal", c_params);
    //Reactor.Variable("Coagulant:external_mass_flow_timeseries")->SetTimeSeries(Coag);
    Coag.writefile(Workingfolder + "coagulant_mfr.csv");
    Stl_element_top.SetProperty("Coagulant:external_mass_flow_timeseries",Workingfolder + "coagulant_mfr.csv");
    */

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

    /*
    CTimeSeries<double> CoagNS2;
    CoagNS2.CreateOUProcess(0,Simulation_time,0.05,1);
    CoagNS2.writefile(Workingfolder + "coagulant_mfr_NS.csv");
    vector<double> c_params2; c_params2.push_back(2.5); c_params2.push_back(0.6);
    CTimeSeries<double> Coag2 = CoagNS2.MapfromNormalScoreToDistribution("lognormal", c_params);
    //Reactor.Variable("Coagulant:external_mass_flow_timeseries")->SetTimeSeries(Coag);
    Coag2.writefile(Workingfolder + "coagulant_mfr.csv");
    Stl_element_bottom.SetProperty("Coagulant:external_mass_flow_timeseries",Workingfolder + "coagulant_mfr.csv");
    */

    Stl_element_bottom.SetVal("Settled_Particles:concentration",0);
    Stl_element_bottom.SetVal("Solids:concentration",0);
    Stl_element_bottom.SetVal("Storage",v_s_b_storage);
    Stl_element_bottom.SetVal("bottom_elevation",v_s_b_bottom_elevation);
    Stl_element_top.SetVal("Volume",v_s_b_volume);
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


    // Reactor Block
    for (int i=0; i<n_tanks; i++)
    {
        Block Reactor;
        Reactor.SetQuantities(system, "Reactor");
        Reactor.SetName("Reactor(" + aquiutils::numbertostring(i+1)+")");
        Reactor.SetType("Reactor");

        Reactor.SetVal("S_S:concentration",v_S_S_concentration);
        Reactor.SetVal("S_O:concentration",v_S_O_concentration);
        Reactor.SetVal("S_M:concentration",v_S_M_concentration);
        Reactor.SetVal("S_NO:concentration",v_S_NO_concentration);
        Reactor.SetVal("S_NH:concentration",v_S_NH_concentration);
        Reactor.SetVal("S_ND:concentration",v_S_ND_concentration);

        Reactor.SetVal("X_S:concentration",v_X_S_concentration);
        Reactor.SetVal("X_p:concentration",v_X_p_concentration);
        Reactor.SetVal("X_BH:concentration",v_X_BH_concentration);
        Reactor.SetVal("X_BM:concentration",v_X_BM_concentration);
        Reactor.SetVal("X_BA:concentration",v_X_BA_concentration);
        Reactor.SetVal("X_ND:concentration",v_X_ND_concentration);

        Reactor.SetVal("Storage",v_t_storage);
        //Reactor.SetVal("Storage",v_r_storage);
        //if (i==0) Reactor.SetVal("constant_inflow",v_r_constant_flow);
        //Reactor.SetVal("constant_inflow",v_r_constant_flow);
        //if (i==4) Reactor.SetVal("S_M:External mass flow time-series",Workingfolder + "S_M_mfr.csv");
        Reactor.SetVal("x",400*(i-n_tanks));
        if (i==0) Reactor.SetVal("y",1200);
        else if (i>0) Reactor.SetVal("y",800);

        system->AddBlock(Reactor,false);

        // Aeration to all tanks
        //system->block("Reactor(" + aquiutils::numbertostring(i+1)+")")->SetProperty("S_O:external_source","Aeration"); // Reactor does not have source "Aeration", so we have to call it!

        // Aeration to selected tanks by aeration_v
        if (aeration_v[i]) system->block("Reactor(" + aquiutils::numbertostring(i+1)+")")->SetProperty("S_O:external_source","Aeration"); // Reactor does not have source "Aeration", so we have to call it!
    }


    // Producing Constant Inflows

    CTimeSeries<double> S_S_inflow_concentration;
    S_S_inflow_concentration.CreateConstant(0,Simulation_time, v_S_S_concentration);
    S_S_inflow_concentration.writefile(Workingfolder + "S_S_constant_inflow_concentration.txt");

    CTimeSeries<double> X_S_inflow_concentration;
    X_S_inflow_concentration.CreateConstant(0,Simulation_time, v_X_S_concentration);
    X_S_inflow_concentration.writefile(Workingfolder + "X_S_constant_inflow_concentration.txt");

    CTimeSeries<double> X_p_inflow_concentration;
    X_p_inflow_concentration.CreateConstant(0,Simulation_time, v_X_p_concentration);
    X_p_inflow_concentration.writefile(Workingfolder + "X_p_constant_inflow_concentration.txt");

    CTimeSeries<double> S_NO_inflow_concentration;
    S_NO_inflow_concentration.CreateConstant(0,Simulation_time, v_S_NO_concentration);
    S_NO_inflow_concentration.writefile(Workingfolder + "S_NO_constant_inflow_concentration.txt");

    CTimeSeries<double> S_NH_inflow_concentration;
    S_NH_inflow_concentration.CreateConstant(0,Simulation_time, v_S_NH_concentration);
    S_NH_inflow_concentration.writefile(Workingfolder + "S_NH_constant_inflow_concentration.txt");

    CTimeSeries<double> S_ND_inflow_concentration;
    S_ND_inflow_concentration.CreateConstant(0,Simulation_time, v_S_ND_concentration);
    S_ND_inflow_concentration.writefile(Workingfolder + "S_ND_constant_inflow_concentration.txt");

    CTimeSeries<double> X_ND_inflow_concentration;
    X_ND_inflow_concentration.CreateConstant(0,Simulation_time, v_X_ND_concentration);
    X_ND_inflow_concentration.writefile("/home/behzad/Projects/ASM_Models/X_ND_constant_inflow_concentration.txt");


    // Constituents Inflow Calculations from Data (Time Variable)


    CTimeSeriesSet<double> Inflow_DeNit(Workingfolder + "Data/DeNit_Influent_Lump.txt",true); // Inflow (Q)
    CTimeSeriesSet<double> Inflow_DeNit_MeOH(Workingfolder + "Data/DeNit_MeOH.txt",true); // Methanol (Externtal S_M)
    CTimeSeriesSet<double> Inflow_DeNit_wasteflow(Workingfolder + "Data/DeNit_wasteflow.txt",true); // Wasteflow (WAS)
    CTimeSeriesSet<double> Inflow_DeNit_returnflow(Workingfolder + "Data/DeNit_returnflow.txt",true); // Returnflow (RAS)


    // Determining Coefficients

    const double c_S_i=p_31; // S_i, 8 //--
    const double c_S_S=1-p_31; // S_S, 8
    const double c_X_S=0.71*p_32*(1-p_2); // X_S, 2
    const double c_X_p=0.71*p_32*p_2; // X_p, 2
    const double c_S_NO=1; // S_NO, 10
    const double c_S_NH=1; // S_NH, 9
    const double c_S_ND=p_36*(1-p_31); // S_ND, 8
    const double c_X_ND=p_3*p_32; // X_ND, 2

    // Inflows

    CTimeSeries<double> Inflow_Q; // Discharge (m3/day)
    CTimeSeries<double> Inflow_WAS; // Wasteflow (WAS) (m3/day)
    CTimeSeries<double> Inflow_RAS; // Returnflow (RAS) (m3/day)

    CTimeSeries<double> Flow_r_r_st; // Reactor to Reactor + Reactor to Settling element top flow (m3/day)
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
    Inflow_WAS=Inflow_DeNit_wasteflow.BTC[0]; // Wasteflow (WAS) (m3/day)
    Inflow_RAS=Inflow_DeNit_returnflow.BTC[0]; // Returnflow (RAS) (m3/day)

    Flow_r_r_st=Inflow_Q+Inflow_RAS;
    Flow_st_sb=Inflow_RAS+Inflow_WAS;
    Flow_st_c=Inflow_Q-Inflow_WAS;

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

    Inflow_Q.writefile(Workingfolder + "Q_time_variable_inflow.txt"); //

    Inflow_WAS.writefile(Workingfolder + "WAS_time_variable_flow.txt"); //
    Inflow_RAS.writefile(Workingfolder + "RAS_time_variable_flow.txt"); //

    Flow_r_r_st.writefile(Workingfolder + "r_r_st_time_variable_flow.txt"); //
    Flow_st_sb.writefile(Workingfolder + "st_sb_time_variable_flow.txt"); //
    Flow_st_c.writefile(Workingfolder + "st_c_time_variable_flow.txt"); //

    Inflow_S_i.writefile(Workingfolder + "S_i_Inflow_concentration.txt");
    Inflow_S_S.writefile(Workingfolder + "S_S_Inflow_concentration.txt");
    Inflow_X_S.writefile(Workingfolder + "X_S_Inflow_concentration.txt");
    Inflow_X_p.writefile(Workingfolder + "X_p_Inflow_concentration.txt");
    Inflow_S_NO.writefile(Workingfolder + "S_NO_Inflow_concentration.txt");
    Inflow_S_NH.writefile(Workingfolder + "S_NH_Inflow_concentration.txt");
    Inflow_S_ND.writefile(Workingfolder + "S_ND_Inflow_concentration.txt");
    Inflow_X_ND.writefile(Workingfolder + "X_ND_Inflow_concentration.txt");

    Inflow_MeOH.writefile(Workingfolder + "S_M_mfr.txt"); // Methanol

    // Assigning Constituents Inflow concentrations to the system (Reactor)

    // Constant inflows by given concentration
/*
    system->block("Reactor")->SetProperty("S_S:constant_inflow_concentration",Workingfolder + "S_S_constant_inflow_concentration.txt");
    system->block("Reactor")->SetProperty("X_S:constant_inflow_concentration",Workingfolder + "X_S_constant_inflow_concentration.txt");
    system->block("Reactor")->SetProperty("X_p:constant_inflow_concentration",Workingfolder + "X_p_constant_inflow_concentration.txt");
    system->block("Reactor")->SetProperty("S_NO:constant_inflow_concentration",Workingfolder + "S_NO_constant_inflow_concentration.txt");
    system->block("Reactor")->SetProperty("S_NH:constant_inflow_concentration",Workingfolder + "S_NH_constant_inflow_concentration.txt");
    system->block("Reactor")->SetProperty("S_ND:constant_inflow_concentration",Workingfolder + "S_ND_constant_inflow_concentration.txt");
    system->block("Reactor")->SetProperty("X_ND:constant_inflow_concentration",Workingfolder + "X_ND_constant_inflow_concentration.txt");
*/

    // Time Variable inflows according DeNite data

    system->block("Reactor(1)")->SetProperty("time_variable_inflow",Workingfolder + "Q_time_variable_inflow.txt"); // Discharge (m3/day)

    system->block("Reactor(1)")->SetProperty("S_S:time_variable_inflow_concentration",Workingfolder + "S_S_Inflow_concentration.txt");
    system->block("Reactor(1)")->SetProperty("X_S:time_variable_inflow_concentration",Workingfolder + "X_S_Inflow_concentration.txt");
    system->block("Reactor(1)")->SetProperty("X_p:time_variable_inflow_concentration",Workingfolder + "X_p_Inflow_concentration.txt");
    system->block("Reactor(1)")->SetProperty("S_NO:time_variable_inflow_concentration",Workingfolder + "S_NO_Inflow_concentration.txt");
    system->block("Reactor(1)")->SetProperty("S_NH:time_variable_inflow_concentration",Workingfolder + "S_NH_Inflow_concentration.txt");
    system->block("Reactor(1)")->SetProperty("S_ND:time_variable_inflow_concentration",Workingfolder + "S_ND_Inflow_concentration.txt");
    system->block("Reactor(1)")->SetProperty("X_ND:time_variable_inflow_concentration",Workingfolder + "X_ND_Inflow_concentration.txt");

    system->block("Reactor(5)")->SetProperty("S_M:external_mass_flow_timeseries",Workingfolder + "S_M_mfr.txt");


    /*
    // Links for Reactor
    for (int i=0; i<n_tanks-1; i++)
    {   Link l_r_st;
        l_r_st.SetQuantities(system, "Fixed flow");
        l_r_st.SetName("Reactor(" + aquiutils::numbertostring(i+1)+"_" + aquiutils::numbertostring(i+2) + ")");
        l_r_st.SetType("Fixed flow");
        l_r_st.SetVal("flow", v_t_t_st_constant_flow);
        system->AddLink(l_r_st, "Reactor(" + aquiutils::numbertostring(i+1)+")", "Reactor(" + aquiutils::numbertostring(i+2)+")", false);
    }

    // Links
    Link l_r_st;
    l_r_st.SetQuantities(system, "Fixed flow");
    l_r_st.SetName("Reactor(" + aquiutils::numbertostring(n_tanks) + ")-Settling element top");
    l_r_st.SetType("Fixed flow");
    l_r_st.SetVal("flow", v_t_t_st_constant_flow);
    system->AddLink(l_r_st, "Reactor(" + aquiutils::numbertostring(n_tanks) + ")", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Fixed flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Fixed flow");
    l_st_c.SetVal("flow", v_st_c_constant_flow);
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface");
    l_st_sb.SetVal("flow", v_st_sb_constant_flow);
    l_st_sb.SetVal("area", v_st_sb_area);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Fixed flow");
    l_sb_r.SetName("Settling element bottom - Reactor(1)");
    l_sb_r.SetType("Fixed flow");
    l_sb_r.SetVal("flow", v_sb_r_constant_flow);
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor(1)", false);

    Link l_sb_was;
    l_sb_was.SetQuantities(system, "Fixed flow");
    l_sb_was.SetName("Settling element bottom - WAS");
    l_sb_was.SetType("Fixed flow");
    l_sb_was.SetVal("flow", v_sb_was_constant_flow);
    system->AddLink(l_sb_was, "Settling element bottom", "WAS", false);
    */

    // Time-Dependent flow Links for Reactor
    for (int i=0; i<n_tanks-1; i++)
    {   Link l_r_st;
        l_r_st.SetQuantities(system, "Time-Dependent flow");
        l_r_st.SetName("Reactor(" + aquiutils::numbertostring(i+1)+"_" + aquiutils::numbertostring(i+2) + ")");
        l_r_st.SetType("Time-Dependent flow");
        l_r_st.SetProperty("flow", Workingfolder + "r_r_st_time_variable_flow.txt");
        system->AddLink(l_r_st, "Reactor(" + aquiutils::numbertostring(i+1)+")", "Reactor(" + aquiutils::numbertostring(i+2)+")", false);
    }

    // Time-Dependent flow Links
    Link l_r_st;
    l_r_st.SetQuantities(system, "Time-Dependent flow");
    l_r_st.SetName("Reactor(" + aquiutils::numbertostring(n_tanks) + ") - Settling element top");
    l_r_st.SetType("Time-Dependent flow");
    l_r_st.SetProperty("flow", Workingfolder + "r_r_st_time_variable_flow.txt");
    system->AddLink(l_r_st, "Reactor(" + aquiutils::numbertostring(n_tanks) + ")", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Time-Dependent flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Time-Dependent flow");
    l_st_c.SetProperty("flow", Workingfolder + "st_c_time_variable_flow.txt");
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_sb_was;
    l_sb_was.SetQuantities(system, "Time-Dependent flow");
    l_sb_was.SetName("Settling element bottom - WAS");
    l_sb_was.SetType("Time-Dependent flow");
    l_sb_was.SetProperty("flow", "/home/behzad/Projects/ASM_Models/WAS_time_variable_flow.txt");
    system->AddLink(l_sb_was, "Settling element bottom", "WAS", false);

    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Time-Dependent flow");
    l_sb_r.SetName("Settling element bottom - Reactor(1)");
    l_sb_r.SetType("Time-Dependent flow");
    l_sb_r.SetProperty("flow", Workingfolder + "RAS_time_variable_flow.txt");
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor(1)", false);

    // Time-Dependent interface Link
    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface (time variable)");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface (time variable)");
    l_st_sb.SetProperty("flow", "/home/behzad/Projects/ASM_Models/st_sb_time_variable_flow.txt");
    l_st_sb.SetVal("area", v_st_sb_area);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    /*
    // Observations
    Observation total_inflow;

    total_inflow.SetQuantities(system, "Observation");
    total_inflow.SetProperty("expression","inflow");
    total_inflow.SetProperty("object","Reactor (1)");
    total_inflow.SetName("Inflow");
    total_inflow.SetType("Observation");
    system->AddObservation(total_inflow,false);


    Observation s_inflow_concentration;

    s_inflow_concentration.SetQuantities(system, "Observation");
    s_inflow_concentration.SetProperty("expression","Solids:inflow_concentration");
    s_inflow_concentration.SetProperty("object","Reactor (1)");
    s_inflow_concentration.SetName("SolidsInflowConcentration");
    s_inflow_concentration.SetType("Observation");
    system->AddObservation(s_inflow_concentration,false);

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









// ----- Producing OUProcessed Random Inflow -----

/*

    CTimeSeries<double> SolidConcentrationNS;
    SolidConcentrationNS.CreateOUProcess(0,Simulation_time,0.05,1);
    SolidConcentrationNS.writefile(Workingfolder + "inflow_concentration_NS.csv");
    vector<double> params; params.push_back(3); params.push_back(1);
    CTimeSeries<double> SolidConcentration = SolidConcentrationNS.MapfromNormalScoreToDistribution("lognormal", params);
//Reactor.Variable("Solids:inflow_concentration")->SetTimeSeries(SolidConcentration);
    SolidConcentration.writefile(Workingfolder + "inflow_concentration.csv");
    Reactor.SetProperty("Solids:inflow_concentration",Workingfolder + "inflow_concentration.csv");

/*
    CTimeSeries<double> InflowNS;
    InflowNS.CreateOUProcess(0,100,0.05,1);
    vector<double> i_params; i_params.push_back(1.5); i_params.push_back(0.7);
    CTimeSeries<double> Inflow = InflowNS.MapfromNormalScoreToDistribution("lognormal", i_params);
    //Reactor.Variable("inflow")->SetTimeSeries(Inflow);
    Inflow.writefile(Workingfolder + "inflow.csv");
    Reactor.SetProperty("inflow",Workingfolder + "inflow.csv");
    */
