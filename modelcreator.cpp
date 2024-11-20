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
    system->GetQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/main_components.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/unsaturated_soil.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Well.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/wastewater.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/mass_transfer.json");
    system->ReadSystemSettingsTemplate("/home/behzad/Projects/OpenHydroQual/resources/settings.json");

    const double Simulation_time=100; // Simulation Time in Days

    //Model Properties
    const double v_settling_vel=10000; // X_b : Unit: m/day
    const double v_S_S_concentration=0;
    const double v_X_b_concentration=1;
    const double v_r_storage=1000;
    const double v_mu=2;
    const double v_Ks=20;
    const double v_Y=0.5;
    const double v_b=0.3;
    const double v_r_constant_flow=800; // Reactor: Constant flow : Unit: m3/day
    const double v_s_t_storage=200; // Settling element top: initial storage
    const double v_s_t_bottom_elevation=1; // Settling element top: bottom elevation
    const double v_s_b_storage=200; // Settling element bottom: initial storage
    const double v_s_b_bottom_elevation=0; // Settling element bottom: bottom elevation
    const double v_r_st_flow=1700; // Link: Reactor to Settling element top: flow
    const double v_st_c_flow=750; // Link: Settling element top to Clarifier: flow
    const double v_st_sb_flow=950; // Link: Settling element top to Settling element bottom: flow
    const double v_st_sb_area=100; // Link: Settling element top to Settling element bottom: area
    const double v_sb_r_flow=900; // Link: Settling element bottom to Reactor: flow
    const double v_sb_was_flow=50; // Link: Settling element bottom to WAS: flow



    //Model Configuration

    //Consistuents
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

    Constituent X_BH;
    X_BH.SetQuantities(system, "Particle");
    X_BH.SetName("X_BH");
    X_BH.SetType("Particle");
    X_BH.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_BH,false);

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

    Constituent X_BM;
    X_BM.SetQuantities(system, "Particle");
    X_BM.SetName("X_BM");
    X_BM.SetType("Particle");
    X_BM.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_BM,false);

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

    Constituent S_ND;
    S_ND.SetQuantities(system, "Constituent");
    S_ND.SetName("S_ND");
    S_ND.SetType("Constituent");
    system->AddConstituent(S_ND,false);

    //Reaction Parameters
    RxnParameter mu;
    mu.SetQuantities(system,"ReactionParameter");
    mu.SetName("mu");
    mu.SetVal("base_value",v_mu);
    system->AddReactionParameter(mu, false);

    RxnParameter Ks;
    Ks.SetQuantities(system,"ReactionParameter");
    Ks.SetName("Ks");
    Ks.SetVal("base_value",v_Ks);
    system->AddReactionParameter(Ks, false);

    RxnParameter Y;
    Y.SetQuantities(system,"ReactionParameter");
    Y.SetName("Y");
    Y.SetVal("base_value",v_Y);
    system->AddReactionParameter(Y, false);

    RxnParameter b;
    b.SetQuantities(system,"ReactionParameter");
    b.SetName("b");
    b.SetVal("base_value",v_b);
    system->AddReactionParameter(b, false);

    Reaction Growth;
    Growth.SetQuantities(system,"Reaction");
    Growth.SetName("Growth");
    Growth.SetProperty("S_S:stoichiometric_constant","(0-1/Y)");
    Growth.SetProperty("X_b:stoichiometric_constant","1");
    Growth.SetProperty("rate_expression","(mu*S_S/(Ks+S_S)*X_b)");
    system->AddReaction(Growth,false);

    Reaction Decay;
    Decay.SetQuantities(system,"Reaction");
    Decay.SetName("Decay");
    Decay.SetProperty("S_S:stoichiometric_constant","1");
    Decay.SetProperty("X_b:stoichiometric_constant","-1");
    Decay.SetProperty("rate_expression","(b*X_b)");
    system->AddReaction(Decay,false);

    //Settling Elements
    Block Stl_element_top;
    Stl_element_top.SetQuantities(system, "Settling element");
    Stl_element_top.SetName("Settling element top");
    Stl_element_top.SetType("Settling element");
    Stl_element_top.SetVal("Coagulant:concentration",0);

    CTimeSeries<double> CoagNS;
    CoagNS.CreateOUProcess(0,Simulation_time,0.05,1);
    CoagNS.writefile("/home/behzad/Projects/ASM_Models/coagulant_mfr_NS.csv");
    vector<double> c_params; c_params.push_back(2.5); c_params.push_back(0.6);
    CTimeSeries<double> Coag = CoagNS.MapfromNormalScoreToDistribution("lognormal", c_params);
    //Reactor.Variable("Coagulant:external_mass_flow_timeseries")->SetTimeSeries(Coag);
    Coag.writefile("/home/behzad/Projects/ASM_Models/coagulant_mfr.csv");
    Stl_element_top.SetProperty("Coagulant:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/coagulant_mfr.csv");

    //Stl_element_top.SetProperty("Coagulant:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/coagulant_mfr.txt");

    Stl_element_top.SetVal("Settled_Particles:concentration",0);
    Stl_element_top.SetVal("Solids:concentration",0);
    Stl_element_top.SetVal("Storage",v_s_t_storage);
    Stl_element_top.SetVal("bottom_elevation",v_s_t_bottom_elevation);
    Stl_element_top.SetVal("x",800);
    Stl_element_top.SetVal("y",600);
    system->AddBlock(Stl_element_top,false);

    Block Stl_element_bottom;
    Stl_element_bottom.SetQuantities(system, "Settling element");
    Stl_element_bottom.SetName("Settling element bottom");
    Stl_element_bottom.SetType("Settling element");
    Stl_element_bottom.SetVal("Coagulant:concentration",0);

    CTimeSeries<double> CoagNS2;
    CoagNS2.CreateOUProcess(0,Simulation_time,0.05,1);
    CoagNS2.writefile("/home/behzad/Projects/ASM_Models/coagulant_mfr_NS.csv");
    vector<double> c_params2; c_params2.push_back(2.5); c_params2.push_back(0.6);
    CTimeSeries<double> Coag2 = CoagNS2.MapfromNormalScoreToDistribution("lognormal", c_params);
    //Reactor.Variable("Coagulant:external_mass_flow_timeseries")->SetTimeSeries(Coag);
    Coag2.writefile("/home/behzad/Projects/ASM_Models/coagulant_mfr.csv");
    Stl_element_bottom.SetProperty("Coagulant:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/coagulant_mfr.csv");

    //Stl_element_bottom.SetProperty("Coagulant:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/coagulant_mfr.txt");

    Stl_element_bottom.SetVal("Settled_Particles:concentration",0);
    Stl_element_bottom.SetVal("Solids:concentration",0);
    Stl_element_bottom.SetVal("Storage",v_s_b_storage);
    Stl_element_bottom.SetVal("bottom_elevation",v_s_b_bottom_elevation);
    Stl_element_bottom.SetVal("x",800);
    Stl_element_bottom.SetVal("y",1000);
    system->AddBlock(Stl_element_bottom,false);

    //Fixed Head Blocks
    Block fh_clarifier;
    fh_clarifier.SetQuantities(system, "fixed_head");
    fh_clarifier.SetName("Clarifier");
    fh_clarifier.SetType("fixed_head");
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
    fh_WAS.SetVal("S_S:concentration",0);
    fh_WAS.SetVal("X_b:concentration",0);
    fh_WAS.SetVal("Storage",100000);
    fh_WAS.SetVal("x",1200);
    fh_WAS.SetVal("y",1000);
    system->AddBlock(fh_WAS,false);

    //Reactor Block
    Block Reactor;
    Reactor.SetQuantities(system, "Reactor");
    Reactor.SetName("Reactor");
    Reactor.SetType("Reactor");
    Reactor.SetProperty("S_S:inflow_concentration","/home/behzad/Projects/ASM_Models/ASM GUI Model/Inflow_S.txt");
    Reactor.SetVal("S_S:concentration",v_S_S_concentration);
    Reactor.SetVal("X_b:concentration",v_X_b_concentration);
    Reactor.SetVal("Storage",v_r_storage);
    Reactor.SetVal("constant_inflow",v_r_constant_flow);
    Reactor.SetVal("x",400);
    Reactor.SetVal("y",800);

    CTimeSeries<double> SolidConcentrationNS;
    SolidConcentrationNS.CreateOUProcess(0,Simulation_time,0.05,1);
    SolidConcentrationNS.writefile("/home/behzad/Projects/ASM_Models/inflow_concentration_NS.csv");
    vector<double> params; params.push_back(3); params.push_back(1);
    CTimeSeries<double> SolidConcentration = SolidConcentrationNS.MapfromNormalScoreToDistribution("lognormal", params);
    //Reactor.Variable("Solids:inflow_concentration")->SetTimeSeries(SolidConcentration);
    SolidConcentration.writefile("/home/behzad/Projects/ASM_Models/inflow_concentration.csv");
    Reactor.SetProperty("Solids:inflow_concentration","/home/behzad/Projects//inflow_concentration.csv");

    //Reactor.SetProperty("Solids:inflow_concentration","/home/behzad/Projects/ASM_Models/inflow_concentration.txt");

    /*CTimeSeries<double> InflowNS;
    InflowNS.CreateOUProcess(0,100,0.05,1);
    vector<double> i_params; i_params.push_back(1.5); i_params.push_back(0.7);
    CTimeSeries<double> Inflow = InflowNS.MapfromNormalScoreToDistribution("lognormal", i_params);
    //Reactor.Variable("inflow")->SetTimeSeries(Inflow);
    Inflow.writefile("/home/behzad/Projects/ASM_Models/inflow.csv");
    Reactor.SetProperty("inflow","/home/behzad/Projects/ASM_Models/inflow.csv");
    */

    // ----- Producing Constant Random Inflow -----

    // Seed the random number generator with the current time
    srand(time(0));

    /*// Generate a random number between 1 and 10
    int r_inflow = rand() % 10 + 1;*/

    // Generate a random double between 0 and 1
    double randomValue = (double)rand() / RAND_MAX;

    // Scale the random value to the range 1-10
    double r_inflow = randomValue * 9 + 1;

    CTimeSeries<double> inflow_timeseries;
    inflow_timeseries.CreateConstant(0,Simulation_time, r_inflow);
    inflow_timeseries.writefile("/home/behzad/Projects/ASM_Models/inflow.csv");

    //Reactor.SetVal("constant_inflow",r_inflow);

    system->AddBlock(Reactor,false);

    //system->block("Reactor (1)")->SetProperty("inflow","/home/behzad/Projects/ASM_Models/inflow.txt");

    //Links
    Link l_r_st;
    l_r_st.SetQuantities(system, "Fixed flow");
    l_r_st.SetName("Reactor - Settling element top");
    l_r_st.SetType("Fixed flow");
    l_r_st.SetVal("flow", v_r_st_flow);
    system->AddLink(l_r_st, "Reactor", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Fixed flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Fixed flow");
    l_st_c.SetVal("flow", v_st_c_flow);
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface");
    l_st_sb.SetVal("flow", v_st_sb_flow);
    l_st_sb.SetVal("area", v_st_sb_area);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Fixed flow");
    l_sb_r.SetName("Settling element bottom - Reactor");
    l_sb_r.SetType("Fixed flow");
    l_sb_r.SetVal("flow", v_sb_r_flow);
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor", false);

    Link l_sb_was;
    l_sb_was.SetQuantities(system, "Fixed flow");
    l_sb_was.SetName("Settling element bottom - WAS");
    l_sb_was.SetType("Fixed flow");
    l_sb_was.SetVal("flow", v_sb_was_flow);
    system->AddLink(l_sb_was, "Settling element bottom", "WAS", false);

    //Observations
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

    system->SetSettingsParameter("simulation_end_time",Simulation_time);
    system->SetSystemSettings();
    cout<<"Populate functions"<<endl;
    system->PopulateOperatorsFunctions();
    cout<<"Variable parents"<<endl;
    system->SetVariableParents();
    return true;
}

















/*
// Simple Model

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
    system->GetQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/main_components.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/unsaturated_soil.json");
    //system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/Well.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/wastewater.json");
    system->AppendQuanTemplate("/home/behzad/Projects/OpenHydroQual/resources/mass_transfer.json");
    system->ReadSystemSettingsTemplate("/home/behzad/Projects/OpenHydroQual/resources/settings.json");

    const double Simulation_time=100; // Simulation Time in Days

    //Model Properties
    const double v_settling_vel=10000; // Unit: m/day
    const double v_S_s_concentration=0;
    const double v_X_b_concentration=1;
    const double v_r_storage=1000;
    const double v_mu=2;
    const double v_Ks=20;
    const double v_Y=0.5;
    const double v_b=0.3;
    const double v_r_constant_flow=800; //Unit: m3/day
    const double v_s_t_storage=200; //Settling element top: initial storage
    const double v_s_t_bottom_elevation=1; //Settling element top: bottom elevation
    const double v_s_b_storage=200; //Settling element bottom: initial storage
    const double v_s_b_bottom_elevation=0; //Settling element bottom: bottom elevation
    const double v_r_st_flow=1700; //Link: Reactor to Settling element top: flow
    const double v_st_c_flow=750; //Link: Settling element top to Clarifier: flow
    const double v_st_sb_flow=950; //Link: Settling element top to Settling element bottom: flow
    const double v_st_sb_area=100; //Link: Settling element top to Settling element bottom: area
    const double v_sb_r_flow=900; //Link: Settling element bottom to Reactor: flow
    const double v_sb_was_flow=50; //Link: Settling element bottom to WAS: flow



    //Model Configuration

    //Consistuents
    Constituent S_s;
    S_s.SetQuantities(system, "Constituent");
    S_s.SetName("S_s");
    S_s.SetType("Constituent");
    system->AddConstituent(S_s,false);

    Constituent X_b;
    X_b.SetQuantities(system, "Particle");
    X_b.SetName("X_b");
    X_b.SetType("Particle");
    X_b.SetVal("settling_velocity",v_settling_vel);
    system->AddConstituent(X_b,false);

    //Reaction Parameters
    RxnParameter mu;
    mu.SetQuantities(system,"ReactionParameter");
    mu.SetName("mu");
    mu.SetVal("base_value",v_mu);
    system->AddReactionParameter(mu, false);

    RxnParameter Ks;
    Ks.SetQuantities(system,"ReactionParameter");
    Ks.SetName("Ks");
    Ks.SetVal("base_value",v_Ks);
    system->AddReactionParameter(Ks, false);

    RxnParameter Y;
    Y.SetQuantities(system,"ReactionParameter");
    Y.SetName("Y");
    Y.SetVal("base_value",v_Y);
    system->AddReactionParameter(Y, false);

    RxnParameter b;
    b.SetQuantities(system,"ReactionParameter");
    b.SetName("b");
    b.SetVal("base_value",v_b);
    system->AddReactionParameter(b, false);

    Reaction Growth;
    Growth.SetQuantities(system,"Reaction");
    Growth.SetName("Growth");
    Growth.SetProperty("S_s:stoichiometric_constant","(0-1/Y)");
    Growth.SetProperty("X_b:stoichiometric_constant","1");
    Growth.SetProperty("rate_expression","(mu*S_s/(Ks+S_s)*X_b)");
    system->AddReaction(Growth,false);

    Reaction Decay;
    Decay.SetQuantities(system,"Reaction");
    Decay.SetName("Decay");
    Decay.SetProperty("S_s:stoichiometric_constant","1");
    Decay.SetProperty("X_b:stoichiometric_constant","-1");
    Decay.SetProperty("rate_expression","(b*X_b)");
    system->AddReaction(Decay,false);

    //Settling Elements
    Block Stl_element_top;
    Stl_element_top.SetQuantities(system, "Settling element");
    Stl_element_top.SetName("Settling element top");
    Stl_element_top.SetType("Settling element");
    Stl_element_top.SetVal("Coagulant:concentration",0);

    CTimeSeries<double> CoagNS;
    CoagNS.CreateOUProcess(0,Simulation_time,0.05,1);
    CoagNS.writefile("/home/behzad/Projects/ASM_Models/coagulant_mfr_NS.csv");
    vector<double> c_params; c_params.push_back(2.5); c_params.push_back(0.6);
    CTimeSeries<double> Coag = CoagNS.MapfromNormalScoreToDistribution("lognormal", c_params);
    //Reactor.Variable("Coagulant:external_mass_flow_timeseries")->SetTimeSeries(Coag);
    Coag.writefile("/home/behzad/Projects/ASM_Models/coagulant_mfr.csv");
    Stl_element_top.SetProperty("Coagulant:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/coagulant_mfr.csv");

    //Stl_element_top.SetProperty("Coagulant:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/coagulant_mfr.txt");

    Stl_element_top.SetVal("Settled_Particles:concentration",0);
    Stl_element_top.SetVal("Solids:concentration",0);
    Stl_element_top.SetVal("Storage",v_s_t_storage);
    Stl_element_top.SetVal("bottom_elevation",v_s_t_bottom_elevation);
    Stl_element_top.SetVal("x",800);
    Stl_element_top.SetVal("y",600);
    system->AddBlock(Stl_element_top,false);

    Block Stl_element_bottom;
    Stl_element_bottom.SetQuantities(system, "Settling element");
    Stl_element_bottom.SetName("Settling element bottom");
    Stl_element_bottom.SetType("Settling element");
    Stl_element_bottom.SetVal("Coagulant:concentration",0);

    CTimeSeries<double> CoagNS2;
    CoagNS2.CreateOUProcess(0,Simulation_time,0.05,1);
    CoagNS2.writefile("/home/behzad/Projects/ASM_Models/coagulant_mfr_NS.csv");
    vector<double> c_params2; c_params2.push_back(2.5); c_params2.push_back(0.6);
    CTimeSeries<double> Coag2 = CoagNS2.MapfromNormalScoreToDistribution("lognormal", c_params);
    //Reactor.Variable("Coagulant:external_mass_flow_timeseries")->SetTimeSeries(Coag);
    Coag2.writefile("/home/behzad/Projects/ASM_Models/coagulant_mfr.csv");
    Stl_element_bottom.SetProperty("Coagulant:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/coagulant_mfr.csv");

    //Stl_element_bottom.SetProperty("Coagulant:external_mass_flow_timeseries","/home/behzad/Projects/ASM_Models/coagulant_mfr.txt");

    Stl_element_bottom.SetVal("Settled_Particles:concentration",0);
    Stl_element_bottom.SetVal("Solids:concentration",0);
    Stl_element_bottom.SetVal("Storage",v_s_b_storage);
    Stl_element_bottom.SetVal("bottom_elevation",v_s_b_bottom_elevation);
    Stl_element_bottom.SetVal("x",800);
    Stl_element_bottom.SetVal("y",1000);
    system->AddBlock(Stl_element_bottom,false);

    //Fixed Head Blocks
    Block fh_clarifier;
    fh_clarifier.SetQuantities(system, "fixed_head");
    fh_clarifier.SetName("Clarifier");
    fh_clarifier.SetType("fixed_head");
    fh_clarifier.SetVal("S_s:concentration",0);
    fh_clarifier.SetVal("X_b:concentration",0);
    fh_clarifier.SetVal("Storage",100000);
    fh_clarifier.SetVal("x",1200);
    fh_clarifier.SetVal("y",600);
    system->AddBlock(fh_clarifier,false);

    Block fh_WAS;
    fh_WAS.SetQuantities(system, "fixed_head");
    fh_WAS.SetName("WAS");
    fh_WAS.SetType("fixed_head");
    fh_WAS.SetVal("S_s:concentration",0);
    fh_WAS.SetVal("X_b:concentration",0);
    fh_WAS.SetVal("Storage",100000);
    fh_WAS.SetVal("x",1200);
    fh_WAS.SetVal("y",1000);
    system->AddBlock(fh_WAS,false);

    //Reactor Block
    Block Reactor;
    Reactor.SetQuantities(system, "Reactor");
    Reactor.SetName("Reactor");
    Reactor.SetType("Reactor");
    Reactor.SetProperty("S_s:inflow_concentration","/home/behzad/Projects/ASM_Models/ASM GUI Model/Inflow_S.txt");
    Reactor.SetVal("S_s:concentration",v_S_s_concentration);
    Reactor.SetVal("X_b:concentration",v_X_b_concentration);
    Reactor.SetVal("Storage",v_r_storage);
    Reactor.SetVal("constant_inflow",v_r_constant_flow);
    Reactor.SetVal("x",400);
    Reactor.SetVal("y",800);

    CTimeSeries<double> SolidConcentrationNS;
    SolidConcentrationNS.CreateOUProcess(0,Simulation_time,0.05,1);
    SolidConcentrationNS.writefile("/home/behzad/Projects/ASM_Models/inflow_concentration_NS.csv");
    vector<double> params; params.push_back(3); params.push_back(1);
    CTimeSeries<double> SolidConcentration = SolidConcentrationNS.MapfromNormalScoreToDistribution("lognormal", params);
    //Reactor.Variable("Solids:inflow_concentration")->SetTimeSeries(SolidConcentration);
    SolidConcentration.writefile("/home/behzad/Projects/ASM_Models/inflow_concentration.csv");
    Reactor.SetProperty("Solids:inflow_concentration","/home/behzad/Projects//inflow_concentration.csv");

    //Reactor.SetProperty("Solids:inflow_concentration","/home/behzad/Projects/ASM_Models/inflow_concentration.txt");

    /*CTimeSeries<double> InflowNS;
    InflowNS.CreateOUProcess(0,100,0.05,1);
    vector<double> i_params; i_params.push_back(1.5); i_params.push_back(0.7);
    CTimeSeries<double> Inflow = InflowNS.MapfromNormalScoreToDistribution("lognormal", i_params);
    //Reactor.Variable("inflow")->SetTimeSeries(Inflow);
    Inflow.writefile("/home/behzad/Projects/ASM_Models/inflow.csv");
    Reactor.SetProperty("inflow","/home/behzad/Projects/ASM_Models/inflow.csv");
    */

/*
    // ----- Producing Constant Random Inflow -----

    // Seed the random number generator with the current time
    srand(time(0));
/*

    // Generate a random number between 1 and 10
    int r_inflow = rand() % 10 + 1;*/
/*
    // Generate a random double between 0 and 1
    double randomValue = (double)rand() / RAND_MAX;

    // Scale the random value to the range 1-10
    double r_inflow = randomValue * 9 + 1;
/*
    CTimeSeries<double> inflow_timeseries;
    inflow_timeseries.CreateConstant(0,Simulation_time, r_inflow);
    inflow_timeseries.writefile("/home/behzad/Projects/ASM_Models/inflow.csv");

    //Reactor.SetVal("constant_inflow",r_inflow);

    system->AddBlock(Reactor,false);

    //system->block("Reactor (1)")->SetProperty("inflow","/home/behzad/Projects/ASM_Models/inflow.txt");

    //Links
    Link l_r_st;
    l_r_st.SetQuantities(system, "Fixed flow");
    l_r_st.SetName("Reactor - Settling element top");
    l_r_st.SetType("Fixed flow");
    l_r_st.SetVal("flow", v_r_st_flow);
    system->AddLink(l_r_st, "Reactor", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Fixed flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Fixed flow");
    l_st_c.SetVal("flow", v_st_c_flow);
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface");
    l_st_sb.SetVal("flow", v_st_sb_flow);
    l_st_sb.SetVal("area", v_st_sb_area);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Fixed flow");
    l_sb_r.SetName("Settling element bottom - Reactor");
    l_sb_r.SetType("Fixed flow");
    l_sb_r.SetVal("flow", v_sb_r_flow);
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor", false);

    Link l_sb_was;
    l_sb_was.SetQuantities(system, "Fixed flow");
    l_sb_was.SetName("Settling element bottom - WAS");
    l_sb_was.SetType("Fixed flow");
    l_sb_was.SetVal("flow", v_sb_was_flow);
    system->AddLink(l_sb_was, "Settling element bottom", "WAS", false);

    //Observations
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

    system->SetSettingsParameter("simulation_end_time",Simulation_time);
    system->SetSystemSettings();
    cout<<"Populate functions"<<endl;
    system->PopulateOperatorsFunctions();
    cout<<"Variable parents"<<endl;
    system->SetVariableParents();
    return true;
}

*/
