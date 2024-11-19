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

    //Model properties
    const double S_s_concentration=0;
    const double X_b_concentration=1;
    const double r_storage=1000;


    //Model configuration
    Constituent S_s;
    S_s.SetQuantities(system, "Constituent");
    S_s.SetName("S_s");
    S_s.SetType("Constituent");
    system->AddConstituent(S_s,false);

    Constituent X_b;
    X_b.SetQuantities(system, "Constituent");
    X_b.SetName("X_b");
    X_b.SetType("Constituent");
    system->AddConstituent(X_b,false);

    /*
    Constituent Settled_Particles;
    Settled_Particles.SetQuantities(system, "Immobile_Constituent");
    Settled_Particles.SetName("Settled_Particles");
    Settled_Particles.SetType("Immobile_Constituent");
    system->AddConstituent(Settled_Particles,false);
     */

    //RxnParameter Rxparam_K;
    //Rxparam_K.SetQuantities(system,"ReactionParameter");
    //Rxparam_K.SetName("K");
    //Rxparam_K.SetVal("base_value",0.1);
    //system->AddReactionParameter(Rxparam_K, false);

    RxnParameter mu;
    mu.SetQuantities(system,"ReactionParameter");
    mu.SetName("mu");
    mu.SetVal("base_value",2);
    system->AddReactionParameter(mu, false);

    RxnParameter Ks;
    Ks.SetQuantities(system,"ReactionParameter");
    Ks.SetName("Ks");
    Ks.SetVal("base_value",20);
    system->AddReactionParameter(Ks, false);

    RxnParameter Y;
    Y.SetQuantities(system,"ReactionParameter");
    Y.SetName("Y");
    Y.SetVal("base_value",0.5);
    system->AddReactionParameter(Y, false);

    RxnParameter b;
    b.SetQuantities(system,"ReactionParameter");
    b.SetName("b");
    b.SetVal("base_value",0.3);
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
    Stl_element_top.SetVal("Storage",20);
    Stl_element_top.SetVal("bottom_elevation",0);
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
    Stl_element_bottom.SetVal("Storage",20);
    Stl_element_bottom.SetVal("bottom_elevation",0);
    Stl_element_bottom.SetVal("x",800);
    Stl_element_bottom.SetVal("y",1000);
    system->AddBlock(Stl_element_bottom,false);

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

    Block Reactor;
    Reactor.SetQuantities(system, "Reactor");
    Reactor.SetName("Reactor");
    Reactor.SetType("Reactor");
    Reactor.SetProperty("S_s:Inflow concentration","/home/behzad/Projects/ASM_Models/ASM GUI Model/Inflow_S.txt");
    Reactor.SetVal("S_s:concentration",S_s_concentration);
    Reactor.SetVal("X_b:concentration",X_b_concentration);
    Reactor.SetVal("Storage",r_storage);
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

    Reactor.SetVal("constant_inflow",r_inflow);

    system->AddBlock(Reactor,false);

    //system->block("Reactor (1)")->SetProperty("inflow","/home/behzad/Projects/ASM_Models/inflow.txt");

    Link l_r_st;
    l_r_st.SetQuantities(system, "Fixed flow");
    l_r_st.SetName("Reactor - Settling element top");
    l_r_st.SetType("Fixed flow");
    l_r_st.SetVal("flow", r_inflow);
    system->AddLink(l_r_st, "Reactor", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Fixed flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Fixed flow");
    l_st_c.SetVal("flow", r_inflow);
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface");
    l_st_sb.SetVal("flow", r_inflow);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Fixed flow");
    l_sb_r.SetName("Settling element bottom - Reactor");
    l_sb_r.SetType("Fixed flow");
    l_sb_r.SetVal("flow", r_inflow);
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor", false);

    Link l_st_was;
    l_st_was.SetQuantities(system, "Fixed flow");
    l_st_was.SetName("Settling element bottom - WAS");
    l_st_was.SetType("Fixed flow");
    l_st_was.SetVal("flow", r_inflow);
    system->AddLink(l_st_was, "Settling element bottom", "WAS", false);

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
