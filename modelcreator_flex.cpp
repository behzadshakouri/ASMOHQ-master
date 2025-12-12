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

#ifdef PowerEdge
    string Workingfolder = "/mnt/3rd900/Projects/ASM_Models/Flex/";
    string OHQ_r_path = "/mnt/3rd900/Projects/OpenHydroQual/resources/";
    CTimeSeriesSet<double> Inflow_DeNit("/mnt/3rd900/Projects/ASM_Models/Flex/Data/DeNit_Influent_Lump.txt",true); // Inflow (Q)
    CTimeSeriesSet<double> Inflow_DeNit_MeOH("/mnt/3rd900/Projects/ASM_Models/Flex/Data/DeNit_MeOH.txt",true); // Methanol (Externtal S_M)
    CTimeSeriesSet<double> Inflow_DeNit_wasteflow("/mnt/3rd900/Projects/ASM_Models/Flex/Data/DeNit_wasteflow.txt",true); // Wasteflow (WAS)
    CTimeSeriesSet<double> Inflow_DeNit_returnflow("/mnt/3rd900/Projects/ASM_Models/Flex/Data/DeNit_returnflow.txt",true); // Returnflow (RAS)
    CTimeSeriesSet<double> DeNit_Temp("/mnt/3rd900/Projects/ASM_Models/Flex/Data/DeNit_Temp.txt",true); // Temperature

#elif Behzad
    string Workingfolder = "/home/behzad/Projects/ASM_Models/Flex/";
    string OHQ_r_path = "/home/behzad/Projects/OpenHydroQual/resources/";
    CTimeSeriesSet<double> Inflow_DeNit("/home/behzad/Projects/ASM_Models/Data/DeNit_Influent_Lump.txt",true); // Inflow (Q)
    CTimeSeriesSet<double> Inflow_DeNit_MeOH("/home/behzad/Projects/ASM_Models/Data/DeNit_MeOH.txt",true); // Methanol (Externtal S_M)
    CTimeSeriesSet<double> Inflow_DeNit_wasteflow("/home/behzad/Projects/ASM_Models/Data/DeNit_wasteflow.txt",true); // Wasteflow (WAS)
    CTimeSeriesSet<double> Inflow_DeNit_returnflow("/home/behzad/Projects/ASM_Models/Data/DeNit_returnflow.txt",true); // Returnflow (RAS)
    CTimeSeriesSet<double> DeNit_Temp("/home/behzad/Projects/ASM_Models/Data/DeNit_Temp.txt",true); // Temperature

#elif Arash
    string Workingfolder = "/home/arash/Projects/ASM_Models/";
    string OHQ_r_path = "/home/arash/Projects/OpenHydroQual/resources/";
    CTimeSeriesSet<double> Inflow_DeNit("/home/arash/Projects/ASM_Models/Data/DeNit_Influent_Lump.txt",true); // Inflow (Q)
    CTimeSeriesSet<double> Inflow_DeNit_MeOH("/home/arash/Projects/ASM_Models/Data/DeNit_MeOH.txt",true); // Methanol (Externtal S_M)
    CTimeSeriesSet<double> Inflow_DeNit_wasteflow("/home/arash/Projects/ASM_Models/Data/DeNit_wasteflow.txt",true); // Wasteflow (WAS)
    CTimeSeriesSet<double> Inflow_DeNit_returnflow("/home/arash/Projects/ASM_Models/Data/DeNit_returnflow.txt",true); // Returnflow (RAS)
    CTimeSeriesSet<double> DeNit_Temp("/home/arash/Projects/ASM_Models/Data/DeNit_Temp.txt",true); // Temperature
#endif

    system->GetQuanTemplate(OHQ_r_path + "main_components.json");
    //system->AppendQuanTemplate("OHQ_r_path + "unsaturated_soil.json");
    //system->AppendQuanTemplate("OHQ_r_path + "Well.json");
    system->AppendQuanTemplate(OHQ_r_path + "wastewater.json");
    system->AppendQuanTemplate(OHQ_r_path + "mass_transfer.json");
    system->ReadSystemSettingsTemplate(OHQ_r_path + "settings.json");

    bool OUP = false; // True for using OUProcess, False for using DeNite Data

    bool St = true; // True for using Simulation Time is Days, False for using Start and End Date

    bool Calibration = true; // True for Calibration

    double Simulation_start_time=40210; // Simulation Start Date
    double Simulation_end_time=40359; // Simulation End Date

    if (OUP)
    {
        Simulation_start_time=40000; // Simulation Start Date
        Simulation_end_time=40700; // Simulation End Date
    }

    const double Simulation_time_Calc = Simulation_end_time - Simulation_start_time;

    const double Simulation_time=Simulation_time_Calc; // Simulation Time in Days

    if (OUP==false) // Set using Simulation Time is Days
        St=false;

    const double dt = 0.1; // OUP Time Step

    const double Initial_time_step = dt; // Observation writing time step

    CTimeSeriesSet<double> ToWrite;

    //-------------------------------------Temp----------------------------------

    if (OUP)
    {
        CTimeSeries<double> OUP_Temp;

        //Temp
        CTimeSeries<double> Fitted_Sin = DeNit_Temp.BTC[0].CreateSinusoidal(sin_T0,sin_a,sin_b);
        CTimeSeries<double> Residual = Fitted_Sin - DeNit_Temp.BTC[0];

        CTimeSeries<double> Residual_normal_score = Residual.ConverttoNormalScore(); // Temperature residual

        Residual_normal_score.writefile(Workingfolder + "Data/Residual_normal_score.txt");

        CTimeSeries<double> Residual_autocorrelation = Residual_normal_score.AutoCorrelation(10,0.5);
        Residual_autocorrelation.writefile(Workingfolder + "Data/Residual_autocorrelation.txt");

        CTimeSeries<double> Residual_CDF = Residual.GetCummulativeDistribution();
        Residual_CDF.writefile(Workingfolder + "Data/Residual_CDF.txt");

        CTimeSeries<double> Residual_inv_cum = Residual_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> Residual_PDF = Residual.distribution(50,0);
        Residual_PDF.writefile(Workingfolder + "Data/Residual_PDF.txt");

        //double Residual_mean = Residual.Log().mean();
        //double Residual_std = Residual.Log().std();

        double Residual_autocorrelation_coeff = Residual_autocorrelation.AutoCorrelationCoeff();

        // Temperature
        CTimeSeries<double> OUP_Residual_NS;
        OUP_Residual_NS.CreateOUProcess(0,Simulation_time_Calc,dt,Residual_autocorrelation_coeff); // This is correct

        OUP_Residual_NS.writefile(Workingfolder + "OUP_Residual_NS_tvf.csv");

        CTimeSeries<double> OUP_Generated_Residuals = OUP_Residual_NS.MapfromNormalScoreToDistribution(Residual_inv_cum);
        CTimeSeries<double> Periodic_Temp = CTimeSeries<double>::CreateSinusoidal(0,Simulation_time_Calc,dt,sin_T0_p,sin_a,sin_b); // Sinusoidal

        OUP_Temp = OUP_Generated_Residuals + Periodic_Temp; // Residual + Sinusoidal

        OUP_Temp.writefile(Workingfolder + "OUP_Temp.csv");

    }


    //----------------------------------------------------------------------------------
    // Model Configuration

    // For Calibration
    if (Calibration)
    {
        Parameter mu_H_cal;
        mu_H_cal.SetQuantities(system,"Parameter");
        mu_H_cal.SetName("mu_H_cal");
        mu_H_cal.SetVal("value",v_mu_H);
        mu_H_cal.SetVal("high",v_mu_H_high);
        mu_H_cal.SetVal("low",v_mu_H_low);
        system->AppendParameter("mu_H_cal", mu_H_cal);

        Parameter b_H_cal;
        b_H_cal.SetQuantities(system,"Parameter");
        b_H_cal.SetName("b_H_cal");
        b_H_cal.SetVal("value",v_b_H);
        b_H_cal.SetVal("high",v_b_H_high);
        b_H_cal.SetVal("low",v_b_H_low);
        system->AppendParameter("b_H_cal", b_H_cal);

        Parameter K_S_cal;
        K_S_cal.SetQuantities(system,"Parameter");
        K_S_cal.SetName("K_S_cal");
        K_S_cal.SetVal("value",v_K_S);
        K_S_cal.SetVal("high",v_K_S_high);
        K_S_cal.SetVal("low",v_K_S_low);
        system->AppendParameter("K_S_cal", K_S_cal);

        Parameter K_MH_cal;
        K_MH_cal.SetQuantities(system,"Parameter");
        K_MH_cal.SetName("K_MH_cal");
        K_MH_cal.SetVal("value",v_K_MH);
        K_MH_cal.SetVal("high",v_K_MH_high);
        K_MH_cal.SetVal("low",v_K_MH_low);
        system->AppendParameter("K_MH_cal", K_MH_cal);

        Parameter K_OH_cal;
        K_OH_cal.SetQuantities(system,"Parameter");
        K_OH_cal.SetName("K_OH_cal");
        K_OH_cal.SetVal("value",v_K_OH);
        K_OH_cal.SetVal("high",v_K_OH_high);
        K_OH_cal.SetVal("low",v_K_OH_low);
        system->AppendParameter("K_OH_cal", K_OH_cal);

        Parameter K_NH_cal;
        K_NH_cal.SetQuantities(system,"Parameter");
        K_NH_cal.SetName("K_NH");
        K_NH_cal.SetVal("value",v_K_NH);
        K_NH_cal.SetVal("high",v_K_NH_high);
        K_NH_cal.SetVal("low",v_K_NH_low);
        system->AppendParameter("K_NH_cal", K_NH_cal);

        Parameter K_NOH_cal;
        K_NOH_cal.SetQuantities(system,"Parameter");
        K_NOH_cal.SetName("K_NOH_cal");
        K_NOH_cal.SetVal("value",v_K_NOH);
        K_NOH_cal.SetVal("high",v_K_NOH_high);
        K_NOH_cal.SetVal("low",v_K_NOH_low);
        system->AppendParameter("K_NOH_cal", K_NOH_cal);

        Parameter eta_g_cal;
        eta_g_cal.SetQuantities(system,"Parameter");
        eta_g_cal.SetName("eta_g_cal");
        eta_g_cal.SetVal("value",v_eta_g);
        eta_g_cal.SetVal("high",v_eta_g_high);
        eta_g_cal.SetVal("low",v_eta_g_low);
        system->AppendParameter("eta_g_cal", eta_g_cal);

        Parameter mu_A_cal;
        mu_A_cal.SetQuantities(system,"Parameter");
        mu_A_cal.SetName("mu_A_cal");
        mu_A_cal.SetVal("value",v_mu_A);
        mu_A_cal.SetVal("high",v_mu_A_high);
        mu_A_cal.SetVal("low",v_mu_A_low);
        system->AppendParameter("mu_A_cal", mu_A_cal);

        Parameter b_A_cal;
        b_A_cal.SetQuantities(system,"Parameter");
        b_A_cal.SetName("b_A_cal");
        b_A_cal.SetVal("value",v_b_A);
        b_A_cal.SetVal("high",v_b_A_high);
        b_A_cal.SetVal("low",v_b_A_low);
        system->AppendParameter("b_A_cal", b_A_cal);

        Parameter K_OA_cal;
        K_OA_cal.SetQuantities(system,"Parameter");
        K_OA_cal.SetName("K_OA_cal");
        K_OA_cal.SetVal("value",v_K_OA);
        K_OA_cal.SetVal("high",v_K_OA_high);
        K_OA_cal.SetVal("low",v_K_OA_low);
        system->AppendParameter("K_OA_cal", K_OA_cal);

        Parameter K_NHA_cal;
        K_NHA_cal.SetQuantities(system,"Parameter");
        K_NHA_cal.SetName("K_NHA_cal");
        K_NHA_cal.SetVal("value",v_K_NHA);
        K_NHA_cal.SetVal("high",v_K_NHA_high);
        K_NHA_cal.SetVal("low",v_K_NHA_low);
        system->AppendParameter("K_NHA_cal", K_NHA_cal);

        Parameter K_NOA_cal;
        K_NOA_cal.SetQuantities(system,"Parameter");
        K_NOA_cal.SetName("K_NOA_cal");
        K_NOA_cal.SetVal("value",v_K_NOA);
        K_NOA_cal.SetVal("high",v_K_NOA_high);
        K_NOA_cal.SetVal("low",v_K_NOA_low);
        system->AppendParameter("K_NOA_cal", K_NOA_cal);

        Parameter mu_M_cal;
        mu_M_cal.SetQuantities(system,"Parameter");
        mu_M_cal.SetName("mu_M_cal");
        mu_M_cal.SetVal("value",v_mu_M);
        mu_M_cal.SetVal("high",v_mu_M_high);
        mu_M_cal.SetVal("low",v_mu_M_low);
        system->AppendParameter("mu_M_cal", mu_M_cal);

        Parameter b_M_cal;
        b_M_cal.SetQuantities(system,"Parameter");
        b_M_cal.SetName("b_M_cal");
        b_M_cal.SetVal("value",v_b_M);
        b_M_cal.SetVal("high",v_b_M_high);
        b_M_cal.SetVal("low",v_b_M_low);
        system->AppendParameter("b_M_cal", b_M_cal);

        Parameter K_MM_cal;
        K_MM_cal.SetQuantities(system,"Parameter");
        K_MM_cal.SetName("K_MM_cal");
        K_MM_cal.SetVal("value",v_K_MM);
        K_MM_cal.SetVal("high",v_K_MM_high);
        K_MM_cal.SetVal("low",v_K_MM_low);
        system->AppendParameter("K_MM_cal", K_MM_cal);

        Parameter K_OM_cal;
        K_OM_cal.SetQuantities(system,"Parameter");
        K_OM_cal.SetName("K_OM_cal");
        K_OM_cal.SetVal("value",v_K_OM);
        K_OM_cal.SetVal("high",v_K_OM_high);
        K_OM_cal.SetVal("low",v_K_OM_low);
        system->AppendParameter("K_OM_cal", K_OM_cal);

        Parameter K_NOM_cal;
        K_NOM_cal.SetQuantities(system,"Parameter");
        K_NOM_cal.SetName("K_NOM");
        K_NOM_cal.SetVal("value",v_K_NOM);
        K_NOM_cal.SetVal("high",v_K_NOM_high);
        K_NOM_cal.SetVal("low",v_K_NOM_low);
        system->AppendParameter("K_NOM_cal", K_NOM_cal);

        Parameter K_a_cal;
        K_a_cal.SetQuantities(system,"Parameter");
        K_a_cal.SetName("K_a_cal");
        K_a_cal.SetVal("value",v_K_a);
        K_a_cal.SetVal("high",v_K_a_high);
        K_a_cal.SetVal("low",v_K_a_low);
        system->AppendParameter("K_a_cal", K_a_cal);

        Parameter K_h_cal;
        K_h_cal.SetQuantities(system,"Parameter");
        K_h_cal.SetName("K_h_cal");
        K_h_cal.SetVal("value",v_K_h);
        K_h_cal.SetVal("high",v_K_h_high);
        K_h_cal.SetVal("low",v_K_h_low);
        system->AppendParameter("K_h_cal", K_h_cal);

        Parameter K_X_cal;
        K_X_cal.SetQuantities(system,"Parameter");
        K_X_cal.SetName("K_X_cal");
        K_X_cal.SetVal("value",v_K_X);
        K_X_cal.SetVal("high",v_K_X_high);
        K_X_cal.SetVal("low",v_K_X_low);
        system->AppendParameter("K_X_cal", K_X_cal);

        Parameter eta_h_cal;
        eta_h_cal.SetQuantities(system,"Parameter");
        eta_h_cal.SetName("eta_h_cal");
        eta_h_cal.SetVal("value",v_eta_h);
        eta_h_cal.SetVal("high",v_eta_h_high);
        eta_h_cal.SetVal("low",v_eta_h_low);
        system->AppendParameter("eta_h_cal", eta_h_cal);

        Parameter i_XB_cal;
        i_XB_cal.SetQuantities(system,"Parameter");
        i_XB_cal.SetName("i_XB_cal");
        i_XB_cal.SetVal("value",v_i_XB);
        i_XB_cal.SetVal("high",v_i_XB_high);
        i_XB_cal.SetVal("low",v_i_XB_low);
        system->AppendParameter("i_XB_cal", i_XB_cal);

        Parameter K_LO2_cal;
        K_LO2_cal.SetQuantities(system,"Parameter");
        K_LO2_cal.SetName("K_LO2_cal");
        K_LO2_cal.SetVal("value",v_K_LO2);
        K_LO2_cal.SetVal("high",v_K_LO2_high);
        K_LO2_cal.SetVal("low",v_K_LO2_low);
        system->AppendParameter("K_LO2_cal", K_LO2_cal);

        Parameter Y_H_cal;
        Y_H_cal.SetQuantities(system,"Parameter");
        Y_H_cal.SetName("Y_H_cal");
        Y_H_cal.SetVal("value",v_Y_H);
        Y_H_cal.SetVal("high",v_Y_H_high);
        Y_H_cal.SetVal("low",v_Y_H_low);
        system->AppendParameter("Y_H_cal", Y_H_cal);

        Parameter Y_M_cal;
        Y_M_cal.SetQuantities(system,"Parameter");
        Y_M_cal.SetName("Y_M_cal");
        Y_M_cal.SetVal("value",v_Y_M);
        Y_M_cal.SetVal("high",v_Y_M_high);
        Y_M_cal.SetVal("low",v_Y_M_low);
        system->AppendParameter("Y_M_cal", Y_M_cal);
    }


    // Consistuents
    Constituent S_i;
    S_i.SetQuantities(system, "Constituent");
    S_i.SetName("S_i");
    S_i.SetType("Constituent");
    system->AddConstituent(S_i,false);

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
    mu_H.SetVal("Arrhenius_factor",v_mu_H_af);
    mu_H.SetVal("reference_temperature",v_rt);
    if (OUP)
    mu_H.SetProperty("temperature",Workingfolder + "OUP_Temp.csv");
    else
    mu_H.SetProperty("temperature",Workingfolder + "Data/DeNit_Temp.txt");
    system->AddReactionParameter(mu_H, false);
    if (Calibration)
    {
    system->SetAsParameter("mu_H","base_value","mu_H_cal");
    system->object("mu_H")->Variable("base_value")->SetParameterAssignedTo("mu_H_cal");
    }

    RxnParameter K_S;
    K_S.SetQuantities(system,"ReactionParameter");
    K_S.SetName("K_S");
    K_S.SetVal("base_value",v_K_S);
    K_S.SetVal("Arrhenius_factor",v_K_S_af);
    K_S.SetVal("reference_temperature",v_rt);
    if (OUP)
    K_S.SetProperty("temperature",Workingfolder + "OUP_Temp.csv");
    else
    K_S.SetProperty("temperature",Workingfolder + "Data/DeNit_Temp.txt");
    system->AddReactionParameter(K_S, false);
    if (Calibration)
    {
    system->SetAsParameter("K_S","base_value","K_S_cal");
    system->object("K_S")->Variable("base_value")->SetParameterAssignedTo("K_S_cal");
    }

    RxnParameter K_MH;
    K_MH.SetQuantities(system,"ReactionParameter");
    K_MH.SetName("K_MH");
    K_MH.SetVal("base_value",v_K_MH);
    system->AddReactionParameter(K_MH, false);
    if (Calibration)
    {
    system->SetAsParameter("K_MH","base_value","K_MH_cal");
    system->object("K_MH")->Variable("base_value")->SetParameterAssignedTo("K_MH_cal");
    }

    RxnParameter K_OH;
    K_OH.SetQuantities(system,"ReactionParameter");
    K_OH.SetName("K_OH");
    K_OH.SetVal("base_value",v_K_OH);
    system->AddReactionParameter(K_OH, false);
    if (Calibration)
    {
    system->SetAsParameter("K_OH","base_value","K_OH_cal");
    system->object("K_OH")->Variable("base_value")->SetParameterAssignedTo("K_OH_cal");
    }

    RxnParameter K_NOH;
    K_NOH.SetQuantities(system,"ReactionParameter");
    K_NOH.SetName("K_NOH");
    K_NOH.SetVal("base_value",v_K_NOH);
    system->AddReactionParameter(K_NOH, false);
    if (Calibration)
    {
    system->SetAsParameter("K_NOH","base_value","K_NOH_cal");
    system->object("K_NOH")->Variable("base_value")->SetParameterAssignedTo("K_NOH_cal");
    }

    RxnParameter eta_g;
    eta_g.SetQuantities(system,"ReactionParameter");
    eta_g.SetName("eta_g");
    eta_g.SetVal("base_value",v_eta_g);
    system->AddReactionParameter(eta_g, false);
    if (Calibration)
    {
    system->SetAsParameter("eta_g","base_value","eta_g_cal");
    system->object("eta_g")->Variable("base_value")->SetParameterAssignedTo("eta_g_cal");
    }

    RxnParameter b_H;
    b_H.SetQuantities(system,"ReactionParameter");
    b_H.SetName("b_H");
    b_H.SetVal("base_value",v_b_H);
    system->AddReactionParameter(b_H, false);
    if (Calibration)
    {
    system->SetAsParameter("b_H","base_value","b_H_cal");
    system->object("b_H")->Variable("base_value")->SetParameterAssignedTo("b_H_cal");
    }

    RxnParameter K_NH;
    K_NH.SetQuantities(system,"ReactionParameter");
    K_NH.SetName("K_NH");
    K_NH.SetVal("base_value",v_K_NH);
    system->AddReactionParameter(K_NH, false);
    if (Calibration)
    {
    system->SetAsParameter("K_NH","base_value","K_NH_cal");
    system->object("K_NH")->Variable("base_value")->SetParameterAssignedTo("K_NH_cal");
    }

    RxnParameter mu_M;
    mu_M.SetQuantities(system,"ReactionParameter");
    mu_M.SetName("mu_M");
    mu_M.SetVal("base_value",v_mu_M);
    mu_M.SetVal("Arrhenius_factor",v_mu_M_af);
    mu_M.SetVal("reference_temperature",v_rt);
    if (OUP)
    mu_M.SetProperty("temperature",Workingfolder + "OUP_Temp.csv");
    else
    mu_M.SetProperty("temperature",Workingfolder + "Data/DeNit_Temp.txt");
    system->AddReactionParameter(mu_M, false);
    if (Calibration)
    {
    system->SetAsParameter("mu_M","base_value","mu_M_cal");
    system->object("mu_M")->Variable("base_value")->SetParameterAssignedTo("mu_M_cal");
    }

    RxnParameter K_MM;
    K_MM.SetQuantities(system,"ReactionParameter");
    K_MM.SetName("K_MM");
    K_MM.SetVal("base_value",v_K_MM);
    system->AddReactionParameter(K_MM, false);
    if (Calibration)
    {
    system->SetAsParameter("K_MM","base_value","K_MM_cal");
    system->object("K_MM")->Variable("base_value")->SetParameterAssignedTo("K_MM_cal");
    }

    RxnParameter K_OM;
    K_OM.SetQuantities(system,"ReactionParameter");
    K_OM.SetName("K_OM");
    K_OM.SetVal("base_value",v_K_OM);
    system->AddReactionParameter(K_OM, false);
    if (Calibration)
    {
    system->SetAsParameter("K_OM","base_value","K_OM_cal");
    system->object("K_OM")->Variable("base_value")->SetParameterAssignedTo("K_OM_cal");
    }

    RxnParameter K_NOM;
    K_NOM.SetQuantities(system,"ReactionParameter");
    K_NOM.SetName("K_NOM");
    K_NOM.SetVal("base_value",v_K_NOM);
    system->AddReactionParameter(K_NOM, false);
    if (Calibration)
    {
    system->SetAsParameter("K_NOM","base_value","K_NOM_cal");
    system->object("K_NOM")->Variable("base_value")->SetParameterAssignedTo("K_NOM_cal");
    }

    RxnParameter b_M;
    b_M.SetQuantities(system,"ReactionParameter");
    b_M.SetName("b_M");
    b_M.SetVal("base_value",v_b_M);
    b_M.SetVal("Arrhenius_factor",v_b_M_af);
    b_M.SetVal("reference_temperature",v_rt);
    if (OUP)
    b_M.SetProperty("temperature",Workingfolder + "OUP_Temp.csv");
    else
    b_M.SetProperty("temperature",Workingfolder + "Data/DeNit_Temp.txt");
    system->AddReactionParameter(b_M, false);
    if (Calibration)
    {
    system->SetAsParameter("b_M","base_value","b_M_cal");
    system->object("b_M")->Variable("base_value")->SetParameterAssignedTo("b_M_cal");
    }

    RxnParameter mu_A;
    mu_A.SetQuantities(system,"ReactionParameter");
    mu_A.SetName("mu_A");
    mu_A.SetVal("base_value",v_mu_A);
    mu_A.SetVal("Arrhenius_factor",v_mu_A_af);
    mu_A.SetVal("reference_temperature",v_rt);
    if (OUP)
    mu_A.SetProperty("temperature",Workingfolder + "OUP_Temp.csv");
    else
    mu_A.SetProperty("temperature",Workingfolder + "Data/DeNit_Temp.txt");
    system->AddReactionParameter(mu_A, false);
    if (Calibration)
    {
    system->SetAsParameter("mu_A","base_value","mu_A_cal");
    system->object("mu_A")->Variable("base_value")->SetParameterAssignedTo("mu_A_cal");
    }

    RxnParameter K_NHA;
    K_NHA.SetQuantities(system,"ReactionParameter");
    K_NHA.SetName("K_NHA");
    K_NHA.SetVal("base_value",v_K_NHA);
    system->AddReactionParameter(K_NHA, false);
    if (Calibration)
    {
    system->SetAsParameter("K_NHA","base_value","K_NHA_cal");
    system->object("K_NHA")->Variable("base_value")->SetParameterAssignedTo("K_NHA_cal");
    }

    RxnParameter K_NOA;
    K_NOA.SetQuantities(system,"ReactionParameter");
    K_NOA.SetName("K_NOA");
    K_NOA.SetVal("base_value",v_K_NOA);
    system->AddReactionParameter(K_NOA, false);
    if (Calibration)
    {
    system->SetAsParameter("K_NOA","base_value","K_NOA_cal");
    system->object("K_NOA")->Variable("base_value")->SetParameterAssignedTo("K_NOA_cal");
    }

    RxnParameter K_OA;
    K_OA.SetQuantities(system,"ReactionParameter");
    K_OA.SetName("K_OA");
    K_OA.SetVal("base_value",v_K_OA);
    system->AddReactionParameter(K_OA, false);
    if (Calibration)
    {
    system->SetAsParameter("K_OA","base_value","K_OA_cal");
    system->object("K_OA")->Variable("base_value")->SetParameterAssignedTo("K_OA_cal");
    }

    RxnParameter b_A;
    b_A.SetQuantities(system,"ReactionParameter");
    b_A.SetName("b_A");
    b_A.SetVal("base_value",v_b_A);
    b_A.SetVal("Arrhenius_factor",v_b_A_af);
    b_A.SetVal("reference_temperature",v_rt);
    if (OUP)
    b_A.SetProperty("temperature",Workingfolder + "OUP_Temp.csv");
    else
    b_A.SetProperty("temperature",Workingfolder + "Data/DeNit_Temp.txt");
    system->AddReactionParameter(b_A, false);
    if (Calibration)
    {
    system->SetAsParameter("b_A","base_value","b_A_cal");
    system->object("b_A")->Variable("base_value")->SetParameterAssignedTo("b_A_cal");
    }

    RxnParameter eta_h;
    eta_h.SetQuantities(system,"ReactionParameter");
    eta_h.SetName("eta_h");
    eta_h.SetVal("base_value",v_eta_h);
    system->AddReactionParameter(eta_h, false);
    if (Calibration)
    {
    system->SetAsParameter("eta_h","base_value","eta_h_cal");
    system->object("eta_h")->Variable("base_value")->SetParameterAssignedTo("eta_h_cal");
    }

    RxnParameter K_h;
    K_h.SetQuantities(system,"ReactionParameter");
    K_h.SetName("K_h");
    K_h.SetVal("base_value",v_K_h);
    K_h.SetVal("Arrhenius_factor",v_K_h_af);
    K_h.SetVal("reference_temperature",v_rt);
    if (OUP)
    K_h.SetProperty("temperature",Workingfolder + "OUP_Temp.csv");
    else
    K_h.SetProperty("temperature",Workingfolder + "Data/DeNit_Temp.txt");
    system->AddReactionParameter(K_h, false);
    if (Calibration)
    {
    system->SetAsParameter("K_h","base_value","K_h_cal");
    system->object("K_h")->Variable("base_value")->SetParameterAssignedTo("K_h_cal");
    }

    RxnParameter K_X;
    K_X.SetQuantities(system,"ReactionParameter");
    K_X.SetName("K_X");
    K_X.SetVal("base_value",v_K_X);
    system->AddReactionParameter(K_X, false);
    if (Calibration)
    {
    system->SetAsParameter("K_X","base_value","K_X_cal");
    system->object("K_X")->Variable("base_value")->SetParameterAssignedTo("K_X_cal");
    }

    RxnParameter K_a;
    K_a.SetQuantities(system,"ReactionParameter");
    K_a.SetName("K_a");
    K_a.SetVal("base_value",v_K_a);
    K_a.SetVal("Arrhenius_factor",v_K_a_af);
    K_a.SetVal("reference_temperature",v_rt);
    if (OUP)
    K_a.SetProperty("temperature",Workingfolder + "OUP_Temp.csv");
    else
    K_a.SetProperty("temperature",Workingfolder + "Data/DeNit_Temp.txt");
    system->AddReactionParameter(K_a, false);
    if (Calibration)
    {
    system->SetAsParameter("K_a","base_value","K_a_cal");
    system->object("K_a")->Variable("base_value")->SetParameterAssignedTo("K_a_cal");
    }

    RxnParameter Y_H;
    Y_H.SetQuantities(system,"ReactionParameter");
    Y_H.SetName("Y_H");
    Y_H.SetVal("base_value",v_Y_H);
    system->AddReactionParameter(Y_H, false);
    if (Calibration)
    {
    system->SetAsParameter("Y_H","base_value","Y_H_cal");
    system->object("Y_H")->Variable("base_value")->SetParameterAssignedTo("Y_H_cal");
    }

    RxnParameter Y_HM;
    Y_HM.SetQuantities(system,"ReactionParameter");
    Y_HM.SetName("Y_HM");
    Y_HM.SetVal("base_value",v_Y_HM);
    system->AddReactionParameter(Y_HM, false);
    //if (Calibration)
    //{
    //system->SetAsParameter("Y_HM","base_value","Y_HM_cal");
    //system->object("Y_HM")->Variable("base_value")->SetParameterAssignedTo("Y_HM_cal");
    //}

    RxnParameter Y_A;
    Y_A.SetQuantities(system,"ReactionParameter");
    Y_A.SetName("Y_A");
    Y_A.SetVal("base_value",v_Y_A);
    system->AddReactionParameter(Y_A, false);
    //if (Calibration)
    //{
    //system->SetAsParameter("Y_A","base_value","Y_A_cal");
    //system->object("Y_A")->Variable("base_value")->SetParameterAssignedTo("Y_A_cal");
    //}

    RxnParameter Y_M;
    Y_M.SetQuantities(system,"ReactionParameter");
    Y_M.SetName("Y_M");
    Y_M.SetVal("base_value",v_Y_M);
    system->AddReactionParameter(Y_M, false);
    if (Calibration)
    {
    system->SetAsParameter("Y_M","base_value","Y_M_cal");
    system->object("Y_M")->Variable("base_value")->SetParameterAssignedTo("Y_M_cal");
    }

    RxnParameter f_p;
    f_p.SetQuantities(system,"ReactionParameter");
    f_p.SetName("f_p");
    f_p.SetVal("base_value",v_f_p);
    system->AddReactionParameter(f_p, false);
    //if (Calibration)
    //{
    //system->SetAsParameter("f_p","base_value","f_p_cal");
    //system->object("f_p")->Variable("base_value")->SetParameterAssignedTo("f_p_cal");
    //}

    RxnParameter i_XB;
    i_XB.SetQuantities(system,"ReactionParameter");
    i_XB.SetName("i_XB");
    i_XB.SetVal("base_value",v_i_XB);
    system->AddReactionParameter(i_XB, false);
    if (Calibration)
    {
    system->SetAsParameter("i_XB","base_value","i_XB_cal");
    system->object("i_XB")->Variable("base_value")->SetParameterAssignedTo("i_XB_cal");
    }

    RxnParameter i_VSSB;
    i_VSSB.SetQuantities(system,"ReactionParameter");
    i_VSSB.SetName("i_VSSB");
    i_VSSB.SetVal("base_value",v_i_VSSB);
    system->AddReactionParameter(i_VSSB, false);
    //if (Calibration)
    //{
    //system->SetAsParameter("VSSB","base_value","VSSB_cal");
    //system->object("VSSB")->Variable("base_value")->SetParameterAssignedTo("VSSB_cal");
    //}

    RxnParameter i_VSSi;
    i_VSSi.SetQuantities(system,"ReactionParameter");
    i_VSSi.SetName("i_VSSi");
    i_VSSi.SetVal("base_value",v_i_VSSi);
    system->AddReactionParameter(i_VSSi, false);
    //if (Calibration)
    //{
    //system->SetAsParameter("VSSi","base_value","VSSi_cal");
    //system->object("VSSi")->Variable("base_value")->SetParameterAssignedTo("VSSi_cal");
    //}

    RxnParameter i_VSSs;
    i_VSSs.SetQuantities(system,"ReactionParameter");
    i_VSSs.SetName("i_VSSs");
    i_VSSs.SetVal("base_value",v_i_VSSs);
    system->AddReactionParameter(i_VSSs, false);
    //if (Calibration)
    //{
    //system->SetAsParameter("VSSs","base_value","VSSs_cal");
    //system->object("VSSs")->Variable("base_value")->SetParameterAssignedTo("VSSs_cal");
    //}

    RxnParameter i_VSSP;
    i_VSSP.SetQuantities(system,"ReactionParameter");
    i_VSSP.SetName("i_VSSP");
    i_VSSP.SetVal("base_value",v_i_VSSP);
    system->AddReactionParameter(i_VSSP, false);
    //if (Calibration)
    //{
    //system->SetAsParameter("VSSP","base_value","VSSP_cal");
    //system->object("VSSP")->Variable("base_value")->SetParameterAssignedTo("VSSP_cal");
    //}

    RxnParameter i_MeOH;
    i_MeOH.SetQuantities(system,"ReactionParameter");
    i_MeOH.SetName("i_MeOH");
    i_MeOH.SetVal("base_value",v_i_MeOH);
    system->AddReactionParameter(i_MeOH, false);
    //if (Calibration)
    //{
    //system->SetAsParameter("MeOH","base_value","MeOH_cal");
    //system->object("MeOH")->Variable("base_value")->SetParameterAssignedTo("MeOH_cal");
    //}



    // Sources
    Source Aeration;
    Aeration.SetQuantities(system, "atmospheric exchange");
    Aeration.SetType("atmospheric exchange");
    Aeration.SetName("Aeration");
    Aeration.SetVal("rate_coefficient",v_K_LO2);
    Aeration.SetVal("saturation",v_a_saturation);
    system->AddSource(Aeration, false);
    if (Calibration)
    {
    // Setting as parameter for calibration
    system->SetAsParameter("Aeration","rate_coefficient","K_LO2_cal");
    system->object("Aeration")->Variable("rate_coefficient")->SetParameterAssignedTo("K_LO2_cal");
    }

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
    Stl_element_top.SetVal("Settled_Particles:concentration",0);
    Stl_element_top.SetVal("Solids:concentration",0);
    Stl_element_top.SetVal("Storage",v_s_t_storage);
    Stl_element_top.SetVal("bottom_elevation",v_s_t_bottom_elevation);
    Stl_element_top.SetVal("Volume",v_s_t_volume);
    Stl_element_top.SetVal("x",600);
    Stl_element_top.SetVal("y",500);
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
    Stl_element_bottom.SetVal("x",600);
    Stl_element_bottom.SetVal("y",1100);
    system->AddBlock(Stl_element_bottom,false);

    // Time-variable Fixed Head Blocks
    Block fh_clarifier;
    fh_clarifier.SetQuantities(system, "time_variable_fixed_head");
    fh_clarifier.SetName("Clarifier");
    fh_clarifier.SetType("time_variable_fixed_head");
    fh_clarifier.SetVal("head",v_c_bottom_elevation);
    fh_clarifier.SetVal("S_S:concentration",0);
    fh_clarifier.SetVal("X_b:concentration",0);
    fh_clarifier.SetVal("Storage",100000);
    fh_clarifier.SetProperty("Dummy_timeseries",Workingfolder + "OUP_Temp.csv");
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
    fh_WAS.SetProperty("Dummy_timeseries",Workingfolder + "OUP_Flow_WAS_tvf.csv");
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
        //if (i==4) Reactor_Flex.SetVal("S_M:External mass flow time-series",Workingfolder + "S_M_mfr.csv");
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

    CTimeSeries<double> S_i_inflow_concentration;
    S_i_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_i_concentration);
    S_i_inflow_concentration.writefile(Workingfolder + "S_i_constant_inflow_concentration.txt");

    CTimeSeries<double> S_S_inflow_concentration;
    S_S_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_S_concentration);
    S_S_inflow_concentration.writefile(Workingfolder + "S_S_constant_inflow_concentration.txt");

    CTimeSeries<double> X_S_inflow_concentration;
    X_S_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_X_S_concentration);
    X_S_inflow_concentration.writefile(Workingfolder + "X_S_constant_inflow_concentration.txt");

    CTimeSeries<double> X_p_inflow_concentration;
    X_p_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_X_p_concentration);
    X_p_inflow_concentration.writefile(Workingfolder + "X_p_constant_inflow_concentration.txt");

    CTimeSeries<double> S_NO_inflow_concentration;
    S_NO_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_NO_concentration);
    S_NO_inflow_concentration.writefile(Workingfolder + "S_NO_constant_inflow_concentration.txt");

    CTimeSeries<double> S_NH_inflow_concentration;
    S_NH_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_NH_concentration);
    S_NH_inflow_concentration.writefile(Workingfolder + "S_NH_constant_inflow_concentration.txt");

    CTimeSeries<double> S_ND_inflow_concentration;
    S_ND_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_S_ND_concentration);
    S_ND_inflow_concentration.writefile(Workingfolder + "S_ND_constant_inflow_concentration.txt");

    CTimeSeries<double> X_ND_inflow_concentration;
    X_ND_inflow_concentration.CreateConstant(0,Simulation_time_Calc, v_X_ND_concentration);
    X_ND_inflow_concentration.writefile(Workingfolder + "X_ND_constant_inflow_concentration.txt");


    // Constituents Inflow Calculations from Data (Time Variable)


    // Determining Coefficients

    const double c_S_i=p_31; // S_i, 8
    const double c_S_S=1-p_31; // S_S, 8
    const double c_X_S=0.71*p_32*(1-p_2); // X_S, 2
    const double c_X_p=0.71*p_32*p_2; // X_p, 2
    const double c_S_NO=1; // S_NO, 10
    const double c_S_NH=1; // S_NH, 9
    const double c_S_ND=p_36*(1-p_31); // S_ND, 8
    const double c_X_ND=p_3*p_32; // X_ND, 2

    if (!OUP)

    {

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
    // Denite Parameters(0-10): Q, TSS, VSS, BOD, sBOD, TP, TSP, TKN, sCOD, NH3, NO3

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

    Inflow_Q.writefile(Workingfolder + "Q_tvif.txt"); //

    Flow_WAS.writefile(Workingfolder + "WAS_tvf.txt"); //
    Flow_RAS.writefile(Workingfolder + "RAS_tvf.txt"); //

    Flow_r_r_st.writefile(Workingfolder + "r_r_st_tvf.txt"); //
    Flow_st_sb.writefile(Workingfolder + "st_sb_tvf.txt"); //
    Flow_st_c.writefile(Workingfolder + "st_c_tvf.txt"); //

    Inflow_S_i.writefile(Workingfolder + "S_i_Inflow_concentration.txt");
    Inflow_S_S.writefile(Workingfolder + "S_S_Inflow_concentration.txt");
    Inflow_X_S.writefile(Workingfolder + "X_S_Inflow_concentration.txt");
    Inflow_X_p.writefile(Workingfolder + "X_p_Inflow_concentration.txt");
    Inflow_S_NO.writefile(Workingfolder + "S_NO_Inflow_concentration.txt");
    Inflow_S_NH.writefile(Workingfolder + "S_NH_Inflow_concentration.txt");
    Inflow_S_ND.writefile(Workingfolder + "S_ND_Inflow_concentration.txt");
    Inflow_X_ND.writefile(Workingfolder + "X_ND_Inflow_concentration.txt");

    Inflow_MeOH.writefile(Workingfolder + "S_M_mfr.txt"); // Methanol

    // Assigning Constituents Inflow concentrations to the system (Reactor_Flex)

    // Constant inflows by given concentration
/*
    system->block("Reactor_Flex")->SetProperty("S_S:constant_inflow_concentration",Workingfolder + "S_S_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("X_S:constant_inflow_concentration",Workingfolder + "X_S_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("X_p:constant_inflow_concentration",Workingfolder + "X_p_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("S_NO:constant_inflow_concentration",Workingfolder + "S_NO_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("S_NH:constant_inflow_concentration",Workingfolder + "S_NH_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("S_ND:constant_inflow_concentration",Workingfolder + "S_ND_constant_inflow_concentration.txt");
    system->block("Reactor_Flex")->SetProperty("X_ND:constant_inflow_concentration",Workingfolder + "X_ND_constant_inflow_concentration.txt");
*/


    // Time Variable inflows according DeNite data

    system->block("Reactor_Flex(1)")->SetProperty("time_variable_inflow",Workingfolder + "Q_tvif.txt"); // Discharge (m3/day)

    system->block("Reactor_Flex(1)")->SetProperty("S_i:time_variable_inflow_concentration",Workingfolder + "S_i_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("S_S:time_variable_inflow_concentration",Workingfolder + "S_S_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("X_S:time_variable_inflow_concentration",Workingfolder + "X_S_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("X_p:time_variable_inflow_concentration",Workingfolder + "X_p_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("S_NO:time_variable_inflow_concentration",Workingfolder + "S_NO_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("S_NH:time_variable_inflow_concentration",Workingfolder + "S_NH_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("S_ND:time_variable_inflow_concentration",Workingfolder + "S_ND_Inflow_concentration.txt");
    system->block("Reactor_Flex(1)")->SetProperty("X_ND:time_variable_inflow_concentration",Workingfolder + "X_ND_Inflow_concentration.txt");

    system->block("Reactor_Flex(5)")->SetProperty("S_M:external_mass_flow_timeseries",Workingfolder + "S_M_mfr.txt");

    // Flex_flow Links for Reactor_Flex
    for (int i=0; i<n_tanks-1; i++)
    {   Link l_r_st;
        l_r_st.SetQuantities(system, "Flex_flow");
        l_r_st.SetName("Reactor_Flex(" + aquiutils::numbertostring(i+1)+"_" + aquiutils::numbertostring(i+2) + ")");
        l_r_st.SetType("Flex_flow");
        //l_r_st.SetProperty("flow", Workingfolder + "r_r_st_tvf.txt");
        l_r_st.SetVal("flow_factor",v_flow_factor_i);
        system->AddLink(l_r_st, "Reactor_Flex(" + aquiutils::numbertostring(i+1)+")", "Reactor_Flex(" + aquiutils::numbertostring(i+2)+")", false);
    }


    // Flex_flow Links
    Link l_r_st;
    l_r_st.SetQuantities(system, "Flex_flow");
    l_r_st.SetName("Reactor_Flex(" + aquiutils::numbertostring(n_tanks) + ") - Settling element top");
    l_r_st.SetType("Flex_flow");
    //l_r_st.SetProperty("flow", Workingfolder + "r_r_st_tvf.txt");
    l_r_st.SetVal("flow_factor",v_flow_factor_i);
    system->AddLink(l_r_st, "Reactor_Flex(" + aquiutils::numbertostring(n_tanks) + ")", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Flex_flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Flex_flow");
    //l_st_c.SetProperty("flow", Workingfolder + "st_c_tvf.txt");
    l_st_c.SetVal("flow_factor",v_flow_factor_o);
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_sb_was;
    l_sb_was.SetQuantities(system, "Flex_flow");
    l_sb_was.SetName("Settling element bottom - WAS");
    l_sb_was.SetType("Flex_flow");
    //l_sb_was.SetProperty("flow", Workingfolder + "WAS_tvf.txt");
    l_sb_was.SetVal("flow_factor",v_flow_factor_o);
    system->AddLink(l_sb_was, "Settling element bottom", "WAS", false);

    // Time-Dependent interface Link
    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface (time variable)");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface (time variable)");
    l_st_sb.SetProperty("flow", Workingfolder + "st_sb_tvf.txt");
    l_st_sb.SetVal("area", v_st_sb_area);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    // Time-Dependent flow Links
    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Time-Dependent flow");
    l_sb_r.SetName("Settling element bottom - Reactor_Flex(1)");
    l_sb_r.SetType("Time-Dependent flow");
    l_sb_r.SetProperty("flow", Workingfolder + "RAS_tvf.txt");
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor_Flex(1)", false);

    }

    if (OUP)

    {
    // Determining Inflows by OUProcess

        // Inflows

        CTimeSeriesSet<double> Inflow_DeNit(Workingfolder + "Data/DeNit_Influent_Lump.txt",true);
        CTimeSeriesSet<double> Inflow_DeNit_wasteflow(Workingfolder + "Data/DeNit_wasteflow.txt",true);
        CTimeSeriesSet<double> Inflow_DeNit_returnflow(Workingfolder + "Data/DeNit_returnflow.txt",true);
        CTimeSeriesSet<double> Inflow_DeNit_MeOH(Workingfolder + "Data/DeNit_MeOH.txt",true);


        CTimeSeries<double> OUP_Inflow_Q; // Discharge (m3/day)

        CTimeSeries<double> OUP_Flow_WAS; // Wasteflow (WAS) (m3/day)
        CTimeSeries<double> OUP_Flow_RAS; // Returnflow (RAS) (m3/day)

        CTimeSeries<double> OUP_Flow_r_r_st; // Reactor_Flex to Reactor_Flex + Reactor_Flex to Settling element top flow (m3/day)
        CTimeSeries<double> OUP_Flow_st_sb; // Settling element top to Settling element bottom flow (m3/day)
        CTimeSeries<double> OUP_Flow_st_c; // Settling element top to Clarifier flow (m3/day)

        CTimeSeries<double> OUP_TSS_Flow; // TSS (1)
        CTimeSeries<double> OUP_VSS_Flow; // VSS (2)
        CTimeSeries<double> OUP_BOD_Flow; // BOD (3)
        CTimeSeries<double> OUP_sCOD_Flow; // sCOD (8)
        CTimeSeries<double> OUP_NH3_Flow; // NH3 (9)
        CTimeSeries<double> OUP_NO3_Flow; // NO3 (10)

        CTimeSeries<double> OUP_MeOH_Flow; // MeOH (Methanol)

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

        //Data analysis       

        //Flow (0)
        CTimeSeries<double> flow_normal_score = Inflow_DeNit.BTC[0].ConverttoNormalScore();
        flow_normal_score.writefile(Workingfolder + "Data/flow_normal_score.txt");

        CTimeSeries<double> flow_autocorrelation = flow_normal_score.AutoCorrelation(10,0.5);
        flow_autocorrelation.writefile(Workingfolder + "Data/flow_autocorrelation.txt");

        CTimeSeries<double> flow_CDF = Inflow_DeNit.BTC[0].GetCummulativeDistribution();
        flow_CDF.writefile(Workingfolder + "Data/flow_CDF.txt");

        CTimeSeries<double> flow_inv_cum = flow_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> flow_PDF = Inflow_DeNit.BTC[0].distribution(50,0);
        flow_PDF.writefile(Workingfolder + "Data/flow_PDF.txt");

        double flow_mean = Inflow_DeNit.BTC[0].Log().mean();
        double flow_std = Inflow_DeNit.BTC[0].Log().std();
        double flow_autocorrelation_coeff = flow_autocorrelation.AutoCorrelationCoeff();

        //TSS (1)
        CTimeSeries<double> TSS_normal_score = Inflow_DeNit.BTC[1].ConverttoNormalScore();
        TSS_normal_score.writefile(Workingfolder + "Data/TSS_normal_score.txt");

        CTimeSeries<double> TSS_autocorrelation = TSS_normal_score.AutoCorrelation(10,0.5);
        TSS_autocorrelation.writefile(Workingfolder + "Data/TSS_autocorrelation.txt");

        CTimeSeries<double> TSS_CDF = Inflow_DeNit.BTC[1].GetCummulativeDistribution();
        TSS_CDF.writefile(Workingfolder + "Data/TSS_CDF.txt");

        CTimeSeries<double> TSS_inv_cum = TSS_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> TSS_PDF = Inflow_DeNit.BTC[1].distribution(50,0);
        TSS_PDF.writefile(Workingfolder + "Data/TSS_PDF.txt");

        double TSS_mean = Inflow_DeNit.BTC[1].Log().mean();
        double TSS_std = Inflow_DeNit.BTC[1].Log().std();
        double TSS_autocorrelation_coeff = TSS_autocorrelation.AutoCorrelationCoeff();

        //VSS (2)
        CTimeSeries<double> VSS_normal_score = Inflow_DeNit.BTC[2].ConverttoNormalScore();
        VSS_normal_score.writefile(Workingfolder + "Data/VSS_normal_score.txt");

        CTimeSeries<double> VSS_autocorrelation = VSS_normal_score.AutoCorrelation(10,0.5);
        VSS_autocorrelation.writefile(Workingfolder + "Data/VSS_autocorrelation.txt");

        CTimeSeries<double> VSS_CDF = Inflow_DeNit.BTC[2].GetCummulativeDistribution();
        VSS_CDF.writefile(Workingfolder + "Data/VSS_CDF.txt");

        CTimeSeries<double> VSS_inv_cum = VSS_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> VSS_PDF = Inflow_DeNit.BTC[2].distribution(50,0);
        VSS_PDF.writefile(Workingfolder + "Data/VSS_PDF.txt");

        double VSS_mean = Inflow_DeNit.BTC[2].Log().mean();
        double VSS_std = Inflow_DeNit.BTC[2].Log().std();
        double VSS_autocorrelation_coeff = VSS_autocorrelation.AutoCorrelationCoeff();

        //BOD (3)
        CTimeSeries<double> BOD_normal_score = Inflow_DeNit.BTC[3].ConverttoNormalScore();
        BOD_normal_score.writefile(Workingfolder + "Data/BOD_normal_score.txt");

        CTimeSeries<double> BOD_autocorrelation = BOD_normal_score.AutoCorrelation(10,0.5);
        BOD_autocorrelation.writefile(Workingfolder + "Data/BOD_autocorrelation.txt");

        CTimeSeries<double> BOD_CDF = Inflow_DeNit.BTC[3].GetCummulativeDistribution();
        BOD_CDF.writefile(Workingfolder + "Data/BOD_CDF.txt");

        CTimeSeries<double> BOD_inv_cum = BOD_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> BOD_PDF = Inflow_DeNit.BTC[3].distribution(50,0);
        BOD_PDF.writefile(Workingfolder + "Data/BOD_PDF.txt");

        double BOD_mean = Inflow_DeNit.BTC[3].Log().mean();
        double BOD_std = Inflow_DeNit.BTC[3].Log().std();
        double BOD_autocorrelation_coeff = BOD_autocorrelation.AutoCorrelationCoeff();

        //sCOD (8)
        CTimeSeries<double> sCOD_normal_score = Inflow_DeNit.BTC[8].ConverttoNormalScore();
        sCOD_normal_score.writefile(Workingfolder + "Data/sCOD_normal_score.txt");

        CTimeSeries<double> sCOD_autocorrelation = sCOD_normal_score.AutoCorrelation(10,0.5);
        sCOD_autocorrelation.writefile(Workingfolder + "Data/sCOD_autocorrelation.txt");

        CTimeSeries<double> sCOD_CDF = Inflow_DeNit.BTC[8].GetCummulativeDistribution();
        sCOD_CDF.writefile(Workingfolder + "Data/sCOD_CDF.txt");

        CTimeSeries<double> sCOD_inv_cum = sCOD_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> sCOD_PDF = Inflow_DeNit.BTC[8].distribution(50,0);
        sCOD_PDF.writefile(Workingfolder + "Data/sCOD_PDF.txt");

        double sCOD_mean = Inflow_DeNit.BTC[8].Log().mean();
        double sCOD_std = Inflow_DeNit.BTC[8].Log().std();
        double sCOD_autocorrelation_coeff = sCOD_autocorrelation.AutoCorrelationCoeff();

        //NH3 (9)
        CTimeSeries<double> NH3_normal_score = Inflow_DeNit.BTC[9].ConverttoNormalScore();
        NH3_normal_score.writefile(Workingfolder + "Data/NH3_normal_score.txt");

        CTimeSeries<double> NH3_autocorrelation = NH3_normal_score.AutoCorrelation(10,0.5);
        NH3_autocorrelation.writefile(Workingfolder + "Data/NH3_autocorrelation.txt");

        CTimeSeries<double> NH3_CDF = Inflow_DeNit.BTC[9].GetCummulativeDistribution();
        NH3_CDF.writefile(Workingfolder + "Data/NH3_CDF.txt");

        CTimeSeries<double> NH3_inv_cum = NH3_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> NH3_PDF = Inflow_DeNit.BTC[9].distribution(50,0);
        NH3_PDF.writefile(Workingfolder + "Data/NH3_PDF.txt");

        double NH3_mean = Inflow_DeNit.BTC[9].Log().mean();
        double NH3_std = Inflow_DeNit.BTC[9].Log().std();
        double NH3_autocorrelation_coeff = NH3_autocorrelation.AutoCorrelationCoeff();

        //NO3 (10)
        CTimeSeries<double> NO3_normal_score = Inflow_DeNit.BTC[10].ConverttoNormalScore();
        NO3_normal_score.writefile(Workingfolder + "Data/NO3_normal_score.txt");

        CTimeSeries<double> NO3_autocorrelation = NO3_normal_score.AutoCorrelation(10,0.5);
        NO3_autocorrelation.writefile(Workingfolder + "Data/NO3_autocorrelation.txt");

        CTimeSeries<double> NO3_CDF = Inflow_DeNit.BTC[10].GetCummulativeDistribution();
        NO3_CDF.writefile(Workingfolder + "Data/NO3_CDF.txt");

        CTimeSeries<double> NO3_inv_cum = NO3_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> NO3_PDF = Inflow_DeNit.BTC[10].distribution(50,0);
        NO3_PDF.writefile(Workingfolder + "Data/NO3_PDF.txt");

        double NO3_mean = Inflow_DeNit.BTC[10].Log().mean();
        double NO3_std = Inflow_DeNit.BTC[10].Log().std();
        double NO3_autocorrelation_coeff = NO3_autocorrelation.AutoCorrelationCoeff();

        //WAS_Flow (0)
        CTimeSeries<double> was_flow_normal_score = Inflow_DeNit_wasteflow.BTC[0].ConverttoNormalScore();
        was_flow_normal_score.writefile(Workingfolder + "Data/was_flow_normal_score.txt");

        CTimeSeries<double> was_flow_autocorrelation = was_flow_normal_score.AutoCorrelation(10,0.5);
        was_flow_autocorrelation.writefile(Workingfolder + "Data/was_flow_autocorrelation.txt");

        CTimeSeries<double> was_flow_CDF = Inflow_DeNit_wasteflow.BTC[0].GetCummulativeDistribution();
        was_flow_CDF.writefile(Workingfolder + "Data/was_flow_CDF.txt");

        CTimeSeries<double> was_flow_inv_cum = was_flow_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> was_flow_PDF = Inflow_DeNit_wasteflow.BTC[0].distribution(50,0);
        was_flow_PDF.writefile(Workingfolder + "Data/was_flow_PDF.txt");

        double was_flow_mean = Inflow_DeNit_wasteflow.BTC[0].Log(1).mean(); // Becuase of Zeros
        double was_flow_std = Inflow_DeNit_wasteflow.BTC[0].Log(1).std(); // Becuase of Zeros
        double was_flow_autocorrelation_coeff = was_flow_autocorrelation.AutoCorrelationCoeff();

        //RAS_Flow (0)
        CTimeSeries<double> ras_flow_normal_score = Inflow_DeNit_returnflow.BTC[0].ConverttoNormalScore();
        ras_flow_normal_score.writefile(Workingfolder + "Data/ras_flow_normal_score.txt");

        CTimeSeries<double> ras_flow_autocorrelation = ras_flow_normal_score.AutoCorrelation(10,0.5);
        ras_flow_autocorrelation.writefile(Workingfolder + "Data/ras_flow_autocorrelation.txt");

        CTimeSeries<double> ras_flow_CDF = Inflow_DeNit_returnflow.BTC[0].GetCummulativeDistribution();
        ras_flow_CDF.writefile(Workingfolder + "Data/ras_flow_CDF.txt");

        CTimeSeries<double> ras_flow_inv_cum = ras_flow_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> ras_flow_PDF = Inflow_DeNit_returnflow.BTC[0].distribution(50,0);
        ras_flow_PDF.writefile(Workingfolder + "Data/ras_flow_PDF.txt");

        double ras_flow_mean = Inflow_DeNit_returnflow.BTC[0].Log().mean();
        double ras_flow_std = Inflow_DeNit_returnflow.BTC[0].Log().std();
        double ras_flow_autocorrelation_coeff = ras_flow_autocorrelation.AutoCorrelationCoeff();

        //MeoH (Methanol)
        CTimeSeries<double> MeOH_normal_score = Inflow_DeNit_MeOH.BTC[0].ConverttoNormalScore();
        MeOH_normal_score.writefile(Workingfolder + "Data/MeOH_normal_score.txt");

        CTimeSeries<double> MeOH_autocorrelation = MeOH_normal_score.AutoCorrelation(10,0.5);
        MeOH_autocorrelation.writefile(Workingfolder + "Data/MeOH_autocorrelation.txt");

        CTimeSeries<double> MeOH_CDF = Inflow_DeNit_MeOH.BTC[0].GetCummulativeDistribution();
        MeOH_CDF.writefile(Workingfolder + "Data/MeOH_CDF.txt");

        CTimeSeries<double> MeOH_inv_cum = MeOH_CDF.inverse_cumulative_uniform(100);

        CTimeSeries<double> MeOH_PDF = Inflow_DeNit_MeOH.BTC[0].distribution(50,0);
        MeOH_PDF.writefile(Workingfolder + "Data/MeOH_PDF.txt");

        double MeOH_mean = Inflow_DeNit_MeOH.BTC[0].Log().mean();
        double MeOH_std = Inflow_DeNit_MeOH.BTC[0].Log().std();
        double MeOH_autocorrelation_coeff = MeOH_autocorrelation.AutoCorrelationCoeff();

        //everything
        CTimeSeriesSet<double> normal_scores = Inflow_DeNit.ConverttoNormalScore();
        normal_scores.writetofile(Workingfolder + "Data/Normal_Scores.txt");

        CTimeSeriesSet<double> autocorrelations = normal_scores.AutoCorrelation(10,0.5);
        autocorrelations.writetofile(Workingfolder + "Data/autocorrelations.txt");

        CTimeSeriesSet<double> CDFs = Inflow_DeNit.GetCummulativeDistribution();
        CDFs.writetofile(Workingfolder + "Data/CDFs.txt");

        CTimeSeriesSet<double> PDFs = Inflow_DeNit.distribution(50,Inflow_DeNit.nvars, 0);
        PDFs.writetofile(Workingfolder + "Data/PDFs.txt");

        CMatrix correlation_matrix = normal_scores.correlation(0,normal_scores.nvars);
        correlation_matrix.writetofile(Workingfolder + "Data/correlation_matrix.txt");

        vector<double> means = Inflow_DeNit.mean(0);
        vector<double> stds = Inflow_DeNit.std(0);

        vector<double> logmeans = Inflow_DeNit.Log().mean(0);
        vector<double> logstds = Inflow_DeNit.Log().std(0);

        CTimeSeries<double> means2 = means;
        means2.writefile(Workingfolder + "Data/means.txt");
        CTimeSeries<double> stds2 = stds;
        stds2.writefile(Workingfolder + "Data/stds.txt");
        CTimeSeries<double> logmeans2 = logmeans;
        logmeans2.writefile(Workingfolder + "Data/logmeans.txt");
        CTimeSeries<double> logstds2 = logstds;
        logstds2.writefile(Workingfolder + "Data/logstds.txt");

    //OUProcess
    CTimeSeries<double> OUP_Inflow_Q_NS; // Discharge (m3/day)
    OUP_Inflow_Q_NS.CreateOUProcess(0,Simulation_time_Calc,dt,flow_autocorrelation_coeff);
    OUP_Inflow_Q_NS.writefile(Workingfolder + "OUP_Inflow_Q_NS_tvif.csv");
    //vector<double> Q_params; Q_params.push_back(flow_mean); Q_params.push_back(flow_std);
    //OUP_Inflow_Q = OUP_Inflow_Q_NS.MapfromNormalScoreToDistribution("lognormal", Q_params);
    OUP_Inflow_Q = OUP_Inflow_Q_NS.MapfromNormalScoreToDistribution(flow_inv_cum);
    OUP_Inflow_Q.writefile(Workingfolder + "OUP_Inflow_Q_tvif.csv");
    ToWrite.append(OUP_Inflow_Q,"Inflow");

    CTimeSeries<double> OUP_Flow_WAS_NS; // Wasteflow (WAS) (m3/day)
    OUP_Flow_WAS_NS.CreateOUProcess(0,Simulation_time_Calc,dt,was_flow_autocorrelation_coeff);
    OUP_Flow_WAS_NS.writefile(Workingfolder + "OUP_Flow_WAS_NS_tvf.csv");
    //vector<double> WAS_params; WAS_params.push_back(was_flow_mean); WAS_params.push_back(was_flow_std);
    //OUP_Flow_WAS = OUP_Flow_WAS_NS.MapfromNormalScoreToDistribution("lognormal", WAS_params);
    OUP_Flow_WAS = OUP_Flow_WAS_NS.MapfromNormalScoreToDistribution(was_flow_inv_cum);
    OUP_Flow_WAS.writefile(Workingfolder + "OUP_Flow_WAS_tvf.csv");
    ToWrite.append(OUP_Flow_WAS,"WAS_flow");

    CTimeSeries<double> OUP_Flow_RAS_NS; // Returnflow (RAS) (m3/day)
    OUP_Flow_RAS_NS.CreateOUProcess(0,Simulation_time_Calc,dt,ras_flow_autocorrelation_coeff);
    OUP_Flow_RAS_NS.writefile(Workingfolder + "OUP_Flow_RAS_NS_tvf.csv");
    //vector<double> RAS_params; RAS_params.push_back(ras_flow_mean); RAS_params.push_back(ras_flow_std);
    //OUP_Flow_RAS = OUP_Flow_RAS_NS.MapfromNormalScoreToDistribution("lognormal", RAS_params);
    OUP_Flow_RAS = OUP_Flow_RAS_NS.MapfromNormalScoreToDistribution(ras_flow_inv_cum);
    OUP_Flow_RAS.writefile(Workingfolder + "OUP_Flow_RAS_tvf.csv");
    ToWrite.append(OUP_Flow_RAS,"RAS_flow");


    OUP_Flow_r_r_st=OUP_Inflow_Q+OUP_Flow_RAS;
    OUP_Flow_st_sb=OUP_Flow_RAS+OUP_Flow_WAS;
    OUP_Flow_st_c=OUP_Inflow_Q-OUP_Flow_WAS;


    OUP_Flow_r_r_st.writefile(Workingfolder + "OUP_r_r_st_tvf.csv"); //
    OUP_Flow_st_sb.writefile(Workingfolder + "OUP_st_sb_tvf.csv"); //
    OUP_Flow_st_c.writefile(Workingfolder + "OUP_st_c_tvf.csv"); //

    //VSS (2)
    CTimeSeries<double> OUP_VSS_Flow_NS; // VSS (2)
    OUP_VSS_Flow_NS.CreateOUProcess(0,Simulation_time_Calc,dt,VSS_autocorrelation_coeff);
    OUP_VSS_Flow_NS.writefile(Workingfolder + "OUP_VSS_Flow_NS_tvf.csv");
    //vector<double> VSS_params; VSS_params.push_back(VSS_mean); VSS_params.push_back(VSS_std);
    //OUP_VSS_Flow = OUP_VSS_Flow_NS.MapfromNormalScoreToDistribution("lognormal", VSS_params);
    OUP_VSS_Flow = OUP_VSS_Flow_NS.MapfromNormalScoreToDistribution(VSS_inv_cum);
    OUP_VSS_Flow.writefile(Workingfolder + "OUP_VSS_Flow_tvf.csv");
    ToWrite.append(OUP_VSS_Flow,"VSS_inflow");


    //sCOD (8)
    CTimeSeries<double> OUP_sCOD_Flow_NS; // sCOD (8)
    OUP_sCOD_Flow_NS.CreateOUProcess(0,Simulation_time_Calc,dt,sCOD_autocorrelation_coeff);
    OUP_sCOD_Flow_NS.writefile(Workingfolder + "OUP_sCOD_Flow_NS_tvf.csv");
    //vector<double> sCOD_params; sCOD_params.push_back(sCOD_mean); sCOD_params.push_back(sCOD_std);
    //OUP_sCOD_Flow = OUP_sCOD_Flow_NS.MapfromNormalScoreToDistribution("lognormal", sCOD_params);
    OUP_sCOD_Flow = OUP_sCOD_Flow_NS.MapfromNormalScoreToDistribution(sCOD_inv_cum);
    OUP_sCOD_Flow.writefile(Workingfolder + "OUP_sCOD_Flow_tvf.csv");
    ToWrite.append(OUP_sCOD_Flow,"sCOD_inflow");

    //NH3 (9)
    CTimeSeries<double> OUP_NH3_Flow_NS; // NH3 (9)
    OUP_NH3_Flow_NS.CreateOUProcess(0,Simulation_time_Calc,dt,NH3_autocorrelation_coeff);
    OUP_NH3_Flow_NS.writefile(Workingfolder + "OUP_NH3_Flow_NS_tvf.csv");
    //vector<double> NH3_params; NH3_params.push_back(NH3_mean); NH3_params.push_back(NH3_std);
    //OUP_NH3_Flow = OUP_NH3_Flow_NS.MapfromNormalScoreToDistribution("lognormal", NH3_params);
    OUP_NH3_Flow = OUP_NH3_Flow_NS.MapfromNormalScoreToDistribution(NH3_inv_cum);
    OUP_NH3_Flow.writefile(Workingfolder + "OUP_NH3_Flow_tvf.csv");
    ToWrite.append(OUP_NH3_Flow,"NH3_inflow");

    //NO3 (10)
    CTimeSeries<double> OUP_NO3_Flow_NS; // NO3 (10)
    OUP_NO3_Flow_NS.CreateOUProcess(0,Simulation_time_Calc,dt,NO3_autocorrelation_coeff);
    OUP_NO3_Flow_NS.writefile(Workingfolder + "OUP_NO3_Flow_NS_tvf.csv");
    //vector<double> NO3_params; NO3_params.push_back(NO3_mean); NO3_params.push_back(NO3_std);
    //OUP_NO3_Flow = OUP_NO3_Flow_NS.MapfromNormalScoreToDistribution("lognormal", NO3_params);
    OUP_NO3_Flow = OUP_NO3_Flow_NS.MapfromNormalScoreToDistribution(NO3_inv_cum);
    OUP_NO3_Flow.writefile(Workingfolder + "OUP_NO3_Flow_tvf.csv");
    ToWrite.append(OUP_NO3_Flow,"NO3_inflow");

    //MeOH (Methanol)
    CTimeSeries<double> OUP_MeOH_Flow_NS; // MeOH (10)
    OUP_MeOH_Flow_NS.CreateOUProcess(0,Simulation_time_Calc,dt,MeOH_autocorrelation_coeff);
    OUP_MeOH_Flow_NS.writefile(Workingfolder + "OUP_MeOH_Flow_NS_tvf.csv");
    //vector<double> MeOH_params; MeOH_params.push_back(MeOH_mean); MeOH_params.push_back(MeOH_std);
    //OUP_MeOH_Flow = OUP_MeOH_Flow_NS.MapfromNormalScoreToDistribution("lognormal", MeOH_params);
    OUP_MeOH_Flow = OUP_MeOH_Flow_NS.MapfromNormalScoreToDistribution(MeOH_inv_cum);
    OUP_MeOH_Flow.writefile(Workingfolder + "OUP_MeoH_Flow_tvf.csv");
    ToWrite.append(OUP_MeOH_Flow,"MeoH_inflow");


        OUP_Inflow_S_i=c_S_i*OUP_sCOD_Flow; // 8

        OUP_Inflow_S_S=c_S_S*OUP_sCOD_Flow; //8
        OUP_Inflow_X_S=c_X_S*OUP_VSS_Flow; //2
        OUP_Inflow_X_p=c_X_p*OUP_VSS_Flow; //2
        OUP_Inflow_S_NO=c_S_NO*OUP_NO3_Flow; //10
        OUP_Inflow_S_NH=c_S_NH*OUP_NH3_Flow; //9
        OUP_Inflow_S_ND=c_S_ND*OUP_sCOD_Flow; //8
        OUP_Inflow_X_ND=c_X_ND*OUP_VSS_Flow; //2

        OUP_Inflow_MeOH=OUP_MeOH_Flow; //MeOH

        //Writing to file
        OUP_Inflow_S_i.writefile(Workingfolder + "OUP_Inflow_S_i.csv");
        OUP_Inflow_S_S.writefile(Workingfolder + "OUP_Inflow_S_S.csv");
        OUP_Inflow_X_S.writefile(Workingfolder + "OUP_Inflow_X_S.csv");
        OUP_Inflow_X_p.writefile(Workingfolder + "OUP_Inflow_X_p.csv");
        OUP_Inflow_S_NO.writefile(Workingfolder + "OUP_Inflow_S_NO.csv");
        OUP_Inflow_S_NH.writefile(Workingfolder + "OUP_Inflow_S_NH.csv");
        OUP_Inflow_S_ND.writefile(Workingfolder + "OUP_Inflow_S_ND.csv");
        OUP_Inflow_X_ND.writefile(Workingfolder + "OUP_Inflow_X_ND.csv");

        OUP_Inflow_MeOH.writefile(Workingfolder + "OUP_Inflow_MeOH.csv");


    // Assigninig Time Variable inflows by OUProcess

    system->block("Reactor_Flex(1)")->SetProperty("time_variable_inflow",Workingfolder + "OUP_Inflow_Q_tvif.csv"); // Discharge (m3/day)

    system->block("Reactor_Flex(1)")->SetProperty("S_i:time_variable_inflow_concentration",Workingfolder + "OUP_Inflow_S_i.csv");
    system->block("Reactor_Flex(1)")->SetProperty("S_S:time_variable_inflow_concentration",Workingfolder + "OUP_Inflow_S_S.csv");
    system->block("Reactor_Flex(1)")->SetProperty("X_S:time_variable_inflow_concentration",Workingfolder + "OUP_Inflow_X_S.csv");
    system->block("Reactor_Flex(1)")->SetProperty("X_p:time_variable_inflow_concentration",Workingfolder + "OUP_Inflow_X_p.csv");
    system->block("Reactor_Flex(1)")->SetProperty("S_NO:time_variable_inflow_concentration",Workingfolder + "OUP_Inflow_S_NO.csv");
    system->block("Reactor_Flex(1)")->SetProperty("S_NH:time_variable_inflow_concentration",Workingfolder + "OUP_Inflow_S_NH.csv");
    system->block("Reactor_Flex(1)")->SetProperty("S_ND:time_variable_inflow_concentration",Workingfolder + "OUP_Inflow_S_ND.csv");
    system->block("Reactor_Flex(1)")->SetProperty("X_ND:time_variable_inflow_concentration",Workingfolder + "OUP_Inflow_X_ND.csv");

    system->block("Reactor_Flex(5)")->SetProperty("S_M:external_mass_flow_timeseries",Workingfolder + "OUP_Inflow_MeOH.csv");


    // Flex_flow Links for Reactor_Flex
    for (int i=0; i<n_tanks-1; i++)
    {   Link l_r_st;
        l_r_st.SetQuantities(system, "Flex_flow");
        l_r_st.SetName("Reactor_Flex(" + aquiutils::numbertostring(i+1)+"_" + aquiutils::numbertostring(i+2) + ")");
        l_r_st.SetType("Flex_flow");
        //l_r_st.SetProperty("flow", Workingfolder + "OUP_r_r_st_tvf.csv");
        l_r_st.SetVal("flow_factor",v_flow_factor_i);
        system->AddLink(l_r_st, "Reactor_Flex(" + aquiutils::numbertostring(i+1)+")", "Reactor_Flex(" + aquiutils::numbertostring(i+2)+")", false);
    }


    // Flex_flow Links
    Link l_r_st;
    l_r_st.SetQuantities(system, "Flex_flow");
    l_r_st.SetName("Reactor_Flex(" + aquiutils::numbertostring(n_tanks) + ") - Settling element top");
    l_r_st.SetType("Flex_flow");
    //l_r_st.SetProperty("flow", Workingfolder + "OUP_r_r_st_tvf.csv");
    l_r_st.SetVal("flow_factor",v_flow_factor_i);
    system->AddLink(l_r_st, "Reactor_Flex(" + aquiutils::numbertostring(n_tanks) + ")", "Settling element top", false);

    Link l_st_c;
    l_st_c.SetQuantities(system, "Flex_flow");
    l_st_c.SetName("Settling element top - Clarifier");
    l_st_c.SetType("Flex_flow");
    //l_st_c.SetProperty("flow", Workingfolder + "OUP_st_c_tvf.csv");
    l_st_c.SetVal("flow_factor",v_flow_factor_o);
    system->AddLink(l_st_c, "Settling element top", "Clarifier", false);

    Link l_sb_was;
    l_sb_was.SetQuantities(system, "Flex_flow");
    l_sb_was.SetName("Settling element bottom - WAS");
    l_sb_was.SetType("Flex_flow");
    //l_sb_was.SetProperty("flow", Workingfolder + "OUP_Flow_WAS_NS_tvf.csv");
    l_sb_was.SetVal("flow_factor",v_flow_factor_o);
    system->AddLink(l_sb_was, "Settling element bottom", "WAS", false);

    // Time-Dependent interface Link
    Link l_st_sb;
    l_st_sb.SetQuantities(system, "Settling element interface (time variable)");
    l_st_sb.SetName("Settling element top - Settling element bottom");
    l_st_sb.SetType("Settling element interface (time variable)");
    l_st_sb.SetProperty("flow", Workingfolder + "OUP_st_sb_tvf.csv");
    l_st_sb.SetVal("area", v_st_sb_area);
    system->AddLink(l_st_sb, "Settling element top", "Settling element bottom", false);

    // Time-Dependent flow Links
    Link l_sb_r;
    l_sb_r.SetQuantities(system, "Time-Dependent flow");
    l_sb_r.SetName("Settling element bottom - Reactor_Flex(1)");
    l_sb_r.SetType("Time-Dependent flow");
    l_sb_r.SetProperty("flow", Workingfolder + "OUP_Flow_RAS_tvf.csv");
    system->AddLink(l_sb_r, "Settling element bottom", "Reactor_Flex(1)", false);

}

    // Observations


    // For Calibration
if (Calibration)
{
    Observation R5_S_NO_inflow_cn; // T3B

    R5_S_NO_inflow_cn.SetQuantities(system, "Observation");
    R5_S_NO_inflow_cn.SetProperty("expression","S_NO:concentration"); // time_variable_inflow_concentration
    R5_S_NO_inflow_cn.SetProperty("object","Reactor_Flex(5)");
    R5_S_NO_inflow_cn.SetName("R5_S_NO_Concentration");
    R5_S_NO_inflow_cn.SetType("Observation");
    system->AddObservation(R5_S_NO_inflow_cn,false);

    Observation R6_S_NO_inflow_cn; // T4

    R6_S_NO_inflow_cn.SetQuantities(system, "Observation");
    R6_S_NO_inflow_cn.SetProperty("expression","S_NO:concentration");
    R6_S_NO_inflow_cn.SetProperty("object","Reactor_Flex(6)");
    R6_S_NO_inflow_cn.SetName("R6_S_NO_Concentration");
    R6_S_NO_inflow_cn.SetType("Observation");
    system->AddObservation(R6_S_NO_inflow_cn,false);

    Observation R8_S_NO_inflow_cn; // T5B (Effluent)

    R8_S_NO_inflow_cn.SetQuantities(system, "Observation");
    R8_S_NO_inflow_cn.SetProperty("expression","S_NO:concentration");
    R8_S_NO_inflow_cn.SetProperty("object","Reactor_Flex(8)");
    R8_S_NO_inflow_cn.SetName("R8_S_NO_Concentration");
    R8_S_NO_inflow_cn.SetType("Observation");
    system->AddObservation(R8_S_NO_inflow_cn,false);

    Observation R4_sCOD_inflow_cn; // T3A

    R4_sCOD_inflow_cn.SetQuantities(system, "Observation");
    R4_sCOD_inflow_cn.SetProperty("expression","S_i:concentration+S_S:concentration+S_M:concentration");
    R4_sCOD_inflow_cn.SetProperty("object","Reactor_Flex(4)");
    R4_sCOD_inflow_cn.SetName("R4_sCOD_Concentration");
    R4_sCOD_inflow_cn.SetType("Observation");
    system->AddObservation(R4_sCOD_inflow_cn,false);

    Observation R5_sCOD_inflow_cn; // T3B

    R5_sCOD_inflow_cn.SetQuantities(system, "Observation");
    R5_sCOD_inflow_cn.SetProperty("expression","S_i:concentration+S_S:concentration+S_M:concentration");
    R5_sCOD_inflow_cn.SetProperty("object","Reactor_Flex(5)");
    R5_sCOD_inflow_cn.SetName("R5_sCOD_Concentration");
    R5_sCOD_inflow_cn.SetType("Observation");
    system->AddObservation(R5_sCOD_inflow_cn,false);

    Observation R6_sCOD_inflow_cn; // T4

    R6_sCOD_inflow_cn.SetQuantities(system, "Observation");
    R6_sCOD_inflow_cn.SetProperty("expression","S_i:concentration+S_S:concentration+S_M:concentration");
    R6_sCOD_inflow_cn.SetProperty("object","Reactor_Flex(6)");
    R6_sCOD_inflow_cn.SetName("R6_sCOD_Concentration");
    R6_sCOD_inflow_cn.SetType("Observation");
    system->AddObservation(R6_sCOD_inflow_cn,false);

    Observation R6_TKN_inflow_cn; // T4

    R6_TKN_inflow_cn.SetQuantities(system, "Observation");
    R6_TKN_inflow_cn.SetProperty("expression","0.086*X_BH:concentration+0.086*X_BM:concentration+0.086*X_BA:concentration+0.086*X_p:concentration+S_NH:concentration+S_ND:concentration+X_ND:concentration");
    R6_TKN_inflow_cn.SetProperty("object","Reactor_Flex(6)");
    R6_TKN_inflow_cn.SetName("R6_TKN_Concentration");
    R6_TKN_inflow_cn.SetType("Observation");
    system->AddObservation(R6_TKN_inflow_cn,false);

    Observation R6_VSS_inflow_cn; // T4

    R6_VSS_inflow_cn.SetQuantities(system, "Observation");
    R6_VSS_inflow_cn.SetProperty("expression","0.556*X_S:concentration+0.704*X_BH:concentration+0.704*X_BM:concentration+0.704*X_BA:concentration+0.704*X_p:concentration");
    R6_VSS_inflow_cn.SetProperty("object","Reactor_Flex(6)");
    R6_VSS_inflow_cn.SetName("R6_VSS_Concentration");
    R6_VSS_inflow_cn.SetType("Observation");
    system->AddObservation(R6_VSS_inflow_cn,false);

    Observation R8_VSS_inflow_cn; // T5B

    R8_VSS_inflow_cn.SetQuantities(system, "Observation");
    R8_VSS_inflow_cn.SetProperty("expression","0.556*X_S:concentration+0.704*X_BH:concentration+0.704*X_BM:concentration+0.704*X_BA:concentration+0.704*X_p:concentration");
    R8_VSS_inflow_cn.SetProperty("object","Reactor_Flex(8)");
    R8_VSS_inflow_cn.SetName("R8_VSS_Concentration");
    R8_VSS_inflow_cn.SetType("Observation");
    system->AddObservation(R8_VSS_inflow_cn,false);

    Observation R5_MeOH_inflow_cn; // T3B

    R5_MeOH_inflow_cn.SetQuantities(system, "Observation");
    R5_MeOH_inflow_cn.SetProperty("expression","S_M:concentration");
    R5_MeOH_inflow_cn.SetProperty("object","Reactor_Flex(5)");
    R5_MeOH_inflow_cn.SetName("R5_MeOH_Concentration");
    R5_MeOH_inflow_cn.SetType("Observation");
    system->AddObservation(R5_MeOH_inflow_cn,false);

    Observation R6_MeOH_inflow_cn; // T4

    R6_MeOH_inflow_cn.SetQuantities(system, "Observation");
    R6_MeOH_inflow_cn.SetProperty("expression","S_M:concentration");
    R6_MeOH_inflow_cn.SetProperty("object","Reactor_Flex(6)");
    R6_MeOH_inflow_cn.SetName("R6_MeOH_Concentration");
    R6_MeOH_inflow_cn.SetType("Observation");
    system->AddObservation(R6_MeOH_inflow_cn,false);

    Observation R8_S_NH_inflow_cn; // T5B

    R8_S_NH_inflow_cn.SetQuantities(system, "Observation");
    R8_S_NH_inflow_cn.SetProperty("expression","S_NH:concentration");
    R8_S_NH_inflow_cn.SetProperty("object","Reactor_Flex(8)");
    R8_S_NH_inflow_cn.SetName("R8_S_NH_Concentration");
    R8_S_NH_inflow_cn.SetType("Observation");
    system->AddObservation(R8_S_NH_inflow_cn,false);
}

    // For FFNWrapper (OUP)
if (OUP)
{
    Observation total_inflow;

    total_inflow.SetQuantities(system, "Observation");
    total_inflow.SetProperty("expression","inflow");
    total_inflow.SetProperty("object","Reactor_Flex(1)");
    total_inflow.SetName("Inflow");
    total_inflow.SetType("Observation");
    system->AddObservation(total_inflow,false);

    Observation WAS_flow;

    WAS_flow.SetQuantities(system, "Observation");
    WAS_flow.SetProperty("expression","Dummy_timeseries");
    WAS_flow.SetProperty("object","WAS");
    WAS_flow.SetName("WAS_Flow");
    WAS_flow.SetType("Observation");
    system->AddObservation(WAS_flow,false);

    /*
    Observation WAS_flow;

    WAS_flow.SetQuantities(system, "Observation");
    WAS_flow.SetProperty("expression","flow");
    WAS_flow.SetProperty("object","Settling element bottom - WAS");
    WAS_flow.SetName("WAS_Flow");
    WAS_flow.SetType("Observation");
    system->AddObservation(WAS_flow,false);
    */

    Observation RAS_flow;

    RAS_flow.SetQuantities(system, "Observation");
    RAS_flow.SetProperty("expression","flow");
    RAS_flow.SetProperty("object","Settling element bottom - Reactor_Flex(1)");
    RAS_flow.SetName("RAS_Flow");
    RAS_flow.SetType("Observation");
    system->AddObservation(RAS_flow,false);

    Observation Temp;

    Temp.SetQuantities(system, "Observation");
    Temp.SetProperty("expression","Dummy_timeseries");
    Temp.SetProperty("object","Clarifier");
    Temp.SetName("Temp");
    Temp.SetType("Observation");
    system->AddObservation(Temp,false);


    Observation R1_VSS_inflow_cn;

    R1_VSS_inflow_cn.SetQuantities(system, "Observation");
    //R1_VSS_inflow_cn.SetProperty("expression","0.556*X_S:time_variable_inflow_concentration+0.704*X_BH:time_variable_inflow_concentration+0.704*X_BM:time_variable_inflow_concentration+0.704*X_BA:time_variable_inflow_concentration+0.704*X_p:time_variable_inflow_concentration");
    R1_VSS_inflow_cn.SetProperty("expression","X_S:time_variable_inflow_concentration/" + QString::number(c_X_S).toStdString());
    R1_VSS_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    R1_VSS_inflow_cn.SetName("R1_VSS_Inflow_Concentration");
    R1_VSS_inflow_cn.SetType("Observation");
    system->AddObservation(R1_VSS_inflow_cn,false);

    Observation R1_sCOD_inflow_cn;

    R1_sCOD_inflow_cn.SetQuantities(system, "Observation");
    //R1_sCOD_inflow_cn.SetProperty("expression","S_i:time_variable_inflow_concentration+S_S:time_variable_inflow_concentration+S_M:time_variable_inflow_concentration");
    R1_sCOD_inflow_cn.SetProperty("expression","S_S:time_variable_inflow_concentration/" + QString::number(c_S_S).toStdString());
    R1_sCOD_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    R1_sCOD_inflow_cn.SetName("R1_sCOD_Inflow_Concentration");
    R1_sCOD_inflow_cn.SetType("Observation");
    system->AddObservation(R1_sCOD_inflow_cn,false);

    Observation R1_NH3_inflow_cn;

    R1_NH3_inflow_cn.SetQuantities(system, "Observation");
    R1_NH3_inflow_cn.SetProperty("expression","S_NH:time_variable_inflow_concentration");
    R1_NH3_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    R1_NH3_inflow_cn.SetName("R1_NH3_Inflow_Concentration");
    R1_NH3_inflow_cn.SetType("Observation");
    system->AddObservation(R1_NH3_inflow_cn,false);

    Observation R1_NO3_inflow_cn;

    R1_NO3_inflow_cn.SetQuantities(system, "Observation");
    R1_NO3_inflow_cn.SetProperty("expression","S_NO:time_variable_inflow_concentration");
    R1_NO3_inflow_cn.SetProperty("object","Reactor_Flex(1)");
    R1_NO3_inflow_cn.SetName("R1_NO3_Inflow_Concentration");
    R1_NO3_inflow_cn.SetType("Observation");
    system->AddObservation(R1_NO3_inflow_cn,false);

    Observation R1_MeOH_inflow_cn;

    R1_MeOH_inflow_cn.SetQuantities(system, "Observation");
    R1_MeOH_inflow_cn.SetProperty("expression","S_M:external_mass_flow_timeseries");
    R1_MeOH_inflow_cn.SetProperty("object","Reactor_Flex(5)");
    R1_MeOH_inflow_cn.SetName("R5_MeOH_Inflow");
    R1_MeOH_inflow_cn.SetType("Observation");
    system->AddObservation(R1_MeOH_inflow_cn,false);

    Observation Stl_t_S_NO_flow_cn;

    Stl_t_S_NO_flow_cn.SetQuantities(system, "Observation");
    Stl_t_S_NO_flow_cn.SetProperty("expression","S_NO:concentration");
    Stl_t_S_NO_flow_cn.SetProperty("object","Settling element top");
    Stl_t_S_NO_flow_cn.SetName("Stl_t_S_NO_Concentration");
    Stl_t_S_NO_flow_cn.SetType("Observation");
    system->AddObservation(Stl_t_S_NO_flow_cn,false);

    Observation Stl_t_S_NH3_flow_cn;

    Stl_t_S_NH3_flow_cn.SetQuantities(system, "Observation");
    Stl_t_S_NH3_flow_cn.SetProperty("expression","S_NH:concentration"); // use NH instead of NH3
    Stl_t_S_NH3_flow_cn.SetProperty("object","Settling element top");
    Stl_t_S_NH3_flow_cn.SetName("Stl_t_S_NH3_Concentration");
    Stl_t_S_NH3_flow_cn.SetType("Observation");
    system->AddObservation(Stl_t_S_NH3_flow_cn,false);

    Observation Stl_t_S_ND_flow_cn;

    Stl_t_S_ND_flow_cn.SetQuantities(system, "Observation");
    Stl_t_S_ND_flow_cn.SetProperty("expression","S_ND:concentration");
    Stl_t_S_ND_flow_cn.SetProperty("object","Settling element top");
    Stl_t_S_ND_flow_cn.SetName("Stl_t_S_ND_Concentration");
    Stl_t_S_ND_flow_cn.SetType("Observation");
    system->AddObservation(Stl_t_S_ND_flow_cn,false);

    Observation Stl_t_sCOD_flow_cn;

    Stl_t_sCOD_flow_cn.SetQuantities(system, "Observation");
    Stl_t_sCOD_flow_cn.SetProperty("expression","S_i:concentration+S_S:concentration+S_M:concentration");
    Stl_t_sCOD_flow_cn.SetProperty("object","Settling element top");
    Stl_t_sCOD_flow_cn.SetName("Stl_t_sCOD_Concentration");
    Stl_t_sCOD_flow_cn.SetType("Observation");
    system->AddObservation(Stl_t_sCOD_flow_cn,false);

    Observation Stl_t_VSS_flow_cn;

    Stl_t_VSS_flow_cn.SetQuantities(system, "Observation");
    Stl_t_VSS_flow_cn.SetProperty("expression","0.556*X_S:concentration+0.704*X_BH:concentration+0.704*X_BM:concentration+0.704*X_BA:concentration+0.704*X_p:concentration");
    Stl_t_VSS_flow_cn.SetProperty("object","Settling element top");
    Stl_t_VSS_flow_cn.SetName("Stl_t_VSS_Concentration");
    Stl_t_VSS_flow_cn.SetType("Observation");
    system->AddObservation(Stl_t_VSS_flow_cn,false);

    Observation Stl_t_TKN_flow_cn;

    Stl_t_TKN_flow_cn.SetQuantities(system, "Observation");
    Stl_t_TKN_flow_cn.SetProperty("expression","0.086*X_BH:concentration+0.086*X_BM:concentration+0.086*X_BA:concentration+0.06*X_p:concentration+S_NH:concentration+S_ND:concentration+X_ND:concentration"); // use NH instead of NH3
    Stl_t_TKN_flow_cn.SetProperty("object","Settling element top");
    Stl_t_TKN_flow_cn.SetName("Stl_t_TKN_Concentration");
    Stl_t_TKN_flow_cn.SetType("Observation");
    system->AddObservation(Stl_t_TKN_flow_cn,false);
}

    // Setting simulation time
    if (St)
    {
        system->SetSettingsParameter("simulation_end_time",Simulation_time);
    }

    else if (!St)
    {
        system->SetSettingsParameter("simulation_start_time",Simulation_start_time);
        system->SetSettingsParameter("simulation_end_time",Simulation_end_time);
    }

        system->SetSettingsParameter("maximum_time_allowed",86400*3);
        system->SetSettingsParameter("maximum_number_of_matrix_inverstions",1e6);

    // Observation writing time step
    system->SetSettingsParameter("initial_time_step",Initial_time_step);


    system->SetSystemSettings();

    cout<<"Populate functions"<<endl;
    system->PopulateOperatorsFunctions();
    cout<<"Variable parents"<<endl;
    system->SetVariableParents();

    /*
    CTimeSeries<double> NOx_out = system->GetObservedOutputs()["Stl_t_S_NO_Concentration"];
    NOx_out.make_uniform(dt);
    ToWrite.append(NOx_out,"Nox_Effluent");
    ToWrite.writetofile("ToEmulate.csv");
    */


    return true;
}


/*
CTimeSeries<double> Temp_normal_score = DeNit_Temp.BTC[0].ConverttoNormalScore();

Temp_normal_score.writefile(Workingfolder + "Data/Temp_normal_score.txt");

CTimeSeries<double> Temp_autocorrelation = Temp_normal_score.AutoCorrelation(10,0.5);
Temp_autocorrelation.writefile(Workingfolder + "Data/Temp_autocorrelation.txt");

CTimeSeries<double> Temp_CDF = DeNit_Temp.BTC[0].GetCummulativeDistribution();
Temp_CDF.writefile(Workingfolder + "Data/Temp_CDF.txt");

CTimeSeries<double> Temp_inv_cum = Temp_CDF.inverse_cumulative_uniform(100);

CTimeSeries<double> Temp_PDF = DeNit_Temp.BTC[0].distribution(50,0);
Temp_PDF.writefile(Workingfolder + "Data/Temp_PDF.txt");

double Temp_mean = DeNit_Temp.BTC[0].Log().mean();
double Temp_std = DeNit_Temp.BTC[0].Log().std();
double Temp_autocorrelation_coeff = Temp_autocorrelation.AutoCorrelationCoeff();

// Temperature
CTimeSeries<double> OUP_Temp_NS;
OUP_Temp_NS.CreateOUProcess(0,Simulation_time_Calc,dt,Temp_autocorrelation_coeff); // This is correct
*/
