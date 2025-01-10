#include "System.h"
#include "Script.h"
#include "qfileinfo.h"
#include "modelcreator.h"
#include "modelcreator_flex.h"
#include "resultgrid.h"
//#include "vtk.h"


int main(int argc, char *argv[])
{

bool Flex = true; // Flex or normal reactor usage

double Realization = 1;

#ifdef Behzad
    string Workingfolder="/home/behzad/Projects/ASM_Models/";
    string Workingfolder_Flex="/home/behzad/Projects/ASM_Models/Flex/";
#else
    string Workingfolder="/home/arash/Projects/ASM_Models/";
    string Workingfolder_Flex="/home/arash/Projects/ASM_Models/Flex/";
#endif
#ifdef  Arash
    string Workingfolder = "/home/arash/Projects/ASM_Models/";
#endif //  Arash
#ifdef Arash_Windows
    string Workingfolder = "C:/Projects/ASM_Models/";
#endif // Arash_Windows

    /*
//Data analysis

    CTimeSeriesSet<double> Inflow_DeNit(Workingfolder + "Data/DeNit_Influent_Lump.txt",true);

//Flow
    CTimeSeries<double> flow_normal_score = Inflow_DeNit.BTC[0].ConverttoNormalScore();
    flow_normal_score.writefile(Workingfolder + "Data/flow_normal_score.txt");

    CTimeSeries<double> flow_autocorrelation = flow_normal_score.AutoCorrelation(10,0.5);
    flow_autocorrelation.writefile(Workingfolder + "Data/flow_autocorrelation.txt");

    CTimeSeries<double> flow_CDF = Inflow_DeNit.BTC[0].GetCummulativeDistribution();
    flow_CDF.writefile(Workingfolder + "Data/flow_CDF.txt");

    CTimeSeries<double> flow_PDF = Inflow_DeNit.BTC[0].distribution(50,0);
    flow_PDF.writefile(Workingfolder + "Data/flow_PDF.txt");

    double flow_mean = exp(Inflow_DeNit.BTC[0].Log().mean());
    double flow_std = Inflow_DeNit.BTC[0].Log().std();
    double flow_autocorrelation_coeff = flow_autocorrelation.AutoCorrelationCoeff();

//TSS
    CTimeSeries<double> TSS_normal_score = Inflow_DeNit.BTC[1].ConverttoNormalScore();
    TSS_normal_score.writefile(Workingfolder + "Data/TSS_normal_score.txt");

    CTimeSeries<double> TSS_autocorrelation = TSS_normal_score.AutoCorrelation(10,0.5);
    TSS_autocorrelation.writefile(Workingfolder + "Data/TSS_autocorrelation.txt");

    CTimeSeries<double> TSS_CDF = Inflow_DeNit.BTC[1].GetCummulativeDistribution();
    TSS_CDF.writefile(Workingfolder + "Data/TSS_CDF.txt");

    CTimeSeries<double> TSS_PDF = Inflow_DeNit.BTC[1].distribution(50,0);
    TSS_PDF.writefile(Workingfolder + "Data/TSS_PDF.txt");

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

*/


    ModelCreator_Flex ModCreate_Flex;

    if (Flex)
    {
    for (int i=0; i<Realization; i++)
    {
        System *system=new System();
        system->Clear();
        cout<<"Creating model "<< i << " ..." <<endl;
        ModCreate_Flex.Create_Flex(system);
        cout<<"Creating model done..." <<endl;

        system->SetWorkingFolder(Workingfolder_Flex);
        system->SetSilent(false);
        cout<<"Saving"<<endl;
        system->SavetoScriptFile(Workingfolder_Flex + "CreatedModel_Flex.ohq");

        cout<<"Solving ..."<<endl;
        system->Solve();

        cout<<"Writing outputs in '"<< system->GetWorkingFolder() + system->OutputFileName() +"'"<<endl;
        CTimeSeriesSet<double> output = system->GetOutputs();
        QString outputfilename = QString::fromStdString(system->OutputFileName()).split(".")[0] +"_" + QString::number(i) + ".txt";
        output.writetofile(system->GetWorkingFolder() + outputfilename.toStdString());

        cout<<"Writing outputs in '"<< system->GetWorkingFolder() + system->ObservedOutputFileName() +"'"<<endl;
        CTimeSeriesSet<double> selectedoutput = system->GetObservedOutputs();
        QString selectedoutputfilename = QString::fromStdString(system->ObservedOutputFileName()).split(".")[0] +"_" + QString::number(i) + ".txt";
        selectedoutput.writetofile(system->GetWorkingFolder() + selectedoutputfilename.toStdString());

        cout<<"Getting results into grid"<<endl;
        ResultGrid resgrid(output,"theta",system);
        //cout<<"Writing VTPs"<<endl;
        //resgrid.WriteToVTP("Moisture_content",system->GetWorkingFolder()+"moisture.vtp");
        delete system;
    }
    }

    ModelCreator ModCreate;

    if (!Flex)
    {
    for (int i=0; i<Realization; i++)
    {
        System *system=new System();
        system->Clear();
        cout<<"Creating model "<< i << " ..." <<endl;
        ModCreate.Create(system);
        cout<<"Creating model done..." <<endl;

        system->SetWorkingFolder(Workingfolder);
        system->SetSilent(false);
        cout<<"Saving"<<endl;
        system->SavetoScriptFile(Workingfolder + "CreatedModel.ohq");

        cout<<"Solving ..."<<endl;
        system->Solve();

        cout<<"Writing outputs in '"<< system->GetWorkingFolder() + system->OutputFileName() +"'"<<endl;
        CTimeSeriesSet<double> output = system->GetOutputs();
        QString outputfilename = QString::fromStdString(system->OutputFileName()).split(".")[0] +"_" + QString::number(i) + ".txt";
        output.writetofile(system->GetWorkingFolder() + outputfilename.toStdString());

        cout<<"Writing outputs in '"<< system->GetWorkingFolder() + system->ObservedOutputFileName() +"'"<<endl;
        CTimeSeriesSet<double> selectedoutput = system->GetObservedOutputs();
        QString selectedoutputfilename = QString::fromStdString(system->ObservedOutputFileName()).split(".")[0] +"_" + QString::number(i) + ".txt";
        selectedoutput.writetofile(system->GetWorkingFolder() + selectedoutputfilename.toStdString());

        //cout<<"Getting results into grid"<<endl;
        //ResultGrid resgrid(output,"theta",system);
        //cout<<"Writing VTPs"<<endl;
        //resgrid.WriteToVTP("Moisture_content",system->GetWorkingFolder()+"moisture.vtp");
        delete system;
    }
    }

    return 0;

}












/*CTimeSeries<double> OUEx1, OUEx2;
    OUEx1.CreateOUProcess(0,10,0.05,0.2);
    OUEx2.CreateOUProcess(0,10,0.05,5);
    vector<double> logparams; logparams.push_back(1); logparams.push_back(0.5);
    vector<double> expparams; expparams.push_back(3);
    CTimeSeries<double> OUEx1Log = OUEx1.MapfromNormalScoreToDistribution("lognormal", logparams);
    CTimeSeries<double> OUEx2Log = OUEx2.MapfromNormalScoreToDistribution("lognormal", logparams);

    CTimeSeries<double> OUEx1Exp = OUEx1.MapfromNormalScoreToDistribution("exp", expparams);
    CTimeSeries<double> OUEx2Exp = OUEx2.MapfromNormalScoreToDistribution("exp", expparams);

    OUEx1Log.writefile("/home/behzad/Projects/Settling_Models/OUEx1Log.txt");
    OUEx2Log.writefile("/home/behzad/Projects/Settling_Models/OUEx2Log.txt");

    OUEx1Exp.writefile("/home/behzad/Projects/Settling_Models/OUEx1Exp.txt");
    OUEx2Exp.writefile("/home/behzad/Projects/Settling_Models/OUEx2Exp.txt");
*/
