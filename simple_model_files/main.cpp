#include "System.h"
#include "Script.h"
#include "qfileinfo.h"
#include "modelcreator.h"
#include "resultgrid.h"
#include "vtk.h"


int main(int argc, char *argv[])
{

    string Workingfolder="/home/behzad/Projects/ASM_Models/";
    ModelCreator ModCreate;
    for (int i=0; i<1; i++)
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

        cout<<"Getting results into grid"<<endl;
        ResultGrid resgrid(output,"theta",system);
        //cout<<"Writing VTPs"<<endl;
        //resgrid.WriteToVTP("Moisture_content",system->GetWorkingFolder()+"moisture.vtp");
        delete system;
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
