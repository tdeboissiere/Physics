
#include <string>
#include <iostream>
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <fstream>
#include <sstream>
#include <string>
#include "TMatrixD.h"
#include "TArrayD.h"
using namespace std;

int Nobs_true;
vector<double> pred_points;
vector<double> s_array;
vector<double> b_array;

double approx_factorial(int n)
{

    if (n<5) 
    {
        return TMath::Log(TMath::Factorial(n)); 
    } 
        
    else 
    {
        double fact = 0.5*TMath::Log(TMath::Pi()) + n*TMath::Log(n/TMath::Exp(1)) + (1./6.)*TMath::Log(8*pow(n,3)+4*pow(n,2)+n+1./30.);
        return fact;
    }
        
}

double likelihood2D(double * x, double * par)
{


    double s = x[0];
    double b = x[1];
    int Nsimu = s_array.size();
    double poisson_part = (s+b) - Nsimu*TMath::Log(s+b) + approx_factorial(Nsimu);
    double prod_part = 0;
    for (int k = 0; k<Nsimu; k++) 
    {
        prod_part += TMath::Log(s/(s+b)*s_array[k] + b/(s+b)*b_array[k]);
    }
    double constrain_part = (1./(sqrt(2*TMath::Pi()*Nobs_true))) * TMath::Exp(-(0.5/Nobs_true)*pow(b-Nobs_true,2));
    return poisson_part - prod_part - TMath::Log(constrain_part);

}

double profile_likelihood(double * x, double * par)
{


    double s = par[0];
    double b = x[0];
    int Nsimu = s_array.size();
    double poisson_part = (s+b) - Nsimu*TMath::Log(s+b) + approx_factorial(Nsimu);
    double prod_part = 0;
    for (int k = 0; k<Nsimu; k++) 
    {
        prod_part += TMath::Log(s/(s+b)*s_array[k] + b/(s+b)*b_array[k]);
    }
    double constrain_part = (1./(sqrt(2*TMath::Pi()*Nobs_true))) * TMath::Exp(-(0.5/Nobs_true)*pow(b-Nobs_true,2));
    return poisson_part - prod_part - TMath::Log(constrain_part);

}


int main() 
{

    //Initialize Random generation
    gRandom->SetSeed(0);
    TRandom3 u(0);

    // Pick an overall bckg normalisation
    int Nexp_bckg = 37468;

    // Pick number of signal events to generate
    int Nsignal = 0;
    
    //Number of simulations
    int simu_number = 50;

    // Read bckg and signal PDF
    TFile f("./ROOT_files/FID837_BDT_bckg_and_sig_4GeV.root", "read");
    TF1* f_bckg_norm = (TF1*)f.Get("bckg_norm");
    TF1* f_signal_norm = (TF1*)f.Get("signal_norm");

    //Output best fit value TH1F 
    TH1F * hbestfit = new TH1F("hbestfit", "hbestfit", 100, Nsignal-5*sqrt(Nsignal), Nsignal + 5*sqrt(Nsignal)); 
    TH1F * hstat = new TH1F("hstat", "hstat", 200, 0,5); 

    vector<double > test_stat_vec;

    for (int nsimu =0; nsimu<simu_number; nsimu++)
    {

        Nobs_true = u.Poisson(Nexp_bckg);

        pred_points.clear();
        s_array.clear();
        b_array.clear();
        // Create simulation of Nobs_true bckg events + desired number of signal events
        for (int k = 0; k<Nobs_true; k++)
        {
            pred_points.push_back(f_bckg_norm->GetRandom());
        }
        for (int k = 0; k<Nsignal; k++)
        {
            pred_points.push_back(f_signal_norm->GetRandom());
        }

        //Fill s_array and b_array 
        for (int k = 0; k<pred_points.size(); k++)
        {
            s_array.push_back(f_signal_norm->Eval(pred_points[k]));
            b_array.push_back(f_bckg_norm->Eval(pred_points[k]));
        }

        double min_x = 0;
        double min_y = 0;

        TF2 * like2D = new TF2("like2D", likelihood2D, 0.1,5000, 37000, 38000,0);
        like2D->GetMinimumXY(min_x, min_y);
        double like2D_optim = like2D->Eval(min_x, min_y);

        TF1 * profile_like = new TF1("profile_like", profile_likelihood, 37000, 38000,1);

        profile_like->SetParameter(0, Nsignal);
        profile_like->SetParameter(0, 125);
        double profile_optim = profile_like->GetMinimum();

        double test_stat = 0;
        if (min_x<125) {test_stat = 2*profile_optim -2*like2D_optim;}
        cout << "Best Nsignal,"<< min_x << ",Best Nbckg," << min_y << ",test stat," << test_stat << endl;
        if (min_x == min_x)
        {
            test_stat_vec.push_back(test_stat);
        }
        // // Draw as a check
        // TCanvas * cc = new TCanvas("cc", "cc");
        // profile_like->Draw();
        // cc->Print("profile.png");
        // delete cc;
        
        // hbestfit->Fill(min_x);
        //hstat->Fill(test_stat);

    }

    ofstream file("./Text_files/test_stat_obs_Nsignal125.txt");
    for (int k = 0; k<test_stat_vec.size(); k++)
    {
        file << test_stat_vec[k] << endl;
    }

    // ofstream file("./Text_files/test_stat_Nsignal200.txt");
    // for (int k = 0; k<test_stat_vec.size(); k++)
    // {
    //     file << test_stat_vec[k] << endl;
    // }
    // TCanvas * cc = new TCanvas("cc", "cc");
    // hstat->Draw();
    // cc->Print("hstat.png");

    return 0;

}

