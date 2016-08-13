
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TMath.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
using namespace std;


TF2* fWIMP2D;

double simulate_events(int Wmass) {

	double Ei, Er;


	//Fill histogram to get maximum (root can't find max for 2D functions, only min
	// Limit to 10 in ionisation to avoid crash. It doesn't influence the contour
	double xmin=2.44 ;double xmax=30.7; double ymin=0 ; double ymax=10.0;
	int binX=2000; int binY=1000;
	float binXf=2000.; float binYf=1000.;

	TH2F* hWIMP2D= new TH2F("h","h",binX, xmin, xmax, binY, ymin, ymax);
	for (int i = 1; i <= binX; i += 1)
	{
		for (int k = 1; k <=binY; k += 1) 
		{
		  hWIMP2D->SetBinContent(i,k,fWIMP2D->Eval(xmin+(xmax-xmin)*i/binXf, ymin+(ymax-ymin)*k/binYf));	
		}
	}

	int max_bin = hWIMP2D->GetMaximumBin();
	int max_bin_x, max_bin_y, max_bin_z;
	hWIMP2D->GetBinXYZ(max_bin, max_bin_x, max_bin_y, max_bin_z);
	double x_max, y_max;
	x_max = hWIMP2D->GetXaxis()->GetBinCenter(max_bin_x);
	y_max = hWIMP2D->GetYaxis()->GetBinCenter(max_bin_y);

	TH1F * h1 = new TH1F("h1", "h1", 10000, 0, 2*fWIMP2D->Eval(x_max, y_max));

	TH2F* hfull= new TH2F("hfull","hfull",binX, xmin, xmax, binY, ymin, ymax);
	TH2F* h95= new TH2F("h95","h95",binX, xmin, xmax, binY, ymin, ymax);



	// First get the histogram which we will integrate to find the 95 % contour
	for (int k = 0; k< 200000; k++) 
	{
		// fWIMP2D->GetRandom2(Ei,Er);
		hWIMP2D->GetRandom2(Ei, Er);
		double PDF_val = fWIMP2D->Eval(Ei, Er);
		h1->Fill(PDF_val);
	}

	// Find the value of the PDF (contour level) such that events with a PDF value between it and the max of the PDF
	// make up 95 % of the signal
	double PDF_cut_val = 0;
	for (int i = 1; i<10000; i++)
	{
		if (h1->Integral(i,10000)/h1->Integral()>0.95 && h1->Integral(i+1,10000)/h1->Integral()<0.95)
		{
			PDF_cut_val = h1->GetBinCenter(i) + 0.5*(h1->GetBinCenter(i+1)-h1->GetBinCenter(i));
			// cout << Wmass <<","<< h1->GetBinCenter(i) + 0.5*(h1->GetBinCenter(i+1)-h1->GetBinCenter(i)) << endl;
		}
	}

	//Run a check to see we get 95% of events indeed
	double counter = 0;
	for (int k = 0; k< 200000; k++) 
	{
		// fWIMP2D->GetRandom2(Ei,Er);
		hWIMP2D->GetRandom2(Ei, Er);
		double PDF_val = fWIMP2D->Eval(Ei, Er);
		if (PDF_val>PDF_cut_val)
		{
			h95->Fill(Ei, Er);
			counter+=1;
		}
		hfull->Fill(Ei,Er);
	}
	cout << Wmass <<","<< PDF_cut_val <<"," << counter/200000. << endl;

	// TFile fout("./ROOT_files/hist2D.root", "recreate");
	// hfull->Write();
	// h95->Write();
	// fout.Close();

	// TFile fout("hist.root", "recreate");
	// h1->Write();
	// fout.Close();
	
	delete hfull;
	delete h95;
	delete h1;
	delete hWIMP2D;

	return 0;

}





int main() {

	TString analysis = TString("3sveto_3sfid");

	//Get WIMP2D functions
	TFile file_WIMP(TString("./ROOT_files/WIMP_PDF2D_") + analysis + TString(".root"), "read");
	//Mass vector 
	vector<int> vmass;
	for (int k = 5; k<= 10; k++) {vmass.push_back(k);}
	// vmass.push_back(5);
	// vmass.push_back(7);
	// vmass.push_back(10); 
	vmass.push_back(15);
	vmass.push_back(20);
	vmass.push_back(25); 



	for(std::vector<int>::iterator it = vmass.begin(); it != vmass.end(); it++)
	{
		TString  s_Wmass= Form("%d", *it);
		fWIMP2D = (TF2*)file_WIMP.Get(TString("WIMP_") + s_Wmass + TString("_GeV"));
		simulate_events(*it);

	}

return 0;




}
