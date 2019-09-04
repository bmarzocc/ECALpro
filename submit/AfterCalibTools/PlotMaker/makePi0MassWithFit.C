#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h    
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++    

#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "./CMS_lumi.h"
//#include "/afs/cern.ch/work/m/mciprian/w_mass_analysis/CMSSW_8_0_25/src/WmassAnalysis/macros/utility.h"
#include "./utility.h"

using namespace std;
using namespace RooFit;

int Xtal_ID[170][360]={0};
int Xtal_Ieta[61200]={0};
int Xtal_Iphi[61200]={0};

int Xtal_Ix[14648]={0};
int Xtal_Iy[14648]={0};
int Xtal_Iz[14648]={0};

static string endcap_ix_iy_zside_ietaRing = "/afs/cern.ch/user/m/mciprian/public/ECALproTools/EE_xyzToEtaRing/eerings_modified.root";
static string deadXtalFileName = "/afs/cern.ch/user/m/mciprian/public/ECALproTools/test_DeadXtal_AlCaP0_Run2017B_3July_upToRun297723/h_DeadXtal.root";
static bool drawAllMassPlot = false;

//=================================================

class vectorManager {

public:
  vectorManager() { };

  vectorManager(const vector<TH1*> & histPtrs,
		const vector<string> & histNames,
		const vector<string> & histLegends
		)
  {
    histPtrs_    = vector<TH1*  >(histPtrs);
    histNames_   = vector<string>(histNames);
    histLegends_ = vector<string>(histLegends);
  };

  ~vectorManager() {};

  vector<TH1*>   getHistPtrs()    const { return histPtrs_;    };
  vector<string> getHistNames()   const { return histNames_;   };
  vector<string> getHistLegends() const { return histLegends_; };

  void addComponent(TH1* histPtr = NULL, const string& histName = "name", const string& histLegend = "leg")  { 
    histPtrs_.push_back(histPtr);
    histNames_.push_back(histName);
    histLegends_.push_back(histLegend);
  };

private:

  vector<TH1*> histPtrs_;
  vector<string> histNames_;
  vector<string> histLegends_;

};


//=================================================

bool noDeadXtalIn3x3matrixSeededByThisXtal(const TH2F* hDeadXtals = NULL, const int x = 1, const int y = 1) {

  // WARNING: it is assumed that the seed is already not adjacent to a gap, therefore we won't have abs(eta)=0 or abs(ieta)=85 or iphi=1 or iphi=360 for the seed 

  int nDeadXtals = 0;

  for (int xspan = x-1; xspan <= x+1 && nDeadXtals == 0; xspan++) {
    for (int yspan = y-1; yspan <= y+1 && nDeadXtals == 0; yspan++) {
      nDeadXtals += (int) (0.5 + hDeadXtals->GetBinContent(xspan,yspan)); // histogram returns float, to avoid bad truncation sum 0.5 and then round to integer 
    }
  }

  return (nDeadXtals == 0) ? true : false;

}

//=================================================

void drawHisto(TH1* hSum = NULL, 
	       const bool isEB = true, 
	       const string& outDir = "./", 
	       const string& hName = "", 
	       const double lumi = 8.6) {

  TGaxis::SetMaxDigits(3); 

  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetRightMargin(0.06);

  hSum->SetStats(0);
  hSum->SetLineColor(kBlack);
  hSum->SetMarkerColor(kBlack);
  hSum->SetMarkerStyle(20);
  hSum->SetMarkerSize(1);

  hSum->SetTitle(0);
  
  hSum->GetXaxis()->SetLabelSize(0.04);
  hSum->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  hSum->GetXaxis()->SetTitleSize(0.05);  
  hSum->GetXaxis()->SetTitleOffset(0.9);
  hSum->GetXaxis()->SetRangeUser(0.05,0.25);

  double maxY = hSum->GetBinContent(hSum->GetMaximumBin());
  hSum->GetYaxis()->SetRangeUser(0.0, 1.2*maxY);
  hSum->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
  hSum->GetYaxis()->SetTitleOffset(1.1);
  hSum->GetYaxis()->SetTitleSize(0.05);
  hSum->Draw("EP");

  canvas->RedrawAxis("sameaxis");
  if (lumi < 1.0) CMS_lumi(canvas,Form("%.2f",lumi),false,false);
  else CMS_lumi(canvas,Form("%.1f",lumi),false,false);
  setTDRStyle();

  string title = "pi0mass";
  title += isEB ? "_EB_" : "_EE_";
  canvas->SaveAs((outDir + title + hName + ".png").c_str());

}

//=================================================   

void doPi0MassWithFit(TH1* h,
		      const string& hName = "h",
		      const bool isEB = true,
		      const string& plotType = "etaring",
		      const string& outDir = "./",
		      const double lumi = 8.6,
		      const string& filename = "root://eoscms//store/group/dpg_ecal/alca_ecalcalib/bmarzocc/CMSPOS2019/AlCaP0_2018_ULrereco_1every2/iter_0/AlCaP0_2018_ULrereco_1every2_epsilonPlots.root"
		      ) 
{

  // force digits on axis (I wanted it only on y axis to trigger exponential notation, but it looks like it is not implemented for a single axis) 
  TGaxis::SetMaxDigits(3); 

  TH1::SetDefaultSumw2();
  gROOT->SetBatch(kTRUE);

  if (isEB) {

    for(int i = 0; i < 61200; i++)
      {
	int det_ID = EBDetId::detIdFromDenseIndex(i);

	EBDetId ebseed(det_ID);
        int ieta = ebseed.ieta();
        int iphi = ebseed.iphi();		
	Xtal_Ieta[i] = ieta;
	Xtal_Iphi[i] = iphi;
	//	cout<<i<<"   "<<ieta<<"  "<<iphi<<endl;
      }

  } else {

    // EE
    for(int i = 0; i < 14648; i++)
      {
	int det_ID = EEDetId::detIdFromDenseIndex(i);

	// TO BE TESTED
	EEDetId eeseed(det_ID);
        int ix = eeseed.ix();
        int iy = eeseed.iy();		       
	int iz = eeseed.zside();		
	Xtal_Ix[i] = ix;
	Xtal_Iy[i] = iy;
	Xtal_Iz[i] = iz;
	//	cout<<i<<"   "<<ix<<"  "<<iy<<endl;
      }

  }
  
  TH2F* hDeadXtalEB = NULL;
  TH2F* hDeadXtalEEm = NULL;
  TH2F* hDeadXtalEEp = NULL;

  TFile* deadXtalFile = TFile::Open(deadXtalFileName.c_str(),"READ");
  if (!deadXtalFile || !deadXtalFile->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<< deadXtalFileName <<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  hDeadXtalEB = (TH2F*) deadXtalFile->Get("rms_EB_r"); // #eta on x #phi on y
  hDeadXtalEEm = (TH2F*) deadXtalFile->Get("rms_EEm");
  hDeadXtalEEp = (TH2F*) deadXtalFile->Get("rms_EEp");

  if (!hDeadXtalEB || hDeadXtalEB == NULL) {
    cout << "Error: histogram rms_EB not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  if (!hDeadXtalEEm || hDeadXtalEEm == NULL) {
    cout << "Error: histogram rms_EEm not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  if (!hDeadXtalEEp || hDeadXtalEEp == NULL) {
    cout << "Error: histogram rms_EEp not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  // needed for EE to convert xyz to ietaRing
  TFile *EEetaRingFile = NULL;
  TH2F *hEE_etaRing = NULL;
  TH2F *hEEp_etaRing = NULL;
  TH2F *hEEm_etaRing = NULL;

  if (not isEB) {

    EEetaRingFile = new TFile((endcap_ix_iy_zside_ietaRing).c_str(),"READ");
    if (!EEetaRingFile || !EEetaRingFile->IsOpen()) {
      cout << "Error: file \"" << endcap_ix_iy_zside_ietaRing << "\" was not opened." << endl;
      exit(EXIT_FAILURE);
    }
    hEEp_etaRing = (TH2F*) EEetaRingFile->Get("hEEp");
    hEEm_etaRing = (TH2F*) EEetaRingFile->Get("hEEm");
    if (!hEEp_etaRing || hEEp_etaRing == NULL) {
      cout << "Error: histogram 'hEEp' not found in file ' " << endcap_ix_iy_zside_ietaRing << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    if (!hEEm_etaRing || hEEm_etaRing == NULL) {
      cout << "Error: histogram 'hEEm' not found in file ' " << endcap_ix_iy_zside_ietaRing << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }

  }

  TFile* f = TFile::Open(filename.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filename<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  f->cd();

  string directoryName = isEB ? "Barrel/" : "Endcap/";
  // bool IamInDirectory = f->cd(directoryName.c_str());
  // if (not IamInDirectory) {
  //   cout<<"Error: I couldn't change directory in file.\nApplication will be terminated."<<endl;
  //   exit(EXIT_FAILURE);
  // }

  TDirectory *dir = NULL;
  dir = f->GetDirectory(directoryName.c_str());
  if (!dir || dir == NULL) {
    cout<<"Error: I couldn't get directory in file.\nApplication will be terminated."<<endl;
    exit(EXIT_FAILURE);
  }

  dir->cd();

  UInt_t nObjectNotFound = 0;

  TH1F* hist = NULL;
  TH1F* hSum = NULL;

  bool isFirstHistogram = true;
  int nHistogramAdded = 0;
  int nTotalObjects = isEB ? 61200 : 14648;
  int nEvents = 0;

  if (plotType == "xtal") {
    
    string histName = isEB ? "epsilon_EB_iR_30003" : "epsilon_EE_iR_6397";
    hist = (TH1F*) dir->Get(histName.c_str());
    if (!hist) {
      cout << "Warning: TH1F object not found in file and plotType = " << plotType << ". Please check. Abort" <<endl;
      exit(EXIT_FAILURE);
    }
    hSum = new TH1F(*((TH1F*) hist->Clone("hSum")));

  } else {

    TIter next(dir->GetListOfKeys());
    TKey *key = NULL;

    while ( (key = (TKey*)next()) ) {

      cout.flush();
      if(nEvents % 50 == 0) cout<<"\r"<<"Crystals processed: "<<double(nEvents)/nTotalObjects*100<<" % ";
      //cout << "entry : " << nEvents << endl;                                                                                                                                
      nEvents++;

      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1F")) continue;

      hist = (TH1F*) key->ReadObj();
      if (!hist) {
	cout << "Warning: TH1F object not found in file. Skipping and going on with next object" <<endl;
	nObjectNotFound++;
	continue;
      }

      //cout << "check" << endl;

      string hname(hist->GetName());
      string hnameTag = isEB ? "epsilon_EB_iR_" : "epsilon_EE_iR_";
      if (hname.find(hnameTag.c_str()) == string::npos) continue;
    
      string fitIndexStr = ""; 
      fitIndexStr.assign(hname, hnameTag.size(), string::npos);
      //cout << "hname " << hname << "     fitIndexStr " << fitIndexStr << endl;
      Int_t fitIndex = std::stoi(fitIndexStr);

      // here we select some crystals with some algorithm
    
      bool conditionFulfilled = false;

      if (isEB) {

	int ieta = Xtal_Ieta[fitIndex];
	int iphi = Xtal_Iphi[fitIndex];

	// hardcoded, implement a flag to choose selection algorithm
	if (plotType == "etaring") {

	  // ieta == -2 and removing crystals near gaps in iphi and removing dead xtals or xtals adjacent to a dead xtals
	  conditionFulfilled = (ieta == -2 && iphi%20 != 0 && (iphi-1)%20 != 0 && noDeadXtalIn3x3matrixSeededByThisXtal(hDeadXtalEB,ieta+86,iphi));

	} else if (plotType == "xtal") {

	  if (ieta == -2 && iphi == 124) conditionFulfilled = true;

	} else {

	  conditionFulfilled = fabs(ieta) < 15 && iphi > 350 && iphi < 360;

	}

      } else {

	int ix = Xtal_Ix[fitIndex];
	int iy = Xtal_Iy[fitIndex];
	int iz = Xtal_Iz[fitIndex];


	if (plotType == "etaring") {
	
	  // etaring value is hardcoded (we used a crystal with ix,iy,iz = 19,83,-1 which is at etaRing=5, so we stick to that ring)
	  if (iz > 0) conditionFulfilled = hEEp_etaRing->GetBinContent(ix,iy) == 5 && noDeadXtalIn3x3matrixSeededByThisXtal(hDeadXtalEEp,ix,iy);
	  else        conditionFulfilled = hEEm_etaRing->GetBinContent(ix,iy) == 5 && noDeadXtalIn3x3matrixSeededByThisXtal(hDeadXtalEEm,ix,iy);

	} else if (plotType == "xtal") {

	  if (ix == 19 && iy == 83 && iz == -1) conditionFulfilled = true;

	} else {

	  int ixFromCenter = (ix < 51) ? (ix - 51) : (ix - 50); // avoid ixFromCenter = 0
	  int iyFromCenter = (iy < 51) ? (iy - 51) : (iy - 50); // avoid ixFromCenter = 0
	  double radius = sqrt(ixFromCenter*ixFromCenter + iyFromCenter*iyFromCenter);
	  //if (sqrt(radius) > 40) histToSum.push_back((TH1F*) hist->Clone());
	  if (radius > 42.0 && radius < 47.0) conditionFulfilled = true;

	}

      }

      if (conditionFulfilled) {
      
	if (isFirstHistogram) {
	  hSum = new TH1F(*((TH1F*) hist->Clone("hSum")));
	  isFirstHistogram = false;
	} else {
	  hSum->Add((TH1F*)hist->Clone());
	}
	nHistogramAdded++;
	//if (nHistogramAdded >= 50) break;
	//if (plotType == "xtal") break;
      
      }

    }

  }
  
  if (plotType != "xtal") {
    if (nHistogramAdded == 0) {
      cout << "Warning: no histogram used. End of programme." <<endl;
      exit(EXIT_FAILURE);
    }
    cout << "Selecting " << nHistogramAdded << " crystals in " << ((isEB) ? "EB" : "EE") << endl; 
  }
  if (nObjectNotFound > 0) cout << nObjectNotFound << " crystals were not found in file" << endl; 
  cout << "hSum->Integral() " << hSum->Integral() << endl;

  // it seems that the first time CMS_lumi is used the settings are screwed up 
  // produce a dummy plot (either do not save it or remove it) 
  //double lumi = 0.18; //in fb-1 
  TCanvas*ctmp = new TCanvas("ctmp","");
  ctmp->cd();
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  htmp1->Fill(0.5);
  htmp1->Draw("H");
  CMS_lumi(ctmp,Form("%.2f",lumi),false,false);
  setTDRStyle();
  delete htmp1;
  delete ctmp;

  if (drawAllMassPlot) drawHisto(hSum, isEB, outDir, hName, lumi);

  //h = new TH1F(hName.c_str(),"",hSum->GetNbinsX(),hSum->GetBinLowEdge(1),hSum->GetBinLowEdge(hSum->GetNbinsX()+1));
  h->SetBins(hSum->GetNbinsX(),hSum->GetBinLowEdge(1),hSum->GetBinLowEdge(hSum->GetNbinsX()+1));
  for (Int_t i = 1; i <= hSum->GetNbinsX(); i++) {
    h->SetBinContent(i, hSum->GetBinContent(i));
    h->SetBinError(i, hSum->GetBinError(i));
  }

  f->Close();
  delete f;
  deadXtalFile->Close();
  delete deadXtalFile;
  if (not isEB) {  
    EEetaRingFile->Close();
    delete EEetaRingFile;
  }

}

//===================================

void makePi0MassWithFit(const bool isEB = true,
			const bool singleXtalOnly = false,
			const string& outDir = "./",
			const double lumi = 58.83,
			const string& dirName = "AlCaP0_2018_ULrereco_1every2",
			const string& filePath = "root://eoscms//store/group/dpg_ecal/alca_ecalcalib/bmarzocc/CMSPOS2019/"
			) 
{

  createPlotDirAndCopyPhp(outDir);

  TH1F* h_xtal_iter0 = new TH1F();
  TH1F* h_xtal_iter2 = new TH1F();
  TH1F* h_xtal_iter4 = new TH1F();
  TH1F* h_etaring_iter0 = new TH1F();
  TH1F* h_etaring_iter2 = new TH1F();
  TH1F* h_etaring_iter4 = new TH1F();

  vectorManager* histData = new vectorManager();
  // (to be fixed) name format must be h_<type>_iter<n>, where type = "xtal" or "etaring", and n = 0,1,...
  histData->addComponent(h_xtal_iter0, "h_xtal_iter0", "1 xtal, 1st iter");
  histData->addComponent(h_xtal_iter2, "h_xtal_iter2", "1 xtal, 3rd iter");
  histData->addComponent(h_xtal_iter4, "h_xtal_iter4", "1 xtal, 5th iter");
  if (not singleXtalOnly) {
    histData->addComponent(h_etaring_iter0, "h_etaring_iter0", "#eta-ring, 1st iter");
    histData->addComponent(h_etaring_iter2, "h_etaring_iter2", "#eta-ring, 3rd iter");
    histData->addComponent(h_etaring_iter4, "h_etaring_iter4", "#eta-ring, 5th iter");
  }

  vector <TH1*> hlist(histData->getHistPtrs());
  vector<string> hname(histData->getHistNames());
  vector<string> legentry(histData->getHistLegends());

  //vector <TH1*> hlist;
  // hlist.push_back(h_xtal_iter0);
  // hlist.push_back(h_xtal_iter2);
  // hlist.push_back(h_xtal_iter4);
  // if (not singleXtalOnly) {
  //   hlist.push_back(h_etaring_iter0);
  //   hlist.push_back(h_etaring_iter2);
  //   hlist.push_back(h_etaring_iter4);
  // }

  //vector<string> hname;
  // hname.push_back("h_xtal_iter0");
  // hname.push_back("h_xtal_iter2");
  // hname.push_back("h_xtal_iter4");
  // if (not singleXtalOnly) {
  //   hname.push_back("h_etaring_iter0");
  //   hname.push_back("h_etaring_iter2");
  //   hname.push_back("h_etaring_iter4");
  // }

  //vector<string> legentry;
  // legentry.push_back("1 xtal, 1st iter");
  // legentry.push_back("1 xtal, 3rd iter");
  // legentry.push_back("1 xtal, 5th iter");
  // if (not singleXtalOnly) {
  //   legentry.push_back("#eta-Ring, 1st iter");
  //   legentry.push_back("#eta-Ring, 3rd iter");
  //   legentry.push_back("#eta-Ring, 5th iter");
  // }


  for (uint i = 0; i < hlist.size(); i++) {

    hlist[i]->SetNameTitle(hname[i].c_str(),"");
    string niter = hname[i].substr(hname[i].find("iter")+4,string::npos);
    string filename = Form("%s%s/iter_%s/%s_epsilonPlots.root",
			   filePath.c_str(),
			   dirName.c_str(),
			   niter.c_str(),
			   dirName.c_str()); 
    string plotType = hname[i].substr(hname[i].find("h_")+2,hname[i].find("_iter")-2);
    cout << "plotType = " << plotType << endl;
    doPi0MassWithFit(hlist[i],hname[i],isEB, plotType, outDir, lumi, filename);
    if (!hlist[i] || hlist[i]==NULL) {
      cout << "Error: hlist[" << i << "] is NULL, please check. Abort" << endl;
      exit(EXIT_FAILURE);
      
    }
    
  }

  cout << "Now going to plot all histograms together" << endl;
  string canvasname = isEB ? "pi0mass_comparison_EB" : "pi0mass_comparison_EE";
  if (singleXtalOnly) canvasname += "_singleXtal";

  //draw_nTH1(hlist,"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.24","a.u.", canvasname, outDir, legentry, "", lumi, 1, false, false);
  draw_nTH1(histData->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.24","a.u.", canvasname, outDir, histData->getHistLegends(), "", lumi, 1, false, false);

}
