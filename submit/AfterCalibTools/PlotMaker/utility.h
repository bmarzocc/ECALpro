#ifndef utility_h
#define utility_h

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file
#include <string>
#include <vector>
#include <map>
#include <iomanip> //for input/output manipulators

#include <algorithm>  // to use the "reverse" function to reverse the order in the array
#include <Rtypes.h> // to use kColor

//ROOT header files
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TTreeIndex.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TStyle.h>
#include <TString.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

using namespace std;

static string PhpToCopy = "/afs/cern.ch/user/m/mciprian/public/index.php";

//======================================================

void createPlotDirAndCopyPhp(const string& outputDIR) {

  if (outputDIR != "./") {
    system(("mkdir -p " + outputDIR).c_str());
    system(("cp "+ PhpToCopy + " " + outputDIR).c_str());
  }

}

//======================================================                                                                                                                     

void myAddOverflowInLastBin(TH1 *h) {

  Int_t lastBinNumber = h->GetNbinsX();
  Int_t overflowBinNumber = 1 + lastBinNumber;
  Double_t lastBinContent = h->GetBinContent(lastBinNumber);
  Double_t overflowBinContent = h->GetBinContent(overflowBinNumber);
  Double_t lastBinError = h->GetBinError(lastBinNumber);
  Double_t overflowBinError = h->GetBinError(overflowBinNumber);

  // add content of overflow bin in last bin and set error as square root of sum of error squares (with the assumption that they are uncorrelated)                           
  h->SetBinContent(lastBinNumber, lastBinContent + overflowBinContent);
  h->SetBinError(lastBinNumber, sqrt(lastBinError * lastBinError + overflowBinError * overflowBinError));
  // deleting content of last bin (safer, since I might be using that bin to add it again somewhere and I want it to be empty)                                               
  h->SetBinContent(overflowBinNumber,0.0);
  h->SetBinError(overflowBinNumber,0.0);

}


//======================================================                                                                                                                     

void myRebinHisto(TH1 *h, const Int_t rebinFactor = 1) {

  if (rebinFactor != 1) {
    h->Rebin(rebinFactor);
    if ( (h->GetNbinsX() % rebinFactor) != 0) myAddOverflowInLastBin(h);
  }

}



//=============================================================

void draw_nTH1(const vector<TH1*>& vecHist1d = {}, 
	       const string& xAxisNameTmp = "", 
	       const string& yAxisName = "Events", 
	       const string& canvasName = "default", 
	       const string& outputDIR = "./", 
	       const vector<string>& vecLegEntry = {""},
	       const string& ratioPadYaxisName = "var/nominal",
	       const Double_t lumi = -1.0, 
	       const Int_t rebinFactor = 1, 
	       const Bool_t drawPlotLogY = true,
	       const Bool_t drawRatioWithNominal = false  // to be implemented
	       ) 
{

  // assume the "nominal histogram is the first one

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 

  string xAxisName = "";
  string separator = "::";
  Bool_t setXAxisRangeFromUser = false;
  Double_t xmin = 0;
  Double_t xmax = 0;

  size_t pos = xAxisNameTmp.find(separator);
  if (pos != string::npos) {
    string xrange = "";
    setXAxisRangeFromUser = true;
    xAxisName.assign(xAxisNameTmp, 0, pos); 
    xrange.assign(xAxisNameTmp, pos + separator.size(), string::npos);
    separator = ",";
    pos = xrange.find(separator);
    string numString = ""; 
    numString.assign(xrange,0,pos);
    xmin = std::stod(numString);
    numString.assign(xrange,pos + separator.size(), string::npos);
    xmax = std::stod(numString);
  } else {
    xAxisName = xAxisNameTmp;
  }

  // cout << "xAxisName = " << xAxisName << "   xmin = " << xmin << "  xmax = " << xmax << endl;

  for (UInt_t i = 0; i < vecHist1d.size(); i++) {
    myRebinHisto(vecHist1d[i],rebinFactor);
    if (yAxisName == "a.u.") vecHist1d[i]->Scale(1./vecHist1d[i]->Integral());
    vecHist1d[i]->SetStats(0);
  }


  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  if (drawRatioWithNominal) canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  TH1* frame =  (TH1*) vecHist1d[0]->Clone("frame");
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->SetStats(0);

  Int_t colorList[] = {kBlack, kBlue, kRed, kGreen+2, kOrange+1, kCyan, kGreen, kCyan+2, kGray+1, kViolet, kYellow+2};
  vector<Int_t> histColor;
  for (UInt_t i = 0; i < vecHist1d.size(); i++) {   // now color are assigned in reverse order (the main contribution is the last object in the sample array)         
    vecHist1d[i]->SetLineColor(colorList[i]);
    vecHist1d[i]->SetLineWidth(2);
    vecHist1d[i]->SetFillColor(0);
  }

  if (drawRatioWithNominal) {
    vecHist1d[0]->GetXaxis()->SetLabelSize(0);
    vecHist1d[0]->GetXaxis()->SetTitle(0);
  } else {
    vecHist1d[0]->GetXaxis()->SetTitle(xAxisName.c_str());
    // vecHist1d[0]->GetXaxis()->SetTitleOffset(0.8);
    vecHist1d[0]->GetXaxis()->SetLabelSize(0.04);
    vecHist1d[0]->GetXaxis()->SetTitleSize(0.05);    
  }
  vecHist1d[0]->GetYaxis()->SetTitle(yAxisName.c_str());
  vecHist1d[0]->GetYaxis()->SetTitleOffset(1.1);
  // vecHist1d[0]->GetYaxis()->SetTitleOffset(0.8);  // was 1.03 without setting also the size
  vecHist1d[0]->GetYaxis()->SetTitleSize(0.05);
  //vecHist1d[0]->GetYaxis()->SetRangeUser(0.0, max(vecHist1d[0]->GetMaximum(),h2->GetMaximum()) * 1.2);

  //////////////////////////////
  // set X and Y axis range

  // search for maximum Y and for minimum > 0 (latter only if using log scale for Y axis
  Double_t maxY = -999.0;
  for (UInt_t i = 0; i < vecHist1d.size(); i++) {
    if ( vecHist1d[i]->GetBinContent(vecHist1d[i]->GetMaximumBin()) > maxY ) maxY = vecHist1d[i]->GetBinContent(vecHist1d[i]->GetMaximumBin());
  }

  Double_t minY = 1e34;

  if (drawPlotLogY) {

    // quick check if there are no empty bins
    for (UInt_t i = 0; i < vecHist1d.size(); i++) {
      if ( vecHist1d[i]->GetBinContent(vecHist1d[i]->GetMinimumBin()) < minY ) minY = vecHist1d[i]->GetBinContent(vecHist1d[i]->GetMinimumBin());
    }

    if (fabs(minY) < 0.00000001) {

      minY = 1e34;
      
      for (UInt_t ihist = 0; ihist < vecHist1d.size(); ihist++) {

	for (Int_t ibin = 0; ibin <= vecHist1d[ihist]->GetNbinsX(); ibin++ ) {
	  if (vecHist1d[ihist]->GetBinContent(ibin) > 0.0000001 && minY > vecHist1d[ihist]->GetBinContent(ibin)) minY = vecHist1d[ihist]->GetBinContent(ibin);
	}      
      
      }

    }

  }

  vecHist1d[0]->GetYaxis()->SetRangeUser(0.0, maxY * 1.2);

  if (setXAxisRangeFromUser) vecHist1d[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  //////////////////////

  vecHist1d[0]->Draw("Hist");
  vecHist1d[0]->SetFillColor(0);
  vecHist1d[0]->SetMarkerStyle(0);
  for (UInt_t i = 1; i < vecHist1d.size(); i++) {
    vecHist1d[i]->Draw("hist same");
  }

  double legLowY = 0.75;
  if (vecHist1d.size() > 4) legLowY = max( 0.5, legLowY - 0.4 * (vecHist1d.size() - 4) );
  TLegend leg (0.58,0.70,0.95,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  for (UInt_t i = 0; i < vecHist1d.size(); i++) {
    leg.AddEntry(vecHist1d[i],vecLegEntry[i].c_str(),"L");
  }
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  //  CMS_lumi(canvas,Form("%.1f",lumi));
  bool cmsPreliminaryIsUp = false;
  if (yAxisName == "a.u.") cmsPreliminaryIsUp = true;
  if (canvasName.find("pi0mass_comparison_E") != string::npos) cmsPreliminaryIsUp = false;

  if (lumi < 0) CMS_lumi(canvas,"",cmsPreliminaryIsUp,false);
  else CMS_lumi(canvas,Form("%.1f",lumi),cmsPreliminaryIsUp,false);
  setTDRStyle();

  if (drawRatioWithNominal) {
    pad2->Draw();
    pad2->cd();
  
    frame->Reset("ICES");
    if (canvasName.find("comparisonMassVariation") != string::npos) {
      frame->GetYaxis()->SetRangeUser(0.99, 1.01);
      /* if      (outputDIR.find("/eta_0/") != string::npos) frame->GetYaxis()->SetRangeUser(0.99, 1.01); */
      /* else if (outputDIR.find("/eta_1/") != string::npos) frame->GetYaxis()->SetRangeUser(0.98, 1.02); */
      /* else if (outputDIR.find("/eta_2/") != string::npos) frame->GetYaxis()->SetRangeUser(0.98, 1.02); */
    } else if (canvasName.find("elescale") != string::npos) {
      frame->GetYaxis()->SetRangeUser(0.99,1.01);
    } else if (canvasName.find("elescale") != string::npos) {
      frame->GetYaxis()->SetRangeUser(0.99,1.01);
    } else frame->GetYaxis()->SetRangeUser(0.9,1.1);
    frame->GetYaxis()->SetNdivisions(5);
    frame->GetYaxis()->SetTitle(ratioPadYaxisName.c_str());
    frame->GetYaxis()->SetTitleOffset(1.2);
    // frame->GetYaxis()->SetTitleSize(0.15);
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetTitle(xAxisName.c_str());
    if (setXAxisRangeFromUser) frame->GetXaxis()->SetRangeUser(xmin,xmax);
    // frame->GetXaxis()->SetTitleOffset(0.8);
    frame->GetXaxis()->SetTitleSize(0.05);

    vector<TH1D*> ratio;
    for (UInt_t ivar = 1; ivar < vecHist1d.size(); ivar++) 
      ratio.push_back( (TH1D*) vecHist1d[ivar]->Clone(Form("ratio_%d",ivar)) );

    TH1D* den_noerr = (TH1D*) vecHist1d[0]->Clone("den_noerr");
    TH1D* den = (TH1D*) vecHist1d[0]->Clone("den");
    for(int iBin = 1; iBin < den->GetNbinsX()+1; iBin++)
      den_noerr->SetBinError(iBin,0.);

    den->Divide(den_noerr);
    den->SetFillColor(kGray);
    frame->Draw();
    den->Draw("E2same");
    for (UInt_t ir = 0; ir < ratio.size(); ir++) {
      ratio[ir]->Divide(den_noerr);
      // ratio[ir]->SetMarkerSize(0.65);
      // ratio[ir]->Draw("EPsame");
      ratio[ir]->SetMarkerStyle(0);
      ratio[ir]->SetLineWidth(2);
      ratio[ir]->Draw("Hist same");
    }
 

    TF1* line = new TF1("horiz_line","1",den->GetXaxis()->GetBinLowEdge(1),den->GetXaxis()->GetBinLowEdge(den->GetNbinsX()+1));
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw("Lsame");
    // for (UInt_t ir = 0; ir < ratio.size(); ir++)
    //   ratio[ir]->Draw("EPsame");
    pad2->RedrawAxis("sameaxis");

  }  // end of ratio plot settings

  if (canvasName.find("tmpToBeRemoved") == string::npos) { 
    canvas->SaveAs((outputDIR + canvasName + ".png").c_str());
    canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str());
  }

  if (drawPlotLogY) {

    if (yAxisName == "a.u.") vecHist1d[0]->GetYaxis()->SetRangeUser(minY*0.05, maxY*100);
    else vecHist1d[0]->GetYaxis()->SetRangeUser(minY*0.05, maxY*100);
    canvas->SetLogy();
    /* if (lumi < 0) CMS_lumi(canvas,"",true,false); */
    /* else CMS_lumi(canvas,Form("%.1f",lumi),true,false); */
    if (canvasName.find("tmpToBeRemoved") == string::npos) { 
      canvas->SaveAs((outputDIR + canvasName + "_logY.png").c_str());
      canvas->SaveAs((outputDIR + canvasName + "_logY.pdf").c_str());
    }
    canvas->SetLogy(0);

  }    

  delete canvas;

}



//=============================================================


#endif
