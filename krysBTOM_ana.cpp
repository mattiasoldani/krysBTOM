// to compile: g++ krysBTOM_ana.cpp -o krysBTOM_ana.out `root-config --cflags --libs`

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <TH1.h>
#include <TFile.h>
#include <math.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TError.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
// global settings ////////////////////////////////////////////////////////////////////////////////

// the only things to be set locally in functions:
// - data structure, which is in krysBTOM_ana() (check treeStructure)
// - several individual-histogram parameters, which are in krysBTOM_ana()
// - also don't hesitate to create new tree aliases (histograms/canvases) in treeConditioning() (krysBTOM_ana())

// base element numbers for each setup section

	const int nSi = 8; // nr. of Si tracking layers
	const int nGonio = 5; // nr. of goniometer DOFs
	const int nDigi = 16; // nr. of digitizer channels
	/*// bent TEST values:
	const int nSi = 8; // nr. of Si tracking layers
	const int nGonio = 3; // nr. of goniometer DOFs
	const int nDigi = 0; // nr. of digitizer channels
	*/
	/*// straight TEST values:
	const int nSi = 8; // nr. of Si tracking layers
	const int nGonio = 5; // nr. of goniometer DOFs
	const int nDigi = 16; // nr. of digitizer channels
	*/

// booleans for histogram and canvas activation (true) or deactivation (false)

	//bool b_DESY22 = false; // DESY22, ADD BOOLEAN!

	bool b_plotLastSpill = true; // general switch for last-spill histograms (creation & drawing in canvases - dedicated tree is not even created)
	bool b_outTracking = false; // general switch for output tracking (controls all the calculations and plots that make use of the output tracking layers)

	bool b_h1_nHit = true; // multiplicities 1D
	bool b_h1_xRaw = true; // beam profiles 1D
	bool b_h1_thIn = true; // input angles 1D (both unaligned & aligned)
	bool b_h1_thOut = false; // output angles 1D (both unaligned & aligned)
	bool b_h1_thDelta = false; // deflection angles 1D
	bool b_h1_time_noCuts = true; // digitizer times 1D, with no cuts
	bool b_h1_ph_noCuts = true; // digitizer PHs 1D, with no cuts
	bool b_h1_phEq_noCuts = false; // digitizer PHs 1D, equalised, with no cuts
	bool b_h1_phEq_counters = true; // digitizer PHs 1D, equalised, only multiplicity counters
		// also check the boolean on the actual availability of these data below
	bool b_h1_phEq_fwdCalo = true; // digitizer PHs 1D, equalised, only forward calorimeter (single channels)
		// also check the boolean on the actual availability of these data below
	bool b_h1_phEq_latCalo = false; // digitizer PHs 1D, equalised, only lateral calorimeter (single channels)
		// also check the boolean on the actual availability of these data below
	bool b_h1_phEq_cry = false; // digitizer PHs 1D, equalised, only read-out crystal
		// also check the boolean on the actual availability of these data below
	bool b_h1_phTot_fwdCalo = true;  // forward calorimeter total PH 1D
		// also check the boolean on the actual availability of these data below
	bool b_h1_phTot_latCalo = false;  // forward calorimeter total PH 1D
		// also check the boolean on the actual availability of these data below
	
	bool b_h2_iStep_thDelta = false; // 2D correlation between deflection angles and step nr.
	bool b_h2_xGonio_thDelta = false; // 2D correlation between deflection angles and gonio. DOFs (alltogether)
	bool b_h2_xCryX_thDelta = false; // 2D correlation between deflection angles and horizontal position @ crystal
	bool b_h2_xCryY_thDelta = false; // 2D correlation between deflection angles and vertical position @ crystal
	bool b_h2_thInX_thDelta = false; // 2D correlation between deflection angles and horizontal input angle @ crystal
	bool b_h2_thInY_thDelta = false; // 2D correlation between deflection angles and vertical input angle @ crystal
	bool b_h2_time_ph_noCuts = true; // digitizer times vs. PHs, with no cuts
	
	bool b_plotsPerSingleStep = false; // general switch for single-step histograms and canvases
		// note: creation of individual objects can be toggled below
		
	bool b_h1_thDelta_SingleStep = false; // single-step deflection angles 1D
	
	bool b_h2_xCryX_thDelta_SingleStep = false; // single-step 2D correlation between deflection angles and horizontal position @ crystal
	bool b_h2_xCryY_thDelta_SingleStep = false; // single-step 2D correlation between deflection angles and vertical position @ crystal
	bool b_h2_thInX_thDelta_SingleStep = false; // single-step 2D correlation between deflection angles and horizontal input angle @ crystal
	bool b_h2_thInY_thDelta_SingleStep = false; // single-step 2D correlation between deflection angles and vertical input angle @ crystal
	
// conditions for booleans, ranges and binnings for histograms

	int nBins_nHit = 10; // nr. of bins in multiplicity plots
	int nBins_xRaw = 100; // nr. of bins in profile plots
	int nBins_thIn = 500; // nr. of bins in input angle plots
	int nBins_thOut = 500; // nr. of bins in output angle plots
		// note: also used for deflection angle binning
	int nBinsX_xGonio_thDelta = 400; // nr. of horizontal bins in deflection angle vs. goniometer DOF plots
	int nBinsX_xCry_thDelta = 200; // nr. of horizontal bins in deflection angle vs. input positions @ crystal plots
	int nBinsX_thIn_thDelta = 100; // nr. of horizontal bins in deflection angle vs. input angles plots

	double r_xCryMin = -1.0; // lower range for crystal spot plots, both hor. & ver.
	double r_xCryMax = 3.0; // upper range for crystal spot plots, both hor. & ver.
	double r_thOut = 1000.0; // range for output angle plots; 
		// note: also used for output & deflection angle cut definition
		// note: also used for deflection angle range
	double cut_xCryXMin = -10.0; // crystal fiducial area boundaries, min. x
	double cut_xCryXMax = 10.0; // crystal fiducial area boundaries, max. x
	double cut_xCryYMin = -10.0; // crystal fiducial area boundaries, min. y
	double cut_xCryYMax = 10.0; // crystal fiducial area boundaries, max. y
	double cut_thIn = 10000.0; // cut on input angles
	double thInCutMultForOut = 1.0; // factor used for inverse cut on deflection angle for crystal surface identification (multiplied by thInCutMultForOut)
	
	int nBins_time = 128; // nr. of bins in digi. time plots
	int nBins_ph = 1024; // nr. of bins on PH plots
	
	double r_timeMax = 512.0; // upper range for digi. time plots
	double r_phMax = 16198.0; // upper range for PH plots
	double cut_digiTimeMin[nDigi] = { // cuts on digi. time, lower value per channel
		0., 0., 0., 0., 0., 0., 0., 0., 
		0., 0., 0., 0., 0., 0., 0., 0.
	};
	double cut_digiTimeMax[nDigi] = { // cuts on digi. time, upper value per channel
		600, 600, 600, 600, 600, 600, 600, 600, 
		600, 600, 600, 600, 600, 600, 600, 600
	};
	double cut_digiPHEqMin[nDigi] = { // cuts on digi. PH (equalised), lower value per channel
		0., 0., 0., 0., 0., 0., 0., 0., 
		0., 0., 0., 0., 0., 0., 0., 0.
	};
	double cut_digiPHEqMax[nDigi] = { // cuts on digi. PH (equalised), upper value per channel
		20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 
		20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 20e3
	};
	/*// bent TEST values:
	double cut_digiTimeMin[nDigi] = {}; // cuts on digi. time, lower value per channel
	double cut_digiTimeMax[nDigi] = {}; // cuts on digi. time, upper value per channel
	double cut_digiPHEqMin[nDigi] = {}; // cuts on digi. PH (equalised), lower value per channel
	double cut_digiPHEqMax[nDigi] = {}; // cuts on digi. PH (equalised), upper value per channel
	*/
	/*// straight TEST values:
	double cut_digiTimeMin[nDigi] = { // cuts on digi. time, lower value per channel
		0., 0., 0., 0., 0., 0., 0., 0., 
		0., 0., 0., 0., 0., 0., 0., 0.
	};
	double cut_digiTimeMax[nDigi] = { // cuts on digi. time, upper value per channel
		600, 600, 600, 600, 600, 600, 600, 600, 
		600, 600, 600, 600, 600, 600, 600, 600
	};
	double cut_digiPHEqMin[nDigi] = { // cuts on digi. PH (equalised), lower value per channel
		0., 0., 0., 0., 0., 0., 0., 0., 
		0., 0., 0., 0., 0., 0., 0., 0.
	};
	double cut_digiPHEqMax[nDigi] = { // cuts on digi. PH (equalised), upper value per channel
		20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 
		20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 20e3, 20e3
	};
	*/
	
// settings on setup & other data conditioning parameters

	int iIn[4] = {0, 1, 2, 3}; // indexes of input tracking layers, always 4 entries: x0, y0, x1, y1
	double zIn[4] = {0, 0, 78.9, 78.9}; // z of input tracking layers, always 4 entries: x0, y0, x1, y1
	/*// bent TEST values:
	int iIn[4] = {0, 1, 2, 3}; // indexes of input tracking layers, always 4 entries: x0, y0, x1, y1
	double zIn[4] = {0, 0, 500.0, 500.0}; // z of input tracking layers, always 4 entries: x0, y0, x1, y1
	*/
	/*// straight TEST values:
	int iIn[4] = {0, 1, 2, 3}; // indexes of input tracking layers, always 4 entries: x0, y0, x1, y1
	double zIn[4] = {0, 0, 41.0, 41.0}; // z of input tracking layers, always 4 entries: x0, y0, x1, y1
	*/

	bool bThOutFromCry = true; // if true (false) output angle is computed w/ (w/o) input projection @ crystal
		// note: it is subject to b_outTracking
		// note: if no output angle analysis is needed, set bThOutFromCry to false and 4 dummy entries in zOut & iOut
	double zCry = 78.9+2.25+24.0; // z of crystal centre
	/*// bent TEST values:
	bool bThOutFromCry = false; // if true (false) output angle is computed w/ (w/o) input projection @ crystal
		// note: it is subject to b_outTracking
		// note: if no output angle analysis is needed, set bThOutFromCry to false and 4 dummy entries in zOut & iOut
	double zCry = 600.0; // z of crystal centre
	*/
	/*// straight TEST values:
	bool bThOutFromCry = true; // if true (false) output angle is computed w/ (w/o) input projection @ crystal
		// note: it is subject to b_outTracking
		// note: if no output angle analysis is needed, set bThOutFromCry to false and 4 dummy entries in zOut & iOut
	double zCry = 63.75; // z of crystal centre
	*/

	int iOut[6] = {
		4, 5
	}; // indexes of output tracking layers, any even nr. of entries: x0, y0, ...
	double zOut[6] = {
		78.9+100.15, 78.9+100.15
	}; // z of output tracking layers, any even nr. of entries: x0, y0, ...
		// note: it is subject to b_outTracking
		// note: if bThOutFromCry layers 0 & 1 are used with projections @ crystal, otherwise layers 0, 1, 2 & 3 are used
		// note: if no output angle analysis is needed, set bThOutFromCry to false and 4 dummy entries in zOut & iOut
	/*// bent TEST values:
	int iOut[6] = {
		4, 5,
		6, 7
	}; // indexes of output tracking layers, any even nr. of entries: x0, y0, ...
	double zOut[6] = {
		700.0, 700.0,
		800.0, 800.0
	}; // z of output tracking layers, any even nr. of entries: x0, y0, ...
	*/
	/*// straight TEST values:
	int iOut[6] = {
		4, 5
	}; // indexes of output tracking layers, any even nr. of entries: x0, y0, ...
	double zOut[6] = {
		106.65, 106.65
	}; // z of output tracking layers, any even nr. of entries: x0, y0, ...
	*/

	double thInCentres[2] = {-1792.0, 5749.0}; // centres of the raw input angle distributions (in urad)
	double thOutCentres[2] = {0.0, 0.0}; // centres of the raw input angle distributions (in urad)
	/*// bent TEST values:
	double thInCentres[2] = {-837.3, -306.1}; // centres of the raw input angle distributions (in urad)
	double thOutCentres[2] = {-1372.0, -725.2}; // centres of the raw input angle distributions (in urad)
	*/
	/*// straight TEST values:
	double thInCentres[2] = {-5561.16, -4270.23}; // centres of the raw input angle distributions (in urad)
	double thOutCentres[2] = {0.0, 0.0}; // centres of the raw input angle distributions (in urad)
	*/

	const char * xGonioCorrect[nGonio] = {"", "", "", "", ""}; 
		// beam parameter (xCryX/Y, thInX/Y and rescaling factors) to be subtracted from each raw goniometer DOF; there must be exactly nGonio entries, if in doubt insert empty string
	/*// bent TEST values:
	const char * xGonioCorrect[nGonio] = {"", "", ""}; 
	*/
	/*// straight TEST values:
	const char * xGonioCorrect[nGonio] = {"", "", "", "", ""}; 
	*/
	
	// equalisation factors for digitizer-related detectors
	// the order is the same as in the raw data, nothing to do with the lists of digitizer sets of channels
	// there must be exactly nDigi entries, if in doubt insert ones
	double digiEq[nDigi] = {
		1, 1, 0.461320, 0.618319, 0.580172, 0.559943, 1, 0.882152,
		0.553153, 1.063839, 1, 1, 1, 1, 1, 1
	};
	/*// bent TEST values:
	double digiEq[nDigi] = {};
	*/
	/*// straight TEST values:
	double digiEq[nDigi] = {
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1
	};
	*/
	
	// digitizer-related detectors availability booleans
	// these are for data availability, also set the boolean for histograms and canvases above
	bool b_counters = true; // all the multiplicity counters & vetoes
	bool b_fwdCalo = true; // forward calorimeter
	bool b_latCalo = false; // lateral calorimeter
	bool b_cry = false; // crystal, if read out
	/*// bent TEST values:
	bool b_counters = false; // all the multiplicity counters & vetoes
	bool b_fwdCalo = false; // forward calorimeter
	bool b_latCalo = false; // lateral calorimeter
	bool b_cry = false; // crystal, if read out
	*/
	/*// straight TEST values:
	bool b_counters = true; // all the multiplicity counters & vetoes
	bool b_fwdCalo = true; // forward calorimeter
	bool b_latCalo = false; // lateral calorimeter
	bool b_cry = false; // crystal, if read out
	*/
	
	// lists of digitizer sets of channels (order matters for plotting)
	vector<int> ls_counters = {0, 1}; // all the multiplicity counters & vetoes
	vector<int> ls_fwdCalo = {2, 3, 4, 5, 6, 7, 8, 9, 10}; // forward calorimeter
	vector<int> ls_latCalo = {}; // lateral calorimeter
	vector<int> ls_cry = {}; // crystal, if read out
	// bent TEST values: whatever
	/*// straight TEST values:
	vector<int> ls_counters = {11, 12}; // all the multiplicity counters & vetoes
	vector<int> ls_fwdCalo = {0, 1, 2, 3, 4, 5, 6, 7, 8}; // forward calorimeter
	vector<int> ls_latCalo = {}; // lateral calorimeter
	vector<int> ls_cry = {}; // crystal, if read out
	*/
	
// miscellaneous graphical settings for canvases

	int cLineTotal = 1; // line color for multi-spill histograms
	int cLineLast = 4; // line color for last-spill histograms
	int cFillTotal = 17; // area color for multi-spill histograms
	
	int optFit = 10001; // code of info to be displayed in fit boxes
	int optStat1D = 1001111; // code of info to be displayed in statistics boxes (1D & default)
	int optStat2D = 1110011; // code of info to be displayed in statistics boxes (2D)
	
	// digitizer-related detector canvas shapes
	int cShape_counters[2] = {1, 2}; // all the multiplicity counters & vetoes
	int cShape_fwdCalo[2] = {3, 3}; // forward calorimeter
	int cShape_latCalo[2] = {0, 0}; // lateral calorimeter
	int cShape_cry[2] = {0, 0}; // crystal, if read out
	// bent TEST values: whatever
	/*// straight TEST values:
	int[2] cShape_counters = {3, 1}; // all the multiplicity counters & vetoes
	int[2] cShape_fwdCalo = {3, 3}; // forward calorimeter
	int[2] cShape_latCalo = {0, 0}; // lateral calorimeter
	int[2] cShape_cry = {0, 0}; // crystal, if read out
	*/
	
// miscellaneous I/O settings
	
	bool bWriteOutTreeFile = false; // write raw tree file on disk?
	
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// function to perform calculations and create new variables in trees /////////////////////////////

void treeConditioning(TTree * tree, string name){
		
// create all derived variables in the trees

	// tracking
	for(int i=0; i<2; i++){
		string direction = i==0 ? Form("X") : Form("Y");
		
		// input raw angles
		tree->SetAlias(Form("thInRaw%s", direction.c_str()), Form("((xRaw[%d] - xRaw[%d]) / (%f)) * 1e6", iIn[2+i], iIn[i], zIn[2+i] - zIn[i]));
		
		// input aligned angles
		tree->SetAlias(Form("thIn%s", direction.c_str()), Form("thInRaw%s - %f", direction.c_str(), thInCentres[i]));
		
		// input track projection to crystal
		tree->SetAlias(Form("xCry%s", direction.c_str()), Form("xRaw[%d] + (%f) * thIn%s * 1e-6", iIn[i], zCry - zIn[i], direction.c_str()));
		
		// output raw angles
		if(b_outTracking){
			tree->SetAlias(Form("thOutRaw%s", direction.c_str()),
				bThOutFromCry ? 
				Form("((xRaw[%d] - xCry%s) / (%f)) * 1e6", iOut[i], direction.c_str(), zOut[i] - zCry) :
				Form("((xRaw[%d] - xRaw[%d]) / (%f)) * 1e6", iOut[2+i], iOut[i], zOut[2+i] - zOut[i])
			);
		};
		
		// output aligned angles
		if(b_outTracking){
			tree->SetAlias(Form("thOut%s", direction.c_str()), Form("thOutRaw%s - %f", direction.c_str(), thOutCentres[i]));
		};
		
		// deflection angles
		if(b_outTracking){
			tree->SetAlias(Form("thDelta%s", direction.c_str()), Form("thOut%s - thIn%s", direction.c_str(), direction.c_str()));
		};
	}
	
	// goniometer data conditioning (manually combine them with input tracking data)
	for(int i=0; i<nGonio; i++){
		string xGonioCorrectTemp = Form("%s", xGonioCorrect[i]);
		tree->SetAlias(Form("xGonio%d", i), xGonioCorrectTemp == Form("") ? Form("xGonioRaw[%d]", i) : Form("xGonioRaw[%d] %s", i, xGonioCorrectTemp.c_str()));
	}
	
	// digitizer-related detectors equalisation
	for(int i=0; i<nDigi; i++){
		if(nDigi!=0){tree->SetAlias(Form("digiPHEq%d", i), Form("%f * digiPH[%d]", digiEq[i], i));};
	}
	
	// calorimeters total signal - forward
	if(b_fwdCalo){
		string strCaloFwdSum = Form("");
		for(auto i=begin(ls_fwdCalo); i!=end(ls_fwdCalo); i++){
			if(i==begin(ls_fwdCalo)){strCaloFwdSum = Form("digiPHEq%d", *i);}
			else{strCaloFwdSum = Form("%s + digiPHEq%d", strCaloFwdSum.c_str(), *i);}
		}
		tree->SetAlias(Form("PHCaloFwd"), Form("%s", strCaloFwdSum.c_str()));
	}
	
	// calorimeters total signal - lateral
	if(b_latCalo){
		string strCaloLatSum = Form("");
		for(auto i=begin(ls_latCalo); i!=end(ls_latCalo); i++){
			if(i==begin(ls_latCalo)){strCaloLatSum = Form("digiPHEq%d", *i);}
			else{strCaloLatSum = Form("%s + digiPHEq%d", strCaloLatSum.c_str(), *i);}
		}
		tree->SetAlias(Form("PHCaloLat"), Form("%s", strCaloLatSum.c_str()));
	}
	
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	// custom tree aliases here below //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	
// if requested, write tree to output file
	if(bWriteOutTreeFile){tree->Write(name.c_str());}
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// function to draw canvases related to specific digitizer-related detectors //////////////////////

	void drawSpecDigiCanvas(TFile * outHistFile, int nEntries, int nEntriesLast, bool b_plotLastSpill, string outPdfPath, bool b, string detNameShort, string detNameLong, int * cShape, vector<int> ls){

		if(b){
			int nCh;
			TH1F *h1TempOut[400];
			
			string canvasName = Form("c_digiPhEq_%s", detNameShort.c_str());
			string canvasTitle = Form("%s equalised PHs", detNameLong.c_str());
			TCanvas * canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
			canvas->Divide(cShape[0], cShape[1]);

			nCh = ls.size();
			int j=0;
			for(auto i=begin(ls); i!=end(ls); i++){
				j+=1;
				canvas->cd(j);
				gPad->SetLogy(); // log scale in y?
				
				h1TempOut[j] = (TH1F*)outHistFile->Get(Form("h1_phEq%d", *i));
				h1TempOut[j]->SetLineColor(cLineTotal);
				h1TempOut[j]->SetFillColor(cFillTotal);
				h1TempOut[j]->Draw("HIST");

				if(b_plotLastSpill&&(nEntriesLast!=0)){
					h1TempOut[nCh+j] = (TH1F*)outHistFile->Get(Form("h1_phEq%d_last", *i));
					h1TempOut[nCh+j]->Scale(h1TempOut[j]->Integral() / h1TempOut[nCh+j]->Integral());
					h1TempOut[nCh+j]->SetLineColor(cLineLast);
					h1TempOut[nCh+j]->Draw("HIST SAME");
				}
			}
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}else{return;}
		
	}
		
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// main ///////////////////////////////////////////////////////////////////////////////////////////

void krysBTOM_ana(){ // this partakes everything that should go in the main (which is at the bottom of the file)
	gErrorIgnoreLevel = kWarning;

// data structure to be used in trees
	string treeStructure = Form("xRaw[%d]/D:nStrHit[%d]/I:nHit[%d]/I:digiBase[%d]/I:digiPH[%d]/I:digiTime[%d]/I:xGonioRaw[%d]/D:iSpill/I:iStep/I:iEv/I:xCharge[%d]/D", nSi, nSi, nSi, nDigi, nDigi, nDigi, nGonio, 4);
	/*// bent TEST values:
	string treeStructure = Form("nHit[%d]/I:xRaw[%d]/D:nStrHit[%d]/I:xGonioRaw[%d]/D:iSpill/I:iStep/I:iEv/I", nSi, nSi+4, nSi+4, nGonio);
	*/
	/*// straight TEST values:
	string treeStructure = Form("epoch/I:iEv/I:xRaw[%d]/D:nHit[%d]/I:xGonioRaw[%d]/D:iSpill/I:iStep/I:digiBase[%d]/I:digiPH[%d]/I:digiTime[%d]/I", nSi, nSi, nGonio, nDigi, nDigi, nDigi);
	*/
	
// read data file list
	vector<string> lsFiles;
	{
		ifstream file0("outtext/lsFiles.conf");
		if (file0.is_open()){for(string line; file0 >> line;){lsFiles.push_back(line);}}
		file0.close();
	}
	
// read data path
	string inFileDataPath;
	{
		ifstream file1("inFileDataPath.conf");
		if (file1.is_open()){for(string line; file1 >> line;){inFileDataPath = line;}}
		file1.close();
	}

// loop on runs --> tree creation w/ raw data

	string outTreeFileName = Form("outtree/treeFile.root");
	TFile * outTreeFile;
	if(bWriteOutTreeFile){outTreeFile = new TFile(outTreeFileName.c_str(), "RECREATE");}
	TTree *tree = new TTree("tree", "tree");
	TTree *treeLast = new TTree("treeLast", "treeLast");
	
	for (int iRun=0; iRun < lsFiles.size(); iRun++){ // loop on runs
		cout << "working on file " << lsFiles[iRun] << "..." << endl;
		string inFileDataName = inFileDataPath + lsFiles[iRun];
		ifstream file(inFileDataName);
		bool bLastRun = iRun==(lsFiles.size()-1);
		if(file.is_open()){
			if (bLastRun){treeLast->ReadFile(inFileDataName.c_str(), treeStructure.c_str());}
			else{tree->ReadFile(inFileDataName.c_str(), treeStructure.c_str());}
		}
		file.close();
	} // end of loop on runs
	
	int nEntries = tree->GetEntries();
	int nEntriesLast = treeLast->GetEntries();
	cout << Form("tree created with %d entries (%d in the last spill)", nEntries, nEntriesLast) << endl;
	
// data conditioning --> new variables in trees
// depends on external function treeConditioning
	string name = Form("tree"), nameLast = Form("treeLast");
	treeConditioning(tree, name);
	if(b_plotLastSpill&&(nEntriesLast!=0)){
		treeConditioning(treeLast, nameLast);
	}
	cout << "tree aliases created, now generating histograms..." << endl;
	
// define booleans

	int iOutLength = sizeof(iOut)/sizeof(iOut[0]);

	// data processing errors on single Si layers & for base track (i.e. all layers used for tracking)
	vector<string> b_xRawErrors;
	for(int i=0; i<nSi; i++){
		b_xRawErrors.push_back(Form("abs(xRaw[%d])<%f", i, 50.0));
	}
	string b_xRawBaseTrackErrors = Form("");
	if(b_outTracking){
		int iBaseTrack[8] = {iIn[0], iIn[1], iIn[2], iIn[3], iOut[0], iOut[1], iOut[2], iOut[3]};
		if(bThOutFromCry){iBaseTrack[6] = iOut[0]; iBaseTrack[7] = iOut[1];}
		for(int i=0; i<8; i++){
			if(i==0){b_xRawBaseTrackErrors = Form("(%s)", b_xRawErrors[iBaseTrack[i]].c_str());}
			else{b_xRawBaseTrackErrors = Form("%s(%s)", b_xRawBaseTrackErrors.c_str(), b_xRawErrors[iBaseTrack[i]].c_str());}
			if(i<8-1){b_xRawBaseTrackErrors = Form("%s&&", b_xRawBaseTrackErrors.c_str());}
		}
	}else{
		int iBaseTrack[4] = {iIn[0], iIn[1], iIn[2], iIn[3]};
		for(int i=0; i<4; i++){
			if(i==0){b_xRawBaseTrackErrors = Form("(%s)", b_xRawErrors[iBaseTrack[i]].c_str());}
			else{b_xRawBaseTrackErrors = Form("%s(%s)", b_xRawBaseTrackErrors.c_str(), b_xRawErrors[iBaseTrack[i]].c_str());}
			if(i<4-1){b_xRawBaseTrackErrors = Form("%s&&", b_xRawBaseTrackErrors.c_str());}
		}
	};

	// multiplicities
	string b_singleHitIn = Form("(nHit[%d]==1)&&(nHit[%d]==1)&&(nHit[%d]==1)&&(nHit[%d]==1)", iIn[0], iIn[1], iIn[2], iIn[3]);
	string b_singleHitOut;
	if(b_outTracking){
		vector<string> b_singleHitOutSingle;
		b_singleHitOut = Form("");
		for(int i=0; i<(bThOutFromCry?2:4); i++){b_singleHitOutSingle.push_back(Form("nHit[%d]==1", iOut[i]));}
		for(int i=0; i<(bThOutFromCry?2:4); i++){
			if(i==0){b_singleHitOut = Form("(%s)", b_singleHitOutSingle[i].c_str());}
			else{b_singleHitOut = Form("%s(%s)", b_singleHitOut.c_str(), b_singleHitOutSingle[i].c_str());}
			if(i<(bThOutFromCry?2:4)-1){b_singleHitOut = Form("%s&&", b_singleHitOut.c_str());}
		}
	}else{
		b_singleHitOut = Form("(true)");
	};
	
	// trajectories (angles & projections)
	string b_inAligned = Form("(abs(thInX)<%f)&&(abs(thInY)<%f)", cut_thIn, cut_thIn);
	string b_inCry = Form("(xCryX>%f)&&(xCryX<%f)&&(xCryY>%f)&&(xCryY<%f)", cut_xCryXMin, cut_xCryXMax, cut_xCryYMin, cut_xCryYMax);
	string b_thOutLim, b_thDeltaHigh;
	if(b_outTracking){
		b_thOutLim = Form("(abs(thOutX)<%f)&&(abs(thOutY)<%f)", r_thOut, r_thOut);
		b_thDeltaHigh = Form("(abs(thDeltaX)>%f)||(abs(thDeltaY)<%f)|", thInCutMultForOut*r_thOut, thInCutMultForOut*r_thOut);
	}else{
		b_thOutLim = Form("(true)");
		b_thDeltaHigh = Form("(true)");
	};
	
	// digitizer-related, single-channel
	vector<string> b_digiTime, b_digiPHEq;
	if(nDigi!=0){
		for(int i=0; i<nDigi; i++){
			b_digiTime.push_back(Form("(digiTime[%d]>%f)&&(digiTime[%d]<%f)", i, cut_digiTimeMin[i], i, cut_digiTimeMax[i]));
			b_digiPHEq.push_back(Form("(digiPHEq%d>%f)&&(digiPHEq%d<%f)", i, cut_digiPHEqMin[i], i, cut_digiPHEqMax[i]));
		}
	}
	
	// digitizer-related, total calorimeters
	string b_digiTimeCaloFwd = Form("((digiTime[0]>%f)&&(digiTime[0]<%f))", -20e3, 20e3);
	string b_digiTimeCaloLat = Form("((digiTime[0]>%f)&&(digiTime[0]<%f))", -20e3, 20e3);
	if(b_fwdCalo){
		for(auto i=begin(ls_fwdCalo); i!=end(ls_fwdCalo); i++){			
			b_digiTimeCaloFwd = Form("%s||(%s)", b_digiTimeCaloFwd.c_str(), b_digiTime[*i].c_str());
		}
	}
	if(b_latCalo){
		for(auto i=begin(ls_latCalo); i!=end(ls_latCalo); i++){
			b_digiTimeCaloLat = Form("%s||(%s)", b_digiTimeCaloLat.c_str(), b_digiTime[*i].c_str());
		}
	}
	
// histogram creation, fit & writing

	int nBins, nBinsX, nBinsY;
	double xRangeLeft, xRangeRight, yRangeLeft, yRangeRight;
	string histName;
	string varName, varNameX, varNameY;
	string histNameAdd = Form("_last");
	string cond, cond1;
	TH1F *h1Temp;
	TH2F *h2Temp;
	TProfile *tpTemp;
	
	string outHistFileName = "outhist/histFile.root";
	TFile *outHistFile = new TFile(outHistFileName.c_str(), "RECREATE");

	TCanvas * canvas0 = new TCanvas("", "");	// this is only needed to silent canvas-related printouts

	// histograms, raw tracking system
	for(int i=0; i<nSi; i++){
		
		// - multiplicities (also last-spill version)
		nBins = nBins_nHit ; xRangeLeft = -0.5 ; xRangeRight = nBins-0.5 ;

		cond = Form("%s", b_xRawErrors[i].c_str()); // set total boolean
		
		varName = Form("nHit[%d]", i);
		histName = Form("h1_nHit%d", i);
		
		if(b_h1_nHit){
			tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, xRangeLeft, xRangeRight), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("Si layer %d multiplicity", i));
			h1Temp->GetXaxis()->SetTitle("nr. of hits");
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, xRangeLeft, xRangeRight), cond.c_str(), "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("Si layer %d multiplicity", i));
				h1Temp->GetXaxis()->SetTitle("nr. of hits");
				h1Temp->Write();			
			}
		}
		
		// - beam profiles (also last-spill version)
		nBins = nBins_xRaw;
		
		varName = Form("xRaw[%d]", i);
		histName = Form("h1_xRaw%d", i);
		
		if(b_h1_xRaw){
			tree->Draw(Form("%s>>%s(%d)", varName.c_str(), histName.c_str(), nBins), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("Si layer %d beam profile", i));
			h1Temp->GetXaxis()->SetTitle("position on layer [cm]");
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d)", varName.c_str(), histName.c_str(), nBins), cond.c_str(), "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("Si layer %d beam profile", i));
				h1Temp->GetXaxis()->SetTitle("position on layer [cm]");
				h1Temp->Write();
			}
		}
	}
	
	// histograms, trajectories
	for(int i=0; i<2; i++){
		string direction = i%2==0 ? Form("X") : Form("Y");
		
		// - input raw angles
		nBins = nBins_thIn; // range settings

		cond = Form("%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str()); // set total boolean
		
		varName = Form("thInRaw%s", direction.c_str());
		histName = Form("h1_thInRaw%s", direction.c_str());
		
		if(b_h1_thIn){
			tree->Draw(Form("%s>>%s(%d)", varName.c_str(), histName.c_str(), nBins), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("%s input raw angle", i==0? "horizontal" : "vertical"));
			h1Temp->GetXaxis()->SetTitle("angle [urad]");
			h1Temp->Fit("gaus", "Q"); // gaussian fit
			h1Temp->Write();
		}
		
		// - input aligned angles (also last-spill version)
		nBins = nBins_thIn;

		cond = Form("%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str()); // set total boolean
		
		varName = Form("thIn%s", direction.c_str());
		histName = Form("h1_thIn%s", direction.c_str());
		
		if(b_h1_thIn){
			tree->Draw(Form("%s>>%s(%d)", varName.c_str(), histName.c_str(), nBins), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("%s input angle", i==0? "horizontal" : "vertical"));
			h1Temp->GetXaxis()->SetTitle("angle [urad]");
			h1Temp->Fit("gaus", "Q"); // gaussian fit
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d)", varName.c_str(), histName.c_str(), nBins), cond.c_str(), "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("%s input angle", i==0? "horizontal" : "vertical"));
				h1Temp->GetXaxis()->SetTitle("angle [urad]");
				h1Temp->Write();
			}
		}
		
		// - output raw angles
		nBins = nBins_thOut; // range settings

		cond = Form("%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitOut.c_str()); // set total boolean
		
		varName = Form("thOutRaw%s", direction.c_str());
		histName = Form("h1_thOutRaw%s", direction.c_str());
		
		if(b_outTracking&&b_h1_thOut){
			tree->Draw(Form("%s>>%s(%d)", varName.c_str(), histName.c_str(), nBins), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("%s output raw angle", i==0? "horizontal" : "vertical"));
			h1Temp->GetXaxis()->SetTitle("angle [urad]");
			h1Temp->Fit("gaus", "Q"); // gaussian fit
			h1Temp->Write();
		}
		
		// - output aligned angles (also last-spill version)
		nBins = nBins_thOut; // range settings

		cond = Form("%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitOut.c_str()); // set total boolean
		
		varName = Form("thOut%s", direction.c_str());
		histName = Form("h1_thOut%s", direction.c_str());

		if(b_outTracking&&b_h1_thOut){
			tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, -r_thOut, r_thOut), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("%s output angle", i==0? "horizontal" : "vertical"));
			h1Temp->GetXaxis()->SetTitle("angle [urad]");
			h1Temp->Fit("gaus", "Q"); // gaussian fit
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, -r_thOut, r_thOut), cond.c_str(), "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("%s output angle", i==0? "horizontal" : "vertical"));
				h1Temp->GetXaxis()->SetTitle("angle [urad]");
				h1Temp->Write();
			}
		}
		
		// - deflection angles 1D (also last-spill version)
		nBins = nBins_thOut; // range settings

		cond = Form("%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_inCry.c_str()); // set total boolean
		
		varName = Form("thDelta%s", direction.c_str());
		histName = Form("h1_thDelta%s", direction.c_str());
		
		if(b_outTracking&&b_h1_thDelta){
			tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, -r_thOut, r_thOut), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("%s deflection", i==0? "horizontal" : "vertical"));
			h1Temp->GetXaxis()->SetTitle("deflection angle [urad]");
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, -r_thOut, r_thOut), cond.c_str(), "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("%s deflection", i==0? "horizontal" : "vertical"));
				h1Temp->GetXaxis()->SetTitle("deflection angle [urad]");
				h1Temp->Write();
			}
		}

		// - deflection angles 1D (one per step)
		nBins = nBins_thOut; // range settings

		varName = Form("thDelta%s", direction.c_str());

		cond = Form("%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_inCry.c_str()); // set total boolean

		if(b_outTracking&&b_h1_thDelta_SingleStep&&b_plotsPerSingleStep){
			for(int j=ceil(tree->GetMinimum("iStep")); j<=ceil(tree->GetMaximum("iStep")); j++){
				cond1 = Form("%s&&(iStep==%d)", cond.c_str(), j); // adding step selection to boolean

				histName = Form("h1_thDelta%s_singleSteps%d", direction.c_str(), j);			
				tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, -r_thOut, r_thOut), cond1.c_str(), "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("%s deflection, step %d", i==0? "horizontal" : "vertical", j));
				h1Temp->GetXaxis()->SetTitle("deflection angle [urad]");
				h1Temp->Write();
			}
		}
	}
	
	// histograms, correlations of output multiplicity (1st one)
	nBinsY = nBins_nHit; // range settings
	
	varNameY = Form("nHit[%d]", iOut[0]);
	
	// - vs. iStep	
	varNameX = Form("iStep");

	cond = Form("%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_thOutLim.c_str()); // set total boolean

	/*// if(b_outTracking&&b_h2_iStep_thDelta){
	if(b_DESY22){ // DESY22, ADD BOOLEAN!
		histName = Form("h2_iStep_nHitOut");
		tree->Draw(Form("%s:%s>>%s(,,,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsY), cond.c_str(), "goff");
cout << "ciao2-" << iOut[0] << endl;
		h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
cout << "ciao3-" << iOut[0] << endl;
		h2Temp->SetTitle(Form("vs. scan step"));
		h2Temp->GetXaxis()->SetTitle("scan step nr.");
		h2Temp->GetYaxis()->SetTitle("nr. of hits");
		h2Temp->Write();
cout << "ciao5-" << iOut[0] << endl;
		
		histName = Form("tp_iStep_nHitOut");
		tree->Draw(Form("%s:%s>>%s", varNameY.c_str(), varNameX.c_str(), histName.c_str()), cond.c_str(), "goff,profs");
		tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
		tpTemp->SetTitle(Form("output multiplicity (1st one) vs. scan step"));
		tpTemp->GetXaxis()->SetTitle("scan step nr.");
		tpTemp->GetYaxis()->SetTitle("nr. of hits");
		tpTemp->Write();
	}*/

	// - vs. goniometer DOFs
	cond = Form("%s&&%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_thOutLim.c_str()); // set total boolean

	for(int i=0; i<nGonio; i++){
		varNameX = Form("xGonio%d", i);
		
		nBinsX = nBinsX_xGonio_thDelta;
		
		/*// if(b_outTracking&&b_h2_xGonio_thDelta){
		if(b_DESY22){ // DESY22, ADD BOOLEAN!
			histName = Form("h2_xGonio%d_nHitOut", i);
			tree->Draw(Form("%s:%s>>%s(%d,,,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, nBinsY), cond.c_str(), "goff");
			h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
			h2Temp->SetTitle(Form("output multiplicity (1st one) vs. goniometer DOF %d", i));
			h2Temp->GetXaxis()->SetTitle("goniometer DOF");
			h2Temp->GetYaxis()->SetTitle("nr. of hits");
			h2Temp->Write();
		
			histName = Form("tp_xGonio%d_nHitOut", i);
			tree->Draw(Form("%s:%s>>%s(%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX), cond.c_str(), "goff,profs");
			tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
			tpTemp->SetTitle(Form("output multiplicity (1st one) vs. goniometer DOF %d", i));
			tpTemp->GetXaxis()->SetTitle("goniometer DOF");
			tpTemp->GetYaxis()->SetTitle("nr. of hits");
			tpTemp->Write();
		}*/
	}

	// histograms, correlations of deflections
	nBinsY = nBins_thOut; // range settings

	bool b_xCryX_thDelta, b_xCryY_thDelta, b_thInX_thDelta, b_thInY_thDelta;

	for(int j=0; j<2; j++){		
		string direction = j%2==0 ? Form("X") : Form("Y");
		string directionName = j%2==0 ? Form("horizontal") : Form("vertical");

		varNameY = Form("thDelta%s", direction.c_str());
		
		// - vs. iStep	
		varNameX = Form("iStep");

		cond = Form("%s&&%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_thOutLim.c_str()); // set total boolean
		
		if(b_outTracking&&b_h2_iStep_thDelta){
			histName = Form("h2_iStep_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(,,,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsY), cond.c_str(), "goff");
			h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
			h2Temp->SetTitle(Form("%s deflection vs. scan step", directionName.c_str()));
			h2Temp->GetXaxis()->SetTitle("scan step nr.");
			h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
			h2Temp->Write();
			
			histName = Form("tp_iStep_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s", varNameY.c_str(), varNameX.c_str(), histName.c_str()), cond.c_str(), "goff,profs");
			tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
			tpTemp->SetTitle(Form("%s deflection vs. scan step", directionName.c_str()));
			tpTemp->GetXaxis()->SetTitle("scan step nr.");
			tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
			tpTemp->Write();
		}
		
		// - vs. goniometer DOFs
		cond = Form("%s&&%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_thOutLim.c_str()); // set total boolean

		for(int i=0; i<nGonio; i++){
			varNameX = Form("xGonio%d", i);
			
			nBinsX = nBinsX_xGonio_thDelta;
			
			if(b_outTracking&&b_h2_xGonio_thDelta){
				histName = Form("h2_xGonio%d_thDelta%s", i, direction.c_str());
				tree->Draw(Form("%s:%s>>%s(%d,,,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, nBinsY), cond.c_str(), "goff");
				h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
				h2Temp->SetTitle(Form("%s deflection vs. goniometer DOF %d", directionName.c_str(), i));
				h2Temp->GetXaxis()->SetTitle("goniometer DOF");
				h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
				h2Temp->Write();
			
				histName = Form("tp_xGonio%d_thDelta%s", i, direction.c_str());
				tree->Draw(Form("%s:%s>>%s(%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX), cond.c_str(), "goff,profs");
				tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
				tpTemp->SetTitle(Form("%s deflection vs. goniometer DOF %d", directionName.c_str(), i));
				tpTemp->GetXaxis()->SetTitle("goniometer DOF");
				tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
				tpTemp->Write();
			}
		}

		// - vs. xCryX (total & one per step)
		cond = Form("%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_thOutLim.c_str()); // set total boolean

		varNameX = Form("xCryX");

		nBinsX = nBinsX_xCry_thDelta;

		if(b_outTracking&&b_h2_xCryX_thDelta){
			histName = Form("h2_xCryX_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(%d,%f,%f,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, r_xCryMin, r_xCryMax, nBinsY), cond.c_str(), "goff");
			h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
			h2Temp->SetTitle(Form("%s deflection vs. horizontal position at crystal", directionName.c_str()));
			h2Temp->GetXaxis()->SetTitle("position [cm]");
			h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
			h2Temp->Write();
			
			histName = Form("tp_xCryX_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(%d,%f,%f)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, r_xCryMin, r_xCryMax), cond.c_str(), "goff,profs");
			tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
			tpTemp->SetTitle(Form("%s deflection vs. horizontal position at crystal", directionName.c_str()));
			tpTemp->GetXaxis()->SetTitle("position [cm]");
			tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
			tpTemp->Write();
		}

		for(int i=ceil(tree->GetMinimum("iStep")); i<=ceil(tree->GetMaximum("iStep")); i++){			
			cond1 = Form("%s&&(iStep==%d)", cond.c_str(), i); // adding step selection to boolean
		
			if(b_h2_xCryX_thDelta_SingleStep&&b_plotsPerSingleStep){
				histName = Form("h2_xCryX_thDelta%s_singleSteps%d", direction.c_str(), i);
				tree->Draw(Form("%s:%s>>%s(%d,%f,%f,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, r_xCryMin, r_xCryMax, nBinsY), cond.c_str(), "goff");
				h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
				h2Temp->SetTitle(Form("%s deflection vs. horizontal position at crystal, step %d", directionName.c_str(), i));
				h2Temp->GetXaxis()->SetTitle("position [cm]");
				h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
				h2Temp->Write();
				
				histName = Form("tp_xCryX_thDelta%s_singleSteps%d", direction.c_str(), i);
				tree->Draw(Form("%s:%s>>%s(%d,%f,%f)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, r_xCryMin, r_xCryMax), cond.c_str(), "goff,profs");
				tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
				tpTemp->SetTitle(Form("%s deflection vs. horizontal position at crystal, step %d", directionName.c_str(), i));
				tpTemp->GetXaxis()->SetTitle("position [cm]");
				tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
				tpTemp->Write();
			}
		}

		// - vs. xCryY (total & one per step)
		cond = Form("%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_thOutLim.c_str()); // set total boolean

		varNameX = Form("xCryY");

		nBinsX = nBinsX_xCry_thDelta;

		if(b_outTracking&&b_h2_xCryY_thDelta){
			histName = Form("h2_xCryY_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(%d,%f,%f,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, r_xCryMin, r_xCryMax, nBinsY), cond.c_str(), "goff");
			h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
			h2Temp->SetTitle(Form("%s deflection vs. vertical position at crystal", directionName.c_str()));
			h2Temp->GetXaxis()->SetTitle("position [cm]");
			h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
			h2Temp->Write();
			
			histName = Form("tp_xCryY_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(%d,%f,%f)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, r_xCryMin, r_xCryMax), cond.c_str(), "goff,profs");
			tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
			tpTemp->SetTitle(Form("%s deflection vs. vertical position at crystal", directionName.c_str()));
			tpTemp->GetXaxis()->SetTitle("position [cm]");
			tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
			tpTemp->Write();
		}

		for(int i=ceil(tree->GetMinimum("iStep")); i<=ceil(tree->GetMaximum("iStep")); i++){			
			cond1 = Form("%s&&(iStep==%d)", cond.c_str(), i); // adding step selection to boolean
		
			if(b_h2_xCryY_thDelta_SingleStep&&b_plotsPerSingleStep){
				histName = Form("h2_xCryY_thDelta%s_singleSteps%d", direction.c_str(), i);
				tree->Draw(Form("%s:%s>>%s(%d,%f,%f,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, r_xCryMin, r_xCryMax, nBinsY), cond.c_str(), "goff");
				h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
				h2Temp->SetTitle(Form("%s deflection vs. vertical position at crystal, step %d", directionName.c_str(), i));
				h2Temp->GetXaxis()->SetTitle("position [cm]");
				h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
				h2Temp->Write();
				
				histName = Form("tp_xCryY_thDelta%s_singleSteps%d", direction.c_str(), i);
				tree->Draw(Form("%s:%s>>%s(%d,%f,%f)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, r_xCryMin, r_xCryMax), cond.c_str(), "goff,profs");
				tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
				tpTemp->SetTitle(Form("%s deflection vs. vertical position at crystal, step %d", directionName.c_str(), i));
				tpTemp->GetXaxis()->SetTitle("position [cm]");
				tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
				tpTemp->Write();
			}
		}
		
		// - vs. thInX (one per step)
		cond = Form("%s&&%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_thOutLim.c_str()); // set total boolean

		varNameX = Form("thInX");

		nBinsX = nBinsX_thIn_thDelta;

		if(b_outTracking&&b_h2_thInX_thDelta){
			histName = Form("h2_thInX_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(%d,,,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, nBinsY), cond1.c_str(), "goff");
			h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
			h2Temp->SetTitle(Form("%s deflection vs. horizontal input angle", directionName.c_str()));
			h2Temp->GetXaxis()->SetTitle("input angle [urad]");
			h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
			h2Temp->Write();
			
			histName = Form("tp_thInX_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX), cond.c_str(), "goff,profs");
			tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
			tpTemp->SetTitle(Form("%s deflection vs. horizontal input angle", directionName.c_str()));
			tpTemp->GetXaxis()->SetTitle("input angle [urad]");
			tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
			tpTemp->Write();
		}

		for(int i=ceil(tree->GetMinimum("iStep")); i<=ceil(tree->GetMaximum("iStep")); i++){			
			cond1 = Form("%s&&(iStep==%d)", cond.c_str(), i); // adding step selection to boolean
		
			if(b_h2_thInX_thDelta_SingleStep&&b_plotsPerSingleStep){
				histName = Form("h2_thInX_thDelta%s_singleSteps%d", direction.c_str(), i);
				tree->Draw(Form("%s:%s>>%s(%d,,,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, nBinsY), cond1.c_str(), "goff");
				h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
				h2Temp->SetTitle(Form("%s deflection vs. horizontal input angle, step %d", directionName.c_str(), i));
				h2Temp->GetXaxis()->SetTitle("input angle [urad]");
				h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
				h2Temp->Write();
				
				histName = Form("tp_thInX_thDelta%s_singleSteps%d", direction.c_str(), i);
				tree->Draw(Form("%s:%s>>%s(%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX), cond.c_str(), "goff,profs");
				tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
				tpTemp->SetTitle(Form("%s deflection vs. horizontal input angle, step %d", directionName.c_str(), i));
				tpTemp->GetXaxis()->SetTitle("input angle [urad]");
				tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
				tpTemp->Write();
			}
		}

		// - vs. thInY (one per step)
		cond = Form("%s&&%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_singleHitOut.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_thOutLim.c_str()); // set total boolean

		varNameX = Form("thInY");

		nBinsX = nBinsX_thIn_thDelta;

		if(b_outTracking&&b_h2_thInY_thDelta){
			histName = Form("h2_thInY_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(%d,,,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, nBinsY), cond1.c_str(), "goff");
			h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
			h2Temp->SetTitle(Form("%s deflection vs. vertical input angle", directionName.c_str()));
			h2Temp->GetXaxis()->SetTitle("input angle [urad]");
			h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
			h2Temp->Write();
			
			histName = Form("tp_thInY_thDelta%s", direction.c_str());
			tree->Draw(Form("%s:%s>>%s(%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX), cond.c_str(), "goff,profs");
			tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
			tpTemp->SetTitle(Form("%s deflection vs. vertical input angle", directionName.c_str()));
			tpTemp->GetXaxis()->SetTitle("input angle [urad]");
			tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
			tpTemp->Write(); 
		}

		for(int i=ceil(tree->GetMinimum("iStep")); i<=ceil(tree->GetMaximum("iStep")); i++){			
			cond1 = Form("%s&&(iStep==%d)", cond.c_str(), i); // adding step selection to boolean
		
			if(b_h2_thInY_thDelta_SingleStep&&b_plotsPerSingleStep){
				histName = Form("h2_thInY_thDelta%s_singleSteps%d", direction.c_str(), i);
				tree->Draw(Form("%s:%s>>%s(%d,,,%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, nBinsY), cond1.c_str(), "goff");
				h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
				h2Temp->SetTitle(Form("%s deflection vs. vertical input angle, step %d", directionName.c_str(), i));
				h2Temp->GetXaxis()->SetTitle("input angle [urad]");
				h2Temp->GetYaxis()->SetTitle("deflection angle [urad]");
				h2Temp->Write();
				
				histName = Form("tp_thInY_thDelta%s_singleSteps%d", direction.c_str(), i);
				tree->Draw(Form("%s:%s>>%s(%d)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX), cond.c_str(), "goff,profs");
				tpTemp = (TProfile*)gDirectory->Get(histName.c_str());
				tpTemp->SetTitle(Form("%s deflection vs. vertical input angle, step %d", directionName.c_str(), i));
				tpTemp->GetXaxis()->SetTitle("input angle [urad]");
				tpTemp->GetYaxis()->SetTitle("deflection angle [urad]");
				tpTemp->Write();
			}
		}
	}
	
	// histograms, digitizer data
	for(int i=0; i<nDigi; i++){
		
		// - times, with no cuts (also last-spill version)
		nBins = nBins_time;

		varName = Form("digiTime[%d]", i);
		histName = Form("h1_time%d_noCuts", i);
		
		if(b_h1_time_noCuts){
			tree->Draw(Form("%s>>%s(%d, 0.0, %f)", varName.c_str(), histName.c_str(), nBins, r_timeMax), "", "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("digi. ch. %d signal (no cuts)", i));
			h1Temp->GetXaxis()->SetTitle("time [sample nr.]");
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d, 0.0, %f)", varName.c_str(), histName.c_str(), nBins, r_timeMax), "", "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("digi. ch. %d signal (no cuts)", i));
				h1Temp->GetXaxis()->SetTitle("time [sample nr.]");
				h1Temp->Write();			
			}
		}
		
		// - PHs, with no cuts (also last-spill version)
		nBins = nBins_ph;
		
		varName = Form("digiPH[%d]", i);
		histName = Form("h1_ph%d_noCuts", i);
		
		if(b_h1_ph_noCuts){
			tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), "", "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("digi. ch. %d signal (no cuts)", i));
			h1Temp->GetXaxis()->SetTitle("PH [ADC]");
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), "", "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("digi. ch. %d signal (no cuts)", i));
				h1Temp->GetXaxis()->SetTitle("PH [ADC]");
				h1Temp->Write();			
			}
		}
		
		// - times vs. PHs, with no cuts
		nBinsX = nBins_time ; nBinsY = nBins_ph ;
		xRangeLeft = 0. ; xRangeRight = r_timeMax ;
		yRangeLeft = 0. ; yRangeRight = r_phMax ;
		
		histName = Form("h2_timeVsPH%d_noCuts", i);
		varNameX = Form("digiTime[%d]", i);
		varNameY = Form("digiPH[%d]", i);
		
		if(b_h2_time_ph_noCuts){
			tree->Draw(Form("%s:%s>>%s(%d,%f,%f,%d,%f,%f)", varNameY.c_str(), varNameX.c_str(), histName.c_str(), nBinsX, xRangeLeft, xRangeRight, nBinsY, yRangeLeft, yRangeRight), "", "goff");
			h2Temp = (TH2F*)gDirectory->Get(histName.c_str());
			h2Temp->SetTitle(Form("digi. ch. %d signal (no cuts)", i));
			h2Temp->GetXaxis()->SetTitle("time [sample nr.]");
			h2Temp->GetYaxis()->SetTitle("PH [ADC]");
			h2Temp->Write();
		}
		
		// - PHs, equalised, with no cuts (also last-spill version)
		nBins = nBins_ph;
		
		varName = Form("digiPHEq%d", i);
		histName = Form("h1_phEq%d_noCuts", i);
		
		if(b_h1_ph_noCuts){
			tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), "", "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("digi. ch. %d equalised signal (no cuts)", i));
			h1Temp->GetXaxis()->SetTitle("PH [ADC]");
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), "", "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("digi. ch. %d equalised signal (no cuts)", i));
				h1Temp->GetXaxis()->SetTitle("PH [ADC]");
				h1Temp->Write();			
			}
		}
		
		// - PHs, equalised, with physics cuts (also last-spill version)
		cond = Form("%s&&%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_digiTime[i].c_str(), b_digiPHEq[i].c_str()); // set total boolean
		
		nBins = nBins_ph;
		
		varName = Form("digiPHEq%d", i);
		histName = Form("h1_phEq%d", i);
		
		if(b_h1_ph_noCuts){
			tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("digi. ch. %d equalised signal", i));
			h1Temp->GetXaxis()->SetTitle("PH [ADC]");
			h1Temp->Write();
			
			if(b_plotLastSpill&&(nEntriesLast!=0)){
				histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
				treeLast->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), cond.c_str(), "goff");
				h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
				h1Temp->SetTitle(Form("digi. ch. %d equalised signal", i));
				h1Temp->GetXaxis()->SetTitle("PH [ADC]");
				h1Temp->Write();			
			}
		}
	}
	
	// histograms, total calorimeters
	// - forward calorimeter total PH, with physics cuts (also last-spill version)
	if(b_fwdCalo&&b_h1_phTot_fwdCalo){		
		cond = Form("%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_digiTimeCaloFwd.c_str()); // set total boolean
		
		nBins = nBins_ph;
		
		varName = Form("PHCaloFwd");
		histName = Form("h1_phTot_fwdCalo");

		tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), cond.c_str(), "goff");
		h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
		h1Temp->SetTitle(Form("forward calorimeter total signal"));
		h1Temp->GetXaxis()->SetTitle("PH [ADC]");
		h1Temp->Write();
		
		if(b_plotLastSpill&&(nEntriesLast!=0)){
			histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
			treeLast->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("forward calorimeter total signal"));
			h1Temp->GetXaxis()->SetTitle("PH [ADC]");
			h1Temp->Write();			
		}
	}
	
	// - lateral calorimeter total PH, with physics cuts (also last-spill version)
	if(b_latCalo&&b_h1_phTot_latCalo){
		cond = Form("%s&&%s&&%s&&%s&&%s", b_xRawBaseTrackErrors.c_str(), b_singleHitIn.c_str(), b_inAligned.c_str(), b_inCry.c_str(), b_digiTimeCaloLat.c_str()); // set total boolean
		
		nBins = nBins_ph;
		
		varName = Form("PHCaloLat");
		histName = Form("h1_phTot_latCalo");
	
		tree->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), cond.c_str(), "goff");
		h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
		h1Temp->SetTitle(Form("lateral calorimeter total signal"));
		h1Temp->GetXaxis()->SetTitle("PH [ADC]");
		h1Temp->Write();
		
		if(b_plotLastSpill&&(nEntriesLast!=0)){
			histName = Form("%s%s", histName.c_str(), histNameAdd.c_str());
			treeLast->Draw(Form("%s>>%s(%d, %f, %f)", varName.c_str(), histName.c_str(), nBins, 0., r_phMax), cond.c_str(), "goff");
			h1Temp = (TH1F*)gDirectory->Get(histName.c_str());
			h1Temp->SetTitle(Form("lateral calorimeter total signal"));
			h1Temp->GetXaxis()->SetTitle("PH [ADC]");
			h1Temp->Write();			
		}
	}
	
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	// custom histograms here below (remember to add canvases below also) //vvvvvvvvvvvvvvvvvvvvvvv
	
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	cout<<"hisograms generated, now generating canvases..."<<endl;
	
// final canvases (& other text info) creation and saving

	TCanvas *canvas;
	TH1F *h1TempOut[400];
	TH2F *h2TempOut[400];
	TProfile *tpTempOut[400];
	TGraph *grTempOut[400];
	string canvasName, canvasTitle;
	
	gStyle->SetOptFit(optFit);
	gStyle->SetOptStat(optStat1D);

	int jMin=ceil(tree->GetMinimum("iStep")) , jMax=ceil(tree->GetMaximum("iStep")) ;
	
	string outPdfPath= Form("outpdf/");

	string outFigsFileName = Form("outfigs/figsFile.root");
	TFile *outFigsFile = new TFile(outFigsFileName.c_str(), "RECREATE");

	// canvases w/ 1D plots
	
	// canvas, multiplicities
	if(b_h1_nHit){
		canvasName = Form("c_nHit");
		canvas = new TCanvas(canvasName.c_str(), "Si layer multiplicities");
		canvas->Divide(2, ceil(nSi/2));

		for(int i=0; i<nSi; i++){
			canvas->cd(i+1);
			gPad->SetLogy(); // log scale in y?

			h1TempOut[i] = (TH1F*)outHistFile->Get(Form("h1_nHit%d", i));
			h1TempOut[i]->SetLineColor(cLineTotal);
			h1TempOut[i]->SetFillColor(cFillTotal);
			h1TempOut[i]->Draw("HIST");

			if(b_plotLastSpill&&(nEntriesLast!=0)){
				h1TempOut[nSi+i] = (TH1F*)outHistFile->Get(Form("h1_nHit%d_last", i));
				h1TempOut[nSi+i]->Scale(h1TempOut[i]->Integral() / h1TempOut[nSi+i]->Integral());
				h1TempOut[nSi+i]->SetLineColor(cLineLast);
				h1TempOut[nSi+i]->Draw("HIST SAME");
			}

		}
		canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
		canvas->Write();
	}
	
	// canvas, beam profiles
	if(b_h1_xRaw){
		canvasName = Form("c_xRaw");
		canvas = new TCanvas(canvasName.c_str(), "Si layer beam profiles");
		canvas->Divide(2, ceil(nSi/2));

		for(int i=0; i<nSi; i++){
			canvas->cd(i+1);
			gPad->SetLogy(); // log scale in y?

			h1TempOut[i] = (TH1F*)outHistFile->Get(Form("h1_xRaw%d", i));
			h1TempOut[i]->SetLineColor(cLineTotal);
			h1TempOut[i]->SetFillColor(cFillTotal);
			h1TempOut[i]->Draw("HIST");

			if(b_plotLastSpill&&(nEntriesLast!=0)){
				h1TempOut[nSi+i] = (TH1F*)outHistFile->Get(Form("h1_xRaw%d_last", i));
				h1TempOut[nSi+i]->Scale(h1TempOut[i]->Integral() / h1TempOut[nSi+i]->Integral());
				h1TempOut[nSi+i]->SetLineColor(cLineLast);
				h1TempOut[nSi+i]->Draw("HIST SAME");
			}

		}
		canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
		canvas->Write();
	}
	
	// canvas, angles
	if(b_h1_thIn||b_h1_thOut||b_h1_thDelta){
		canvasName = Form("c_th");
		canvas = new TCanvas(canvasName.c_str(), "angles");
		canvas->Divide(2, 3);

		gStyle->SetOptStat(optStat1D);

		int k=0;
		for(int i=0; i<6; i++){
			canvas->cd(k+1);
			gPad->SetLogy(); // log scale in y?

			string direction = i%2==0 ? Form("X") : Form("Y");
			string stage;
			bool b_local;
			if(i<2){
				stage = Form("In");
				b_local = b_h1_thIn;
			}
			else if((i>=2)&&(i<4)){
				stage = Form("Out");
				b_local = b_h1_thOut&&b_outTracking;
			}
			else {
				stage = Form("Delta");
				b_local = b_h1_thDelta&&b_outTracking;
			}

			if(b_local){
				h1TempOut[i] = (TH1F*)outHistFile->Get(Form("h1_th%s%s", stage.c_str(), direction.c_str()));
				h1TempOut[i]->SetLineColor(cLineTotal);
				h1TempOut[i]->SetFillColor(cFillTotal);
				h1TempOut[i]->Draw("HIST");

				if(b_plotLastSpill&&(nEntriesLast!=0)){
					h1TempOut[2+i] = (TH1F*)outHistFile->Get(Form("h1_th%s%s_last", stage.c_str(), direction.c_str()));
					h1TempOut[2+i]->Scale(h1TempOut[i]->Integral() / h1TempOut[2+i]->Integral());
					h1TempOut[2+i]->SetLineColor(cLineLast);
					h1TempOut[2+i]->Draw("HIST SAME");
				}
				
				k+=1;
			}

		}
		canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
		canvas->Write();
	}

	// canvas(es), scan single-step deflections
	if(b_outTracking&&b_h1_thDelta_SingleStep&&b_plotsPerSingleStep){
		for(int jCanv=0; jCanv<=ceil((jMax-jMin+1)/6); jCanv++){
			canvasName = Form("c_thDelta_singleSteps%d", jCanv);
			canvasTitle = Form("single-step deflections, subcanvas %d", jCanv);
			canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
			canvas->Divide(2, 6);
			for(int j=0; j<6; j++){
				for(int k=0; k<2; k++){
					canvas->cd(2*j+k+1);
					gPad->SetLogy(); // log scale in y?

					string direction = k%2==0 ? Form("X") : Form("Y");

					if(jMin+6*jCanv+j<=jMax){
						h1TempOut[2*j+k%2] = (TH1F*)outHistFile->Get(Form("h1_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						h1TempOut[2*j+k%2]->SetLineColor(cLineTotal);
						h1TempOut[2*j+k%2]->SetFillColor(cFillTotal);
						h1TempOut[2*j+k%2]->Draw("HIST");

						if(b_plotLastSpill&&(nEntriesLast!=0)){
							h1TempOut[13+2*j+k%2] = (TH1F*)outHistFile->Get(Form("h1_thDelta%s_last", direction.c_str()));
							h1TempOut[13+2*j+k%2]->Scale(h1TempOut[2*j+k%2]->Integral() / h1TempOut[13+2*j+k%2]->Integral());
							h1TempOut[13+2*j+k%2]->SetLineColor(cLineLast);
							h1TempOut[13+2*j+k%2]->Draw("HIST SAME");
						}
					}

				}
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}

	// canvases w/ 2D plots
	gStyle->SetOptStat(optStat2D);
	
	/*// canvas, goniometer correlations of output multiplicity (1st one)
	// if(b_outTracking&&(b_h2_iStep_thDelta||b_h2_xGonio_thDelta)){
	if(b_DESY22){ // DESY22, ADD BOOLEAN!
		canvasName = Form("c_xGonio_nHitOut");
		canvasTitle = Form("goniometer correlations of output multiplicity (1st one)");
		canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
		canvas->Divide(2, nGonio+1);

		int k=0;
		for(int i=0; i<nGonio+1; i++){
			string dof;
			string grVarName;
			bool b_local;
			if(i==0){
				dof = Form("iStep");
				grVarName = Form("scan step nr.");
				// b_local = b_h2_iStep_thDelta;
				b_local = b_DESY22; // DESY22, ADD BOOLEAN!
			}
			else {
				dof = Form("xGonio%d", i-1);
				grVarName = Form("goniometer DOF");
				// b_local = b_h2_xGonio_thDelta;
				b_local = b_DESY22; // DESY22, ADD BOOLEAN!
			}

			if(b_local){
				h2TempOut[i] = (TH2F*)outHistFile->Get(Form("h2_%s_nHitOut", dof.c_str()));
				tpTempOut[i] = (TProfile*)outHistFile->Get(Form("tp_%s_nHitOut", dof.c_str()));
				tpTempOut[i]->SetLineColor(2);

				canvas->cd(2*k+1);
				h2TempOut[i]->Draw("COL");
				tpTempOut[i]->Draw("L SAME");
				
				canvas->cd(2*k+2);
				double xErr[tpTempOut[i]->GetNbinsX()], yErr[tpTempOut[i]->GetNbinsX()];
				for(int j=0; j<tpTempOut[i]->GetNbinsX(); j++){
					xErr[j] = h2TempOut[i]->GetXaxis()->GetBinCenter(j);
					yErr[j] = abs(tpTempOut[i]->GetBinError(j));
				}
				grTempOut[i] = new TGraph(tpTempOut[i]->GetNbinsX(), xErr, yErr);
				grTempOut[i]->SetTitle("errors associated to TProfile (to the left)");
				grTempOut[i]->GetXaxis()->SetTitle(grVarName.c_str());
				grTempOut[i]->GetYaxis()->SetTitle("error on output multiplicity from TProfile [urad]");
				grTempOut[i]->SetLineColor(2);
				grTempOut[i]->SetMarkerColor(2);
				grTempOut[i]->Draw("A P L");
				
				k+=1;
			}
		}
		canvas->Draw();
		canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
		canvas->Write();
	}*/

	// canvas, goniometer correlations of deflections
	if(b_outTracking&&(b_h2_iStep_thDelta||b_h2_xGonio_thDelta)){
		for(int i=0; i<2; i++){
			string direction = i%2==0 ? Form("X") : Form("Y");
			string directionName = i%2==0 ? Form("horizontal") : Form("vertical");
			canvasName = Form("c_xGonio_thDelta%s", direction.c_str());
			canvasTitle = Form("goniometer correlations of %s deflection", directionName.c_str());
			canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
			canvas->Divide(2, nGonio+1);

			int k=0;
			for(int i=0; i<nGonio+1; i++){
				string dof;
				string grVarName;
				bool b_local;
				if(i==0){
					dof = Form("iStep");
					grVarName = Form("scan step nr.");
					b_local = b_h2_iStep_thDelta;
				}
				else {
					dof = Form("xGonio%d", i-1);
					grVarName = Form("goniometer DOF");
					b_local = b_h2_xGonio_thDelta;
				}

				if(b_local){
					h2TempOut[i] = (TH2F*)outHistFile->Get(Form("h2_%s_thDelta%s", dof.c_str(), direction.c_str()));
					tpTempOut[i] = (TProfile*)outHistFile->Get(Form("tp_%s_thDelta%s", dof.c_str(), direction.c_str()));
					tpTempOut[i]->SetLineColor(2);

					canvas->cd(2*k+1);
					h2TempOut[i]->Draw("COL");
					tpTempOut[i]->Draw("L SAME");
					
					canvas->cd(2*k+2);
					double xErr[tpTempOut[i]->GetNbinsX()], yErr[tpTempOut[i]->GetNbinsX()];
					for(int j=0; j<tpTempOut[i]->GetNbinsX(); j++){
						xErr[j] = h2TempOut[i]->GetXaxis()->GetBinCenter(j);
						yErr[j] = abs(tpTempOut[i]->GetBinError(j));
					}
					grTempOut[i] = new TGraph(tpTempOut[i]->GetNbinsX(), xErr, yErr);
					grTempOut[i]->SetTitle("errors associated to TProfile (to the left)");
					grTempOut[i]->GetXaxis()->SetTitle(grVarName.c_str());
					grTempOut[i]->GetYaxis()->SetTitle("error on deflection angle from TProfile [urad]");
					grTempOut[i]->SetLineColor(2);
					grTempOut[i]->SetMarkerColor(2);
					grTempOut[i]->Draw("A P L");
					
					k+=1;
				}
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}

	// canvas, scan correlations between deflections and input beam features (xCryX/Y, thInX/Y)
	if(b_outTracking&&(b_h2_xCryX_thDelta||b_h2_xCryY_thDelta||b_h2_thInX_thDelta||b_h2_thInY_thDelta)){
		for(int l=0; l<2; l++){
			string direction = l%2==0 ? Form("X") : Form("Y");
			string directionName = l%2==0 ? Form("horizontal") : Form("vertical");
			canvasName = Form("c_inputBeam_thDelta%s", direction.c_str());
			canvasTitle = Form("input beam correlations of %s deflection", directionName.c_str());
			canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
			canvas->Divide(2, 4);

			int k=0;
			for(int i=0; i<4; i++){
				string varNameX;
				string VarNameY;
				string grVarName;
				bool b_local;
				if(i==0){
					varNameX = Form("xCryX");
					varNameY = Form("thDelta%s", direction.c_str());
					b_local = b_h2_xCryX_thDelta;
					grVarName = Form("position [cm]");
				}
				else if(i==1){
					varNameX = Form("xCryY");
					varNameY = Form("thDelta%s", direction.c_str());
					b_local = b_h2_xCryY_thDelta;
					grVarName = Form("position [cm]");
				}
				else if(i==2){
					varNameX = Form("thInX");
					varNameY = Form("thDelta%s", direction.c_str());
					b_local = b_h2_thInX_thDelta;
					grVarName = Form("input angle [urad]");
				}
				else{
					varNameX = Form("thInY");
					varNameY = Form("thDelta%s", direction.c_str());
					b_local = b_h2_thInY_thDelta;
					grVarName = Form("input angle [urad]");
				}

				if(b_local){
					h2TempOut[i] = (TH2F*)outHistFile->Get(Form("h2_%s_%s", varNameX.c_str(), varNameY.c_str()));
					tpTempOut[i] = (TProfile*)outHistFile->Get(Form("tp_%s_%s", varNameX.c_str(), varNameY.c_str()));
					tpTempOut[i]->SetLineColor(2);

					canvas->cd(2*k+1);
					h2TempOut[i]->Draw("COL");
					tpTempOut[i]->Draw("L SAME");
					
					canvas->cd(2*k+2);
					double xErr[tpTempOut[i]->GetNbinsX()], yErr[tpTempOut[i]->GetNbinsX()];
					for(int j=0; j<tpTempOut[i]->GetNbinsX(); j++){
						xErr[j] = h2TempOut[i]->GetXaxis()->GetBinCenter(j);
						yErr[j] = abs(tpTempOut[i]->GetBinError(j));
					}
					grTempOut[i] = new TGraph(tpTempOut[i]->GetNbinsX(), xErr, yErr);
					grTempOut[i]->SetTitle("errors associated to TProfile (to the left)");
					grTempOut[i]->GetXaxis()->SetTitle(grVarName.c_str());
					grTempOut[i]->GetYaxis()->SetTitle("error on deflection angle from TProfile [urad]");
					grTempOut[i]->SetLineColor(2);
					grTempOut[i]->SetMarkerColor(2);
					grTempOut[i]->Draw("A P L");

					k+=1;
				}
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}

	// canvas(es), scan single-step correlations between deflections and horizontal beam profile (xCryX)
	for(int k=0; k<2; k++){
		if(b_outTracking&&b_h2_xCryX_thDelta_SingleStep&&b_plotsPerSingleStep){
			for(int jCanv=0; jCanv<=ceil((jMax-jMin+1)/6); jCanv++){
				string direction = k%2==0 ? Form("X") : Form("Y");
				string directionName = k%2==0 ? Form("horizontal") : Form("vertical");
				canvasName = Form("c_xCryX_thDelta%s_singleSteps%d", direction.c_str(), jCanv);
				canvasTitle = Form("single-step horizontal beam profile correlations of %s deflection, canvas %d", directionName.c_str(), jCanv);
				canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
				canvas->Divide(2, 3);
				for(int j=0; j<6; j++){
					canvas->cd(j+1);

					if(jMin+6*jCanv+j<=jMax){
						h2TempOut[j] = (TH2F*)outHistFile->Get(Form("h2_xCryX_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						tpTempOut[j] = (TProfile*)outHistFile->Get(Form("tp_xCryX_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						tpTempOut[j]->SetLineColor(2);

						h2TempOut[j]->Draw("COL");
						tpTempOut[j]->Draw("SAME");
					}
				}
				canvas->Draw();
				canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
				canvas->Write();
			}
		}
	}

	// canvas(es), scan single-step correlations between deflections and vertical beam profile (xCryY)
	for(int k=0; k<2; k++){
		if(b_outTracking&&b_h2_xCryY_thDelta_SingleStep&&b_plotsPerSingleStep){
			for(int jCanv=0; jCanv<=ceil((jMax-jMin+1)/6); jCanv++){
				string direction = k%2==0 ? Form("X") : Form("Y");
				string directionName = k%2==0 ? Form("horizontal") : Form("vertical");
				canvasName = Form("c_xCryY_thDelta%s_singleSteps%d", direction.c_str(), jCanv);
				canvasTitle = Form("single-step vertical beam profile correlations of %s deflection, canvas %d", directionName.c_str(), jCanv);
				canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
				canvas->Divide(2, 3);
				for(int j=0; j<6; j++){
					canvas->cd(j+1);

					if(jMin+6*jCanv+j<=jMax){
						h2TempOut[j] = (TH2F*)outHistFile->Get(Form("h2_xCryY_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						tpTempOut[j] = (TProfile*)outHistFile->Get(Form("tp_xCryY_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						tpTempOut[j]->SetLineColor(2);

						h2TempOut[j]->Draw("COL");
						tpTempOut[j]->Draw("SAME");
					}
				}
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}

	// canvas(es), scan single-step correlations between deflections and horizontal input angle (thInX)
	for(int k=0; k<2; k++){
		if(b_outTracking&&b_h2_thInX_thDelta_SingleStep&&b_plotsPerSingleStep){
			for(int jCanv=0; jCanv<=ceil((jMax-jMin+1)/6); jCanv++){
				string direction = k%2==0 ? Form("X") : Form("Y");
				string directionName = k%2==0 ? Form("horizontal") : Form("vertical");
				canvasName = Form("c_thInX_thDelta%s_singleSteps%d", direction.c_str(), jCanv);
				canvasTitle = Form("single-step horizontal input angle correlations of %s deflection, canvas %d", directionName.c_str(), jCanv);
				canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
				canvas->Divide(2, 3);
				for(int j=0; j<6; j++){
					canvas->cd(j+1);

					if(jMin+6*jCanv+j<=jMax){
						h2TempOut[j] = (TH2F*)outHistFile->Get(Form("h2_thInX_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						tpTempOut[j] = (TProfile*)outHistFile->Get(Form("tp_thInX_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						tpTempOut[j]->SetLineColor(2);

						h2TempOut[j]->Draw("COL");
						tpTempOut[j]->Draw("SAME");
					}
				}
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}

	// canvas(es), scan single-step correlations between deflections and vertical input angle (thInY)
	for(int k=0; k<2; k++){
		if(b_outTracking&&b_h2_thInY_thDelta_SingleStep&&b_plotsPerSingleStep){
			for(int jCanv=0; jCanv<=ceil((jMax-jMin+1)/6); jCanv++){
				string direction = k%2==0 ? Form("X") : Form("Y");
				string directionName = k%2==0 ? Form("horizontal") : Form("vertical");
				canvasName = Form("c_thInY_thDelta%s_singleSteps%d", direction.c_str(), jCanv);
				canvasTitle = Form("single-step vertical input angle correlations of %s deflection, canvas %d", directionName.c_str(), jCanv);
				canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
				canvas->Divide(2, 3);
				for(int j=0; j<6; j++){
					canvas->cd(j+1);

					if(jMin+6*jCanv+j<=jMax){
						h2TempOut[j] = (TH2F*)outHistFile->Get(Form("h2_thInY_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						tpTempOut[j] = (TProfile*)outHistFile->Get(Form("tp_thInY_thDelta%s_singleSteps%d", direction.c_str(), jMin+6*jCanv+j));
						tpTempOut[j]->SetLineColor(2);

						h2TempOut[j]->Draw("COL");
						tpTempOut[j]->Draw("SAME");
					}
				}
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}
	
	// canvas(es), digitizer times with no cuts
	if(b_h1_time_noCuts){
		for(int jCanv=0; jCanv<ceil(nDigi/8); jCanv++){
			canvasName = Form("c_digiTime_noCuts%d", jCanv);
			canvasTitle = Form("digitizer times, no cuts, subcanvas %d", jCanv);
			canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
			canvas->Divide(2, 4);
			for(int j=0; j<8; j++){
				canvas->cd(j+1);
				// gPad->SetLogy(); // log scale in y?

				if(jCanv*8+j<nDigi){
					h1TempOut[j] = (TH1F*)outHistFile->Get(Form("h1_time%d_noCuts", 8*jCanv+j));
					h1TempOut[j]->SetLineColor(cLineTotal);
					h1TempOut[j]->SetFillColor(cFillTotal);
					h1TempOut[j]->Draw("HIST");

					if(b_plotLastSpill&&(nEntriesLast!=0)){
						h1TempOut[8+j] = (TH1F*)outHistFile->Get(Form("h1_time%d_noCuts_last", 8*jCanv+j));
						h1TempOut[8+j]->Scale(h1TempOut[j]->Integral() / h1TempOut[8+j]->Integral());
						h1TempOut[8+j]->SetLineColor(cLineLast);
						h1TempOut[8+j]->Draw("HIST SAME");
					}
				}
					
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}
	
	// canvas(es), digitizer PHs with no cuts
	if(b_h1_ph_noCuts){
		for(int jCanv=0; jCanv<ceil(nDigi/8); jCanv++){
			canvasName = Form("c_digiPh_noCuts%d", jCanv);
			canvasTitle = Form("digitizer PHs, no cuts, subcanvas %d", jCanv);
			canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
			canvas->Divide(2, 4);
			for(int j=0; j<8; j++){
				canvas->cd(j+1);
				gPad->SetLogy(); // log scale in y?

				if(jCanv*8+j<nDigi){
					h1TempOut[j] = (TH1F*)outHistFile->Get(Form("h1_ph%d_noCuts", 8*jCanv+j));
					h1TempOut[j]->SetLineColor(cLineTotal);
					h1TempOut[j]->SetFillColor(cFillTotal);
					h1TempOut[j]->Draw("HIST");

					if(b_plotLastSpill&&(nEntriesLast!=0)){
						h1TempOut[8+j] = (TH1F*)outHistFile->Get(Form("h1_ph%d_noCuts_last", 8*jCanv+j));
						h1TempOut[8+j]->Scale(h1TempOut[j]->Integral() / h1TempOut[8+j]->Integral());
						h1TempOut[8+j]->SetLineColor(cLineLast);
						h1TempOut[8+j]->Draw("HIST SAME");
					}
				}
					
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}
	
	// canvas(es), digitizer times vs. PHs with no cuts
	if(b_h2_time_ph_noCuts){
		for(int jCanv=0; jCanv<ceil(nDigi/8); jCanv++){
			canvasName = Form("c_digiTimeVsPh_noCuts%d", jCanv);
			canvasTitle = Form("digitizer times vs PHs, no cuts, subcanvas %d", jCanv);
			canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
			canvas->Divide(2, 4);
			for(int j=0; j<8; j++){
				canvas->cd(j+1);

				if(jCanv*8+j<nDigi){
					h2TempOut[j] = (TH2F*)outHistFile->Get(Form("h2_timeVsPH%d_noCuts", 8*jCanv+j));
					h2TempOut[j]->SetLineColor(cLineTotal);
					h2TempOut[j]->SetFillColor(cFillTotal);
					h2TempOut[j]->Draw("COL");
				}
					
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}
	
	// canvas(es), digitizer equalised PHs with no cuts
	if(b_h1_phEq_noCuts){
		for(int jCanv=0; jCanv<ceil(nDigi/8); jCanv++){
			canvasName = Form("c_digiPhEq_noCuts%d", jCanv);
			canvasTitle = Form("digitizer equalised PHs, no cuts, subcanvas %d", jCanv);
			canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
			canvas->Divide(2, 4);
			for(int j=0; j<8; j++){
				canvas->cd(j+1);
				gPad->SetLogy(); // log scale in y?

				if(jCanv*8+j<nDigi){
					h1TempOut[j] = (TH1F*)outHistFile->Get(Form("h1_phEq%d_noCuts", 8*jCanv+j));
					h1TempOut[j]->SetLineColor(cLineTotal);
					h1TempOut[j]->SetFillColor(cFillTotal);
					h1TempOut[j]->Draw("HIST");

					if(b_plotLastSpill&&(nEntriesLast!=0)){
						h1TempOut[8+j] = (TH1F*)outHistFile->Get(Form("h1_phEq%d_noCuts_last", 8*jCanv+j));
						h1TempOut[8+j]->Scale(h1TempOut[j]->Integral() / h1TempOut[8+j]->Integral());
						h1TempOut[8+j]->SetLineColor(cLineLast);
						h1TempOut[8+j]->Draw("HIST SAME");
					}
				}
					
			}
			canvas->Draw();
			canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
			canvas->Write();
		}
	}
	
	// canvas, specific detector in digitizer: multiplicity counters
	if(b_counters&&b_h1_phEq_counters){
		drawSpecDigiCanvas(outHistFile, nEntries, nEntriesLast, b_plotLastSpill, outPdfPath, b_counters, "counters", "multiplicity counters", cShape_counters, ls_counters);
	}
	
	// canvas, specific detector in digitizer: forward calorimeter (single channels)
	if(b_fwdCalo&&b_h1_phEq_fwdCalo){
		drawSpecDigiCanvas(outHistFile, nEntries, nEntriesLast, b_plotLastSpill, outPdfPath, b_fwdCalo, "fwdCalo", "forward calorimeter", cShape_fwdCalo, ls_fwdCalo);
	}
	
	// canvas, specific detector in digitizer: lateral calorimeter (single channels)
	if(b_latCalo&&b_h1_phEq_latCalo){
		drawSpecDigiCanvas(outHistFile, nEntries, nEntriesLast, b_plotLastSpill, outPdfPath, b_latCalo, "latCalo", "lateral calorimeter", cShape_latCalo, ls_latCalo);
	}
	
	// canvas, specific detector in digitizer: read-out crystal
	if(b_cry&&b_h1_phEq_cry){
		drawSpecDigiCanvas(outHistFile, nEntries, nEntriesLast, b_plotLastSpill, outPdfPath, b_cry, "cry", "crystal", cShape_cry, ls_cry);
	}
	
	// canvas, forward calorimeter total signal
	if(b_fwdCalo&&b_h1_phTot_fwdCalo){
		canvasName = Form("c_phTot_fwdCalo");
		canvas = new TCanvas(canvasName.c_str(), "forward calorimeter total signal");
		canvas->Divide(1, 2);

		for(int j=0; j<2; j++){
			canvas->cd(j+1);
			if(j!=0){gPad->SetLogy();}
			h1TempOut[0] = (TH1F*)outHistFile->Get(Form("h1_phTot_fwdCalo"));
			h1TempOut[0]->SetLineColor(cLineTotal);
			h1TempOut[0]->SetFillColor(cFillTotal);
			h1TempOut[0]->Draw("HIST");

			if(b_plotLastSpill&&(nEntriesLast!=0)){
				h1TempOut[1] = (TH1F*)outHistFile->Get(Form("h1_phTot_fwdCalo_last"));
				h1TempOut[1]->Scale(h1TempOut[0]->Integral() / h1TempOut[1]->Integral());
				h1TempOut[1]->SetLineColor(cLineLast);
				h1TempOut[1]->Draw("HIST SAME");
			}
		}

		canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
		canvas->Write();
	}
	
	// canvas, lateral calorimeter total signal
	if(b_latCalo&&b_h1_phTot_latCalo){
		canvasName = Form("c_phTot_latCalo");
		canvas = new TCanvas(canvasName.c_str(), "lateral calorimeter total signal");
		canvas->Divide(1, 2);

		for(int j=0; j<2; j++){
			canvas->cd(j+1);
			if(j!=0){gPad->SetLogy();}
			h1TempOut[0] = (TH1F*)outHistFile->Get(Form("h1_phTot_latCalo"));
			h1TempOut[0]->SetLineColor(cLineTotal);
			h1TempOut[0]->SetFillColor(cFillTotal);
			h1TempOut[0]->Draw("HIST");

			if(b_plotLastSpill&&(nEntriesLast!=0)){
				h1TempOut[1] = (TH1F*)outHistFile->Get(Form("h1_phTot_latCalo_last"));
				h1TempOut[1]->Scale(h1TempOut[0]->Integral() / h1TempOut[1]->Integral());
				h1TempOut[1]->SetLineColor(cLineLast);
				h1TempOut[1]->Draw("HIST SAME");
			}
		}

		canvas->SaveAs(Form("%s%s.pdf", outPdfPath.c_str(), canvasName.c_str()));
		canvas->Write();
	}
	
	//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	// custom canvases here below (remember to add histograms above also) //vvvvvvvvvvvvvvvvvvvvvvv
	
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	
// close all output files

	if(bWriteOutTreeFile){
		cout << Form("raw tree output file written: %s", outTreeFileName.c_str()) << endl;
		outTreeFile->Close();
	}
	
	cout << Form("histogram output file written: %s", outHistFileName.c_str()) << endl;
	outHistFile->Close();
	
	cout << Form("canvas output file written: %s", outFigsFileName.c_str()) << endl;
	cout << Form("also single PDFs written in %s", outPdfPath.c_str()) << endl;
	outFigsFile->Close();

// done
	cout<<"done w/ analysis!"<<endl;
	
}

int main(){krysBTOM_ana(); return 0;}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////