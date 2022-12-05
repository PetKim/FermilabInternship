#Lots of help by Dr. A. Nikolakopoulos

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>

void flux_0_CCQE()
{
  //Reading ROOT file and setting pointer
  TFile *f = new TFile("CRPASuSAv2Hybrid_AR40_flx0_flat.root", "read");

  TTree *tree = (TTree*)f->Get("FlatTree_VARS"); //Getting the tree

  //WHY set the pointer events to a new TH1F ? should just declare the pointer
  //TH1F *events = new TH1F("event", "Event", 600, 0, 6); //here previously
  TH1F *events;
  events = (TH1F*)f->Get("FlatTree_EVT");  //Getting the histogram in ROOT file

  int entries = tree->GetEntries(); //Number of entries in tree
  Float_t e, a; //Setting variables for branches to be assigned to

//NEW: here we set a flag for CCQE events , events which GENIE categorises as 1-nucleon knockout
//If we only add these events to the historgram, this is a subset of the total event distribution, so the rest of code remains the same.
  bool CCQEflag; //New flag

  //Accessing necessary branches and assigning to variables
  tree->SetBranchAddress("ELep", &e);
  tree->SetBranchAddress("CosLep", &a);
  tree->SetBranchAddress("flagCCQE", &CCQEflag); //Set the branch to "flagCCQE"

  TH1F *lepen = new TH1F("lepen", "Lepton Energy", 100, 0, 5000);

  //Moved this loop to the bottom, so that we can fill both histograms at once
//  for(int i = 0; i < entries; i++)
//  {
//    tree->GetEntry(i);
//
//    lepen->Fill(e);
//  }

  //Plot for bins
  //Adding the bins - inspired by Alexis -> I changed it, if awe have a constant bin-size on the x and y axis we can use TH2D, which is easier to deal with.
  //
  //TH2Poly *en_an2 = new TH2Poly();
  //en_an2->SetName("en_an2");
  //en_an2->SetTitle("Partition");

  //In this case the cos\theta bins are all the same ? So why the need of a 2-d histogram ? Can simply do 1-d but only add the \cos\theta < -0.5
  //Maybe a first try to later expand to more bins ?
//  cout << en_an2->AddBin(0, -1, 250, -0.5);
//  cout << en_an2->AddBin(250, -1, 500, -0.5);
//  cout << en_an2->AddBin(500, -1, 750, -0.5);
//  cout << en_an2->AddBin(750, -1, 1000, -0.5);
//  cout << en_an2->AddBin(1000, -1, 1250, -0.5);
//  cout << en_an2->AddBin(1250, -1, 1500, -0.5);
//  cout << en_an2->AddBin(1500, -1, 1750, -0.5);
//  cout << en_an2->AddBin(1750, -1, 2000, -0.5);
//  cout << en_an2->AddBin(1000, -1, 1250, -0.5);
//  cout << en_an2->AddBin(1250, -1, 1500, -0.5);
//  cout << en_an2->AddBin(1500, -1, 1750, -0.5);
//  cout << en_an2->AddBin(1750, -1, 2000, -0.5);
//  cout << en_an2->AddBin(2000, -1, 2250, -0.5);
//  cout << en_an2->AddBin(2250, -1, 2500, -0.5);
//  cout << en_an2->AddBin(2500, -1, 2750, -0.5);
//  cout << en_an2->AddBin(2750, -1, 3000, -0.5);
//  cout << en_an2->AddBin(3000, -1, 3250, -0.5);
//  cout << en_an2->AddBin(3250, -1, 3500, -0.5);
//  cout << en_an2->AddBin(3500, -1, 3750, -0.5);
//  cout << en_an2->AddBin(3750, -1, 4000, -0.5);
//  cout << en_an2->AddBin(4000, -1, 4250, -0.5);
//  cout << en_an2->AddBin(4250, -1, 4500, -0.5);
//  cout << en_an2->AddBin(4500, -1, 4750, -0.5);
//  cout << en_an2->AddBin(4750, -1, 5000, -0.5);

  int n_bins_cos = 20; //With this the bin_width will be 0.10
  double cos_low = -1;
  double cos_up = 1;

  int n_bins_E = 50.;
  double E_low = 105.;
  double E_high = 3000.; //Going up to 3000 instead of larger

  TH2D *en_an2 = new TH2D("2d_E_cos", "Microboone_2d_inclusive", n_bins_cos, cos_low,cos_up, n_bins_E, E_low, E_high);

  //Variables used to calculate the differential cross section
  double n_total = entries; //Total number is just number of entries

  //Unit of totcross ? put a note to remember here : 1e-43 cm^2
  double totcross = events->Integral(); //Integral of histogram



   //FIlling both histograms
   for(int i = 0; i < entries; i++)
   {
     tree->GetEntry(i);

     //Include only the subset of events labeled as CCQE!
     if(CCQEflag == true)
     {
	     lepen->Fill(e);

	     en_an2->Fill(a,e);
     }

   }

  //This is ok, because all your bins have the same width (250 MeV), but would break down if the widths are not the same! I put it in the loop to be more general
  //double bin_width = lepen->GetXaxis()->GetBinWidth(1);


  //The problem: you do a loop over events previously to make a histogram of lepton energy, this is fine, but then why do a new loop over the tree ? this will fill the histogram too many times

  //Loop to pull out values at each bin
  //No need for this: we already filled the histograms above
//  for(int k = 0; k < entries; k++)
 // {
  //  tree->GetEntry(k);

//Name the output file
//
    std::string fileName = "El_bins_flux_0_CCQE.txt";
    std::fstream myFile;
    myFile.open(fileName.c_str(), std::ios::out | std::ios::trunc);

    int n_bins_lepen = lepen->GetXaxis()->GetNbins();
    //Watchout! the bins in the histogram work like this:
    //ibin == 0 -> underflow bin (the number of events that are below the lower bound)
    //ibin == n_bins_lepen + 1 -> overflow bin (the number of events above the upper bound)
    //We will ignore the under and overflow, hence go from ibin = 1 to ibin = n_bins_lepep
    for(int ibin = 1; ibin < n_bins_lepen+1; ibin++)
    {
      double N_bin = lepen->GetBinContent(ibin);
      double bin_width = lepen->GetXaxis()->GetBinWidth(ibin);
      double diff_cross = N_bin/n_total*totcross/bin_width; //Differential cross-section

      //Print it together with the bin edges
      double E_up = lepen->GetXaxis()->GetBinUpEdge(ibin);
      double E_low  = E_up - bin_width;
      myFile << E_low << "  " << E_up << "  " << diff_cross << endl;

      //cout << lepen->GetBinContent(l) << endl;
      //en_an2->Fill(diff_cross, a);
    }

    myFile.close();

 // Now the same for the 2-d histogram
    fileName = "2d_bins_flux_0_CCQE.txt";
    myFile.open(fileName.c_str(), std::ios::out | std::ios::trunc);

    n_bins_cos = en_an2->GetXaxis()->GetNbins();
    n_bins_E = en_an2->GetYaxis()->GetNbins();
   //Same thing with the bins here, both in cosine and lepton energy
   for( int icos = 1; icos < n_bins_cos+1; icos++)
   {
	//Get cosine binwidth and edges:
	 double cos_width = en_an2->GetXaxis()->GetBinWidth(icos);
   double cos_up = en_an2->GetXaxis()->GetBinUpEdge(icos);
	 double cos_low = cos_up - cos_width;

	   //Now for this cosine bin we output the energy distribution:
	 for(int iE = 1; iE < n_bins_E+1; iE++)
	  {
		double E_width = en_an2->GetYaxis()->GetBinWidth(iE);
		double E_up = en_an2->GetYaxis()->GetBinUpEdge(iE);
		double E_low = E_up - E_width;
	  double N_bin = en_an2->GetBinContent(icos,iE);
		double bin_width = cos_width*E_width;
	  double diff_cross = N_bin/n_total*totcross/bin_width; //Differential cross-section

      		myFile << cos_low << "  " << cos_up << "  " << E_low << "  " << E_up << "  " << diff_cross << endl;
	}
	//I add two blank lines between every cosine block because its easier for making plots for me
	myFile << endl;
	myFile << endl;


   }

   myFile.close();

  //Delete the pointers lepen and en_an2
  delete lepen;
  delete en_an2;

  //Delete the TFile pointer
  delete f;

  //fileName = "flux_0_part.txt"; //For later

}
