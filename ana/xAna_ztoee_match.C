#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TAxis.h>
#include <math.h>
#include <string>
using namespace std;
void efferr(float nsig, float ntotal, float factor = 1)
{
    float eff = nsig / ntotal;
    float err = sqrt((1 - eff) * eff / ntotal);
    cout << "efficiency = " << eff * factor << " +- " << err * factor << endl;
}
bool pt_greater(const TLorentzVector a, const TLorentzVector b)
{
    double A = a.Pt();
    double B = b.Pt();
    return (A > B);
}
double median_value(vector<double> tmpvector)
{
    float med_value = 0.0;
    sort(tmpvector.begin(), tmpvector.end());

    if (tmpvector.size() % 2 == 0) // even
    {
        med_value = (tmpvector[tmpvector.size() / 2 - 1] + tmpvector[tmpvector.size() / 2]) / 2;
    }
    else // odd
    {
        med_value = tmpvector[tmpvector.size() / 2];
    }
    return (med_value);
}
double mean_value(const vector<double> &ttmpvector)
{
    float sum_value = 0.0;
    for (float x : ttmpvector)
    {
        sum_value += x;
    }
    return (sum_value / ttmpvector.size());
}
float cal_dphi(float phi1, float phi2)
{
    float dphi = phi1 - phi2;
    while (dphi >= TMath::Pi())
        dphi -= 2 * TMath::Pi();
    while (dphi < -TMath::Pi())
        dphi += 2 * TMath::Pi();
    return TMath::Abs(dphi);
}
// void xAna_ztoee(string inputtxtFilename = "root://se01.grid.nchc.org.tw//dpm/grid.nchc.org.tw/home/cms/store/user/syu/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Imp_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-ext2/210929_173914/0004/NCUGlobalTuples_4485.root", string outputfile = "/eos/user/k/kuanyu/data/DY_inclusive/DY_inclusive_9423.root")
// void xAna_ztoee()
/// afs/cern.ch/user/k/kuanyu/test/tmp150.txt
void xAna_ztoee_match(string inputtxtFilename = "/afs/cern.ch/user/k/kuanyu/test/tmp150.txt", string outputfile = "./test_ee_Mx150_match.root")
{
    // string inputFile(inputtxtFilename.data());
    cout << "inputtxtFilename = " << inputtxtFilename << endl;
    ifstream flist(inputtxtFilename.data());
    string inputFile;

    fstream fin(inputtxtFilename, ios::in);
    if (!fin)
    {
        // cout << "bug" << endl;
    }
    string tmps;
    int line = 0;
    while (getline(fin, tmps))
    {
        line++;
    }
    fin.close();
    cout << "line" << line << endl;

    //------------------
    // Create histrogram
    //------------------
    TH1D *h_genee_event = new TH1D("h_genee_event", "gen events", 5, 0, 5);
    h_genee_event->Sumw2();

    TH1D *h_recoee_event = new TH1D("h_recoee_event", "reco events", 5, 0, 5);
    h_recoee_event->Sumw2();

    TH1D *h_totevent = new TH1D("h_totevent", "total events", 5, 0, 5);
    h_totevent->Sumw2();

    TH1F *h_HT_eventCout = new TH1F("h_HT_eventCout", "", 10, 0, 10);
    h_HT_eventCout->SetYTitle("N event");
    h_HT_eventCout->Sumw2();

    TH1F *h_ee_npass = new TH1F("h_ee_npass", "", 10, 0, 10);
    h_ee_npass->SetXTitle("npass");
    h_ee_npass->Sumw2();

    TH1D *gen_chi2numb = new TH1D("gen_chi2numb", "Chi2", 5, 0, 5);
    TH1D *gen_dquarknumb = new TH1D("gen_dquarknumb", "d quark", 10, 0, 10);
    TH1D *gen_eenumber = new TH1D("gen_eenumber", "e", 10, 0, 10);
    TH1D *match_dquarknumb = new TH1D("match_dquarknumb", "d quark", 10, 0, 10);

    TH1D *h_ele_n = new TH1D("h_ele_n", "Nr of electron", 10, 0, 10);
    h_ele_n->GetXaxis()->SetTitle("Number of electron");
    h_ele_n->GetYaxis()->SetTitle("Number of Events");
    h_ele_n->Sumw2();

    TH1D *h_tau_n = new TH1D("h_tau_n", "Nr of tau", 10, 0, 10);
    h_tau_n->GetXaxis()->SetTitle("Number of tau");
    h_tau_n->GetYaxis()->SetTitle("Number of Events");
    h_tau_n->Sumw2();

    TH1D *h_mu_n = new TH1D("h_mu_n", "Nr of Muon", 10, 0, 10);
    h_mu_n->GetXaxis()->SetTitle("Number of muon");
    h_mu_n->GetYaxis()->SetTitle("Number of Events");
    h_mu_n->Sumw2();

    TH1D *Z_eemass = new TH1D("Z_eemass", "Z->ee", 150, 0, 150);
    Z_eemass->GetXaxis()->SetTitle("Z mass");
    Z_eemass->GetYaxis()->SetTitle("Number of Events");
    Z_eemass->Sumw2();

    TH1D *dilepton_pT = new TH1D("dilepton_pT", "dilepton_pT", 100, 0, 1000);
    dilepton_pT->GetXaxis()->SetTitle("dilepton pT");
    dilepton_pT->GetYaxis()->SetTitle("Number of Events");
    dilepton_pT->Sumw2();

    TH1D *dilepton_pT_after_ledptcut = new TH1D("dilepton_pT_after_ledptcut", "dilepton_pT_after_ledptcut", 100, 0, 1000);
    dilepton_pT_after_ledptcut->GetXaxis()->SetTitle("dilepton pT");
    dilepton_pT_after_ledptcut->GetYaxis()->SetTitle("Number of Events");
    dilepton_pT_after_ledptcut->Sumw2();

    TH1D *dilepton_pT_after_dilepmasscut = new TH1D("dilepton_pT_after_dilepmasscut", "dilepton_pT_after_dilepmasscut", 100, 0, 1000);
    dilepton_pT_after_dilepmasscut->GetXaxis()->SetTitle("dilepton pT");
    dilepton_pT_after_dilepmasscut->GetYaxis()->SetTitle("Number of Events");
    dilepton_pT_after_dilepmasscut->Sumw2();

    TH1D *dilepton_pT_after_extralepcut = new TH1D("dilepton_pT_after_extralepcut", "dilepton_pT_after_extralepcut", 100, 0, 1000);
    dilepton_pT_after_extralepcut->GetXaxis()->SetTitle("dilepton pT");
    dilepton_pT_after_extralepcut->GetYaxis()->SetTitle("Number of Events");
    dilepton_pT_after_extralepcut->Sumw2();

    TH1D *dilepton_pT_after_vetotaucut = new TH1D("dilepton_pT_after_vetotaucut", "dilepton_pT_after_extralepcut", 100, 0, 1000);
    dilepton_pT_after_vetotaucut->GetXaxis()->SetTitle("dilepton pT");
    dilepton_pT_after_vetotaucut->GetYaxis()->SetTitle("Number of Events");
    dilepton_pT_after_vetotaucut->Sumw2();

    TH1D *dilepton_pT_after_nJetcut = new TH1D("dilepton_pT_after_nJetcut", "dilepton_pT_after_nJetcut", 100, 0, 1000);
    dilepton_pT_after_nJetcut->GetXaxis()->SetTitle("dilepton pT");
    dilepton_pT_after_nJetcut->GetYaxis()->SetTitle("Number of Events");
    dilepton_pT_after_nJetcut->Sumw2();

    TH1D *h_njet = new TH1D("h_njet", "Nr of jets pass jet preselection", 15, 0, 15);
    h_njet->GetXaxis()->SetTitle("Number of Jets");
    h_njet->GetYaxis()->SetTitle("Number of Events");
    h_njet->Sumw2();

    TH1D *h_njet_pass_trks = new TH1D("h_njet_pass_trks", "Nr of jets pass jet preselection + track selection", 15, 0, 15);
    h_njet_pass_trks->GetXaxis()->SetTitle("Number of Jets");
    h_njet_pass_trks->GetYaxis()->SetTitle("Number of Events");
    h_njet_pass_trks->Sumw2();

    TH1D *h_jet_csv = new TH1D("h_jet_csv", "", 120, 11, 1.);
    h_jet_csv->GetXaxis()->SetTitle("JetCSV");
    h_jet_csv->GetYaxis()->SetTitle("nJet");
    h_jet_csv->Sumw2();

    TH1D *h_jet_rank = new TH1D("h_jet_rank", "Nr of rank Matched jets", 15, 0, 15);
    h_jet_rank->GetXaxis()->SetTitle("Jet_rank");
    h_jet_rank->GetYaxis()->SetTitle("Number of Events");
    h_jet_rank->Sumw2();

    TH1D *h_trk_npass = new TH1D("h_trk_npass", "", 10, 0, 10);
    h_trk_npass->Sumw2();

    TH1D *h_trkpt = new TH1D("h_trkpt", "", 1000, 0, 1000);
    h_trkpt->Sumw2();

    //----------------------
    // Volid Tree variable
    //----------------------
    Int_t I_event;
    Int_t I_weight;
    ULong64_t I_eventID;
    Int_t I_tot_Recoevt;
    Int_t I_tot_gencoevt;
    Int_t I_Sumeventweight;
    Float_t f_HT;
    Float_t f_Met;
    Double_t d_Met_leptons_ratio;
    Float_t f_dileptonmass;
    Float_t f_dileptonPT;

    vector<float> v_MetdeltaPhi;
    /*electron match*/
    vector<int> v_match_parid;
    vector<int> v_match_momparid;
    //-----------------
    // Jets variable
    //-----------------
    Int_t I_nThinJets;
    vector<float> f_JetEta;
    vector<float> f_JetPt;
    vector<float> f_matchJet_PT;
    vector<float> f_matchJet_Eta;
    vector<float> f_thinjetCSV;
    Float_t f_leadJeteta;
    Float_t f_leadJetpt;
    Int_t I_leadJetflavor;
    Float_t f_minalphaJeteta;
    Float_t f_minalphaJetpt;
    Int_t I_minJethadronflavor;
    Int_t I_minJetpartonflavor;
    Int_t I_nJets;

    //-----------------
    // Tracks variable
    //-----------------

    vector<int> v_N_Tracks;
    vector<float> v_TrackPT;
    vector<double> v_IP2D;
    vector<float> v_IP2DxTrackPT;
    vector<float> v_TrackEta;
    vector<float> v_Trackdr;
    vector<int> v_Trackindex;
    vector<double> v_Chi3DlogPaper;
    vector<double> v_Chi3DPaper;
    vector<double> v_Chi3Dlog;
    vector<double> v_Chi3D;
    vector<double> v_Chi2Dlog;

    vector<double> v_Median_log3DIPsig;
    vector<double> v_Mean_log3DIPsig;
    vector<double> v_Median_log2DIPsig;
    vector<double> v_Median_2DIPsig;
    Double_t f_alphamax;
    Double_t f_alphamin;
    vector<int> v_N_Trk_cut3Dsig;
    vector<float> v_fakeJetEta;
    vector<float> v_fakeJetPt;
    vector<float> v_fakeJetMass;
    vector<float> v_fakeJetCSV;
    vector<float> v_fakealpha;
    vector<float> v_fakealpha2;
    vector<float> v_fakealpha3;
    vector<float> v_fakealpha4;
    vector<float> v_fakeJethadronflavor;
    vector<float> v_fakeJetpartonflavor;
    vector<float> v_fakeJet_ledLepdr;
    vector<float> v_fakeJet_subLepdr;

    vector<float> v_bflavorIP2D;
    vector<float> v_bflavor3Dsig;
    vector<float> v_cflavorIP2D;
    vector<float> v_cflavor3Dsig;
    vector<float> v_lightIP2D;
    vector<float> v_light3Dsig;

    TTree *T_tree = new TTree("T_tree", "Tree");
    T_tree->Branch("I_event", &I_event);
    T_tree->Branch("I_weight", &I_weight);
    T_tree->Branch("I_eventID", &I_eventID);
    T_tree->Branch("I_tot_gencoevt", &I_tot_gencoevt);
    T_tree->Branch("I_tot_Recoevt", &I_tot_Recoevt);
    T_tree->Branch("I_Sumeventweight", &I_Sumeventweight);
    T_tree->Branch("f_Met", &f_Met);
    T_tree->Branch("d_Met_leptons_ratio", &d_Met_leptons_ratio);
    T_tree->Branch("f_HT", &f_HT);
    T_tree->Branch("f_dileptonmass", &f_dileptonmass);
    T_tree->Branch("f_dileptonPT", &f_dileptonPT);
    T_tree->Branch("v_MetdeltaPhi", &v_MetdeltaPhi);
    T_tree->Branch("v_match_parid", &v_match_parid);
    T_tree->Branch("v_match_momparid", &v_match_momparid);
    T_tree->Branch("I_nThinJets", &I_nThinJets);
    T_tree->Branch("f_JetEta", &f_JetEta);
    T_tree->Branch("f_JetPt", &f_JetPt);
    T_tree->Branch("f_matchJet_PT", &f_matchJet_PT);
    T_tree->Branch("f_matchJet_Eta", &f_matchJet_Eta);
    T_tree->Branch("f_thinjetCSV", &f_thinjetCSV);
    T_tree->Branch("f_leadJetpt", &f_leadJetpt);
    T_tree->Branch("f_leadJeteta", &f_leadJeteta);
    T_tree->Branch("I_leadJetflavor", &I_leadJetflavor);
    T_tree->Branch("f_minalphaJetpt", &f_minalphaJetpt);
    T_tree->Branch("f_minalphaJeteta", &f_minalphaJeteta);
    T_tree->Branch("I_minJethadronflavor", &I_minJethadronflavor);
    T_tree->Branch("I_minJetpartonflavor", &I_minJetpartonflavor);
    T_tree->Branch("I_nJets", &I_nJets);
    T_tree->Branch("v_N_Tracks", &v_N_Tracks);
    T_tree->Branch("v_TrackPT", &v_TrackPT);
    T_tree->Branch("v_IP2DxTrackPT", &v_IP2DxTrackPT);
    T_tree->Branch("v_TrackEta", &v_TrackEta);
    T_tree->Branch("v_Trackdr", &v_Trackdr);
    T_tree->Branch("v_Trackindex", &v_Trackindex);
    T_tree->Branch("v_IP2D", &v_IP2D);
    T_tree->Branch("v_Chi2Dlog", &v_Chi2Dlog);
    T_tree->Branch("v_Chi3Dlog", &v_Chi3Dlog);
    T_tree->Branch("v_Chi3D", &v_Chi3D);
    T_tree->Branch("v_Chi3DlogPaper", &v_Chi3DlogPaper);
    T_tree->Branch("v_Chi3DPaper", &v_Chi3DPaper);
    T_tree->Branch("v_Median_log3DIPsig", &v_Median_log3DIPsig);
    T_tree->Branch("v_Mean_log3DIPsig", &v_Mean_log3DIPsig);
    T_tree->Branch("v_Median_log2DIPsig", &v_Median_log2DIPsig);
    T_tree->Branch("v_Median_2DIPsig", &v_Median_2DIPsig);
    T_tree->Branch("f_alphamax", &f_alphamax);
    T_tree->Branch("f_alphamin", &f_alphamin);
    T_tree->Branch("v_N_Trk_cut3Dsig", &v_N_Trk_cut3Dsig);
    T_tree->Branch("v_fakeJetPt", &v_fakeJetPt);
    T_tree->Branch("v_fakeJetEta", &v_fakeJetEta);
    T_tree->Branch("v_fakeJetMass", &v_fakeJetMass);
    T_tree->Branch("v_fakeJetCSV", &v_fakeJetCSV);
    T_tree->Branch("v_fakealpha", &v_fakealpha);
    T_tree->Branch("v_fakealpha2", &v_fakealpha2);
    T_tree->Branch("v_fakealpha3", &v_fakealpha3);
    T_tree->Branch("v_fakealpha4", &v_fakealpha4);
    T_tree->Branch("v_fakeJethadronflavor", &v_fakeJethadronflavor);
    T_tree->Branch("v_fakeJetpartonflavor", &v_fakeJetpartonflavor);
    T_tree->Branch("v_fakeJet_ledLepdr", &v_fakeJet_ledLepdr);
    T_tree->Branch("v_fakeJet_subLepdr", &v_fakeJet_subLepdr);
    T_tree->Branch("v_bflavorIP2D", &v_bflavorIP2D);
    T_tree->Branch("v_bflavor3Dsig", &v_bflavor3Dsig);
    T_tree->Branch("v_cflavorIP2D", &v_cflavorIP2D);
    T_tree->Branch("v_cflavor3Dsig", &v_cflavor3Dsig);
    T_tree->Branch("v_lightIP2D", &v_lightIP2D);
    T_tree->Branch("v_light3Dsig", &v_light3Dsig);

    int filenumber = 0;
    int max_filenumber = line;
    while (getline(flist, inputFile))
    {
        cout << "inputFile = " << inputFile << endl;
        // if (filenumber < max_filenumber - 1)
        //{
        //     cout << "file number = " << filenumber << endl;
        //     cout << "process = "
        //          << "[" << filenumber * 100.0 / (max_filenumber - 1) << "%"
        //          << "]"
        //          << "\n"
        //          << endl;
        // }
        // else
        //{
        //     cout << "file number = " << filenumber << endl;
        //     cout << "finish = "
        //          << "[" << filenumber * 100.0 / (max_filenumber - 1) << "%"
        //          << "]" << endl;
        // }
        filenumber++;
        bool BKG2016MC = true;
        //  Identify input file type
        if (((inputFile).find("DYJetsToLL") != std::string::npos) ||
            ((inputFile).find("ZZTo") != std::string::npos) ||
            ((inputFile).find("WWTo") != std::string::npos) ||
            ((inputFile).find("WZTo") != std::string::npos) ||
            ((inputFile).find("top_5f") != std::string::npos) ||
            ((inputFile).find("TTTo") != std::string::npos) ||
            ((inputFile).find("TTWJetsTo") != std::string::npos) ||
            ((inputFile).find("TTZTo") != std::string::npos) ||
            ((inputFile).find("WWZ_") != std::string::npos) ||
            ((inputFile).find("WZZ_") != std::string::npos) ||
            ((inputFile).find("ZZZ_") != std::string::npos))
        {
            BKG2016MC = true;
        }

        // cout << "BKG2016MC before = " << BKG2016MC << endl;
        TreeReader data(inputFile.data());
        for (Long64_t jEntry = 0; jEntry < data.GetEntriesFast(); jEntry++)
        // for (Long64_t jEntry = 0; jEntry < 5; jEntry++)
        {
            // if (jEntry % 2000 == 0)
            //{
            //     fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
            // }
            //-------------------
            //  void some variable
            //-------------------
            Int_t neeRecoTotal = 0;
            Int_t neeGenevt = 0;
            Int_t SumeventWeight = 0;
            //  clear Tree vector for each event
            f_JetEta.clear();
            f_JetPt.clear();
            f_matchJet_PT.clear();
            f_matchJet_Eta.clear();
            f_thinjetCSV.clear();
            v_Chi3Dlog.clear();
            v_Chi3D.clear();
            v_Chi3DlogPaper.clear();
            v_Chi3DPaper.clear();
            v_Chi2Dlog.clear();
            v_IP2D.clear();
            v_Trackindex.clear();
            v_N_Tracks.clear();
            v_N_Trk_cut3Dsig.clear();
            v_TrackPT.clear();
            v_IP2DxTrackPT.clear();
            v_TrackEta.clear();
            v_fakeJetPt.clear();
            v_fakeJetEta.clear();
            v_fakeJetMass.clear();
            v_fakeJetCSV.clear();
            v_fakealpha.clear();
            v_fakealpha2.clear();
            v_fakealpha3.clear();
            v_fakealpha4.clear();
            v_fakeJethadronflavor.clear();
            v_fakeJetpartonflavor.clear();
            v_fakeJet_ledLepdr.clear();
            v_fakeJet_subLepdr.clear();
            v_Trackdr.clear();
            v_MetdeltaPhi.clear();
            v_match_parid.clear();
            v_match_momparid.clear();
            v_Median_log3DIPsig.clear();
            v_Mean_log3DIPsig.clear();
            v_Median_log2DIPsig.clear();
            v_Median_2DIPsig.clear();
            v_bflavorIP2D.clear();
            v_bflavor3Dsig.clear();
            v_cflavorIP2D.clear();
            v_cflavor3Dsig.clear();
            v_lightIP2D.clear();
            v_light3Dsig.clear();

            data.GetEntry(jEntry);
            Float_t mcWeight = data.GetFloat("mcWeight");
            Double_t eventWeight = mcWeight;
            if (eventWeight > 0)
            {
                eventWeight = 1;
            }
            else if (eventWeight < 0)
            {
                eventWeight = -1;
            }
            else
            {
                eventWeight = 1;
            }
            SumeventWeight += eventWeight;
            //-----------------------------------
            // For inclusive sample event counter
            //-----------------------------------
            Float_t HT = data.GetFloat("HT");
            if (HT < 70)
            {
                h_HT_eventCout->Fill(1, eventWeight);
            }
            else if (HT >= 70 && HT < 100)
            {
                h_HT_eventCout->Fill(2, eventWeight);
            }
            else if (HT >= 100 && HT < 200)
            {
                h_HT_eventCout->Fill(3, eventWeight);
            }
            else if (HT >= 200 && HT < 400)
            {
                h_HT_eventCout->Fill(4, eventWeight);
            }
            else if (HT >= 400 && HT < 600)
            {
                h_HT_eventCout->Fill(5, eventWeight);
            }
            else if (HT >= 600 && HT < 800)
            {
                h_HT_eventCout->Fill(6, eventWeight);
            }
            else if (HT >= 800 && HT < 1200)
            {
                h_HT_eventCout->Fill(7, eventWeight);
            }
            else if (HT >= 1200 && HT < 2500)
            {
                h_HT_eventCout->Fill(8, eventWeight);
            }
            else if (HT >= 2500)
            {
                h_HT_eventCout->Fill(9, eventWeight);
            }
            //---------------------------
            // Store Total event number
            //---------------------------
            h_totevent->Fill(1.0, eventWeight);

            // h_nlototevent->Fill(eventWeight, eventWeight);

            bool matchee = false;
            vector<TLorentzVector> myEles;
            vector<TLorentzVector> dquark;
            vector<TLorentzVector> chi2s;
            dquark.clear();
            myEles.clear();
            chi2s.clear();
            // 0. check the generator-level information and make sure there is a Z->e+e-
            Int_t nGenPar = data.GetInt("nGenPar");
            TClonesArray *genParP4 = (TClonesArray *)data.GetPtrTObject("genParP4");
            Int_t *genParId = data.GetPtrInt("genParId");
            Int_t *genParSt = data.GetPtrInt("genParSt");
            Int_t *genMomParId = data.GetPtrInt("genMomParId");
            for (int ig = 0; ig < nGenPar; ig++)
            {
                TLorentzVector *thisGen = (TLorentzVector *)genParP4->At(ig);
                int pid = genParId[ig];
                int mompid = genMomParId[ig];
                int status = genParSt[ig];
                if (abs(pid) == 11 && mompid == 23)
                {
                    matchee = true;
                    myEles.push_back(*thisGen);
                }
                // chi2Id 18
                if (abs(pid) == 1 && abs(mompid) == 18)
                {
                    dquark.push_back(*thisGen);
                }
                if (abs(pid) == 18)
                {
                    chi2s.push_back(*thisGen);
                }
            }
            gen_dquarknumb->Fill(dquark.size(), eventWeight);
            gen_chi2numb->Fill(chi2s.size(), eventWeight);
            gen_eenumber->Fill(myEles.size(), eventWeight);

            if (BKG2016MC)
            {
                matchee = true;
            }
            if (matchee)
            {
                neeGenevt++;

                // 1. electron
                int nEle = data.GetInt("nEle");
                TClonesArray *eleP4 = (TClonesArray *)data.GetPtrTObject("eleP4");
                vector<bool> &eleIsPassLoose = *((vector<bool> *)data.GetPtr("eleIsPassLoose"));
                vector<bool> &eleIsPassMedium = *((vector<bool> *)data.GetPtr("eleIsPassMedium"));
                vector<bool> &eleIsPassVeto = *((vector<bool> *)data.GetPtr("eleIsPassVeto"));
                float *eleChHadIso = data.GetPtrFloat("eleChHadIso");
                float *eleNeHadIso = data.GetPtrFloat("eleNeHadIso");
                float *eleGamIso = data.GetPtrFloat("eleGamIso");
                float *elePUPt = data.GetPtrFloat("elePUPt");
                vector<TLorentzVector> goodElectrons;
                goodElectrons.clear();
                vector<int> vetoee;
                vetoee.clear();
                for (int ie = 0; ie < nEle; ie++)
                {
                    TLorentzVector *myEle = (TLorentzVector *)eleP4->At(ie);
                    double elepT = myEle->Pt();
                    if (myEle->Pt() <= 20)
                        continue;
                    if (fabs(myEle->Eta()) >= 2.4)
                        continue;
                    if (!eleIsPassMedium[ie])
                        continue;

                    goodElectrons.push_back(*myEle);
                } // End of loop nEle event

                h_ele_n->Fill(goodElectrons.size(), eventWeight);
                // Sort electron by PT
                sort(goodElectrons.begin(), goodElectrons.end(), pt_greater);

                // 2. Muon
                int nMu = data.GetInt("nMu");
                TClonesArray *muP4 = (TClonesArray *)data.GetPtrTObject("muP4");
                vector<bool> &isTightMuon = *((vector<bool> *)data.GetPtr("isTightMuon"));
                vector<bool> &isSoftMuon = *((vector<bool> *)data.GetPtr("isSoftMuon"));
                float *muChHadIso = data.GetPtrFloat("muChHadIso");
                float *muNeHadIso = data.GetPtrFloat("muNeHadIso");
                float *muGamIso = data.GetPtrFloat("muGamIso");
                float *muPUPt = data.GetPtrFloat("muPUPt");
                int *muTrkLayers = data.GetPtrInt("muTrkLayers");
                vector<TLorentzVector> goodmuons;
                goodmuons.clear();
                for (int imu = 0; imu < nMu; imu++)
                {
                    TLorentzVector *myMu = (TLorentzVector *)muP4->At(imu);
                    double mupT = myMu->Pt();
                    if (myMu->Pt() <= 20)
                        continue;
                    if (fabs(myMu->Eta()) >= 2.4)
                        continue;
                    double myMuIso = (muChHadIso[imu] + max(0., muNeHadIso[imu] + muGamIso[imu] - 0.5 * muPUPt[imu])) / mupT;
                    if (myMuIso > 0.15) //~~~!!!!
                        continue;
                    if (muTrkLayers[imu] < 5)
                        continue;
                    if (!isSoftMuon[imu])
                        continue;
                    goodmuons.push_back(*myMu);
                } // End of Muon loop
                h_mu_n->Fill(goodmuons.size(), eventWeight);

                bool recoeeEvent = false;
                bool recouuEvent = false;
                if (goodmuons.size() == goodElectrons.size())
                {
                    continue;
                }
                if (goodElectrons.size() >= 2 && goodmuons.size() < 2)
                {
                    recoeeEvent = true;
                }
                if (goodmuons.size() >= 2 && goodElectrons.size() < 2)
                {
                    recouuEvent = true;
                }
                if (recoeeEvent)
                {
                    // 3. has a good vertex
                    int nVtx = data.GetInt("nVtx");
                    // 4. veto tau
                    int nTau = data.GetInt("HPSTau_n");
                    TClonesArray *tauP4 = (TClonesArray *)data.GetPtrTObject("HPSTau_4Momentum");
                    vector<bool> &disc_decayModeFindingNewDMs = *((vector<bool> *)data.GetPtr("disc_decayModeFindingNewDMs"));
                    vector<bool> &disc_byVTightIsolationMVA3newDMwLT = *((vector<bool> *)data.GetPtr("disc_byVTightIsolationMVA3newDMwLT"));
                    vector<TLorentzVector> goodtau;
                    goodtau.clear();
                    for (int it = 0; it < nTau; it++)
                    {
                        TLorentzVector *myTau = (TLorentzVector *)tauP4->At(it);
                        if (myTau->Pt() <= 20)
                            continue;
                        if (fabs(myTau->Eta()) >= 2.5)
                            continue;
                        if (!disc_decayModeFindingNewDMs[it])
                            continue;
                        if (!disc_byVTightIsolationMVA3newDMwLT[it])
                            continue;
                        goodtau.push_back(*myTau);
                    } // End of tau loop
                      // 5. For Met
                    Float_t met = data.GetFloat("pfMetCorrPt");
                    // 6. Z boson
                    TLorentzVector Z_boson_ee;
                    Z_boson_ee = goodElectrons[0] + goodElectrons[1];
                    float PDGZmass = 91.1876;
                    float dilepPt = Z_boson_ee.Pt();
                    float dilepMass = Z_boson_ee.M();
                    /*Study Match ee particle ID*/
                    for (int iee = 0; iee < goodElectrons.size(); iee++)
                    {
                        double eleEta = goodElectrons[iee].Eta();
                        double elePhi = goodElectrons[iee].Phi();
                        for (int ig = 0; ig < nGenPar; ig++)
                        {
                            TLorentzVector *thisGen = (TLorentzVector *)genParP4->At(ig);
                            int pid = genParId[ig];
                            int mompid = genMomParId[ig];
                            double genParEta = thisGen->Eta();
                            double genParPhi = thisGen->Phi();
                            /*Calculate delta R*/
                            double diff_Eta = eleEta - genParEta;
                            double diff_Phi = elePhi - genParPhi;
                            double deltaR = sqrt(pow(diff_Eta, 2) + pow(diff_Phi, 2));
                            if (deltaR < 0.1)
                            {
                                v_match_parid.push_back(pid);
                                v_match_momparid.push_back(mompid);
                            }
                            continue;
                        }
                        continue;
                    } // End of matching ele loop
                    
                    /*Fill Tree event variable*/
                    f_HT = HT;
                    f_Met = met;
                    f_dileptonPT = dilepPt;
                    f_dileptonmass = dilepMass;
                    T_tree->Fill();
                } // recoEE
            }     // GenEE
        }         // End of loop over entries
    }             // End of loop all files
    // out Tree branches
    TFile *outFile = new TFile(outputfile.c_str(), "RECREATE");
    outFile->cd();
    T_tree->Write();
    outFile->mkdir("Event_Variable", "Event_Variable")->cd();
    h_totevent->Write();
    h_genee_event->Write();
    h_recoee_event->Write();
    h_ele_n->Write();
    h_mu_n->Write();
    h_tau_n->Write();
    Z_eemass->Write();
    dilepton_pT->Write();
    dilepton_pT_after_ledptcut->Write();
    dilepton_pT_after_dilepmasscut->Write();
    dilepton_pT_after_extralepcut->Write();
    dilepton_pT_after_vetotaucut->Write();
    dilepton_pT_after_nJetcut->Write();
    h_ee_npass->Write();
    gen_chi2numb->Write();
    gen_dquarknumb->Write();
    gen_eenumber->Write();
    match_dquarknumb->Write();
    h_HT_eventCout->Write();
    outFile->cd("/");
    outFile->mkdir("Jet_Variable", "Jet_Variable")->cd();
    h_jet_rank->Write();
    h_njet->Write();
    h_njet_pass_trks->Write();
    h_jet_csv->Write();
    outFile->cd("/");
    outFile->mkdir("Track_Variable", "Track_Variable")->cd();
    h_trk_npass->Write();
    h_trkpt->Write();
    outFile->cd("/");
    outFile->Close();
}
