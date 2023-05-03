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
void xAna_ztoee_2017(string inputtxtFilename = "/afs/cern.ch/user/k/kuanyu/ZtoLL/conder_script/2016BkgMc/filelist_submit/inclusive/2017/mydata.txtah", string outputtxtFilename = "/eos/user/k/kuanyu/2017test/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_9.root")
{
    // cout << "inputtxtFilename = " << inputtxtFilename << endl;
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

    TH1D *h_totevent = new TH1D("h_totevent", "total events", 5, 0, 5);
    h_totevent->Sumw2();

    TH1F *h_total_mcweight = new TH1F("h_total_mcweight", "total MC events", 5, 0, 5);
    h_total_mcweight->Sumw2();

    TH1F *h_total_mcweight_new = new TH1F("h_total_mcweight_new", "total MC events of all files of certain process", 5, 0, 5);
    h_total_mcweight_new->Sumw2();

    TH1D *h_HT_eventCount = new TH1D("h_HT_eventCount", "", 10, 0, 10);
    h_HT_eventCount->SetYTitle("N event");
    h_HT_eventCount->Sumw2();

    TH1F *h_ee_npass = new TH1F("h_ee_npass", "", 10, 0, 10);
    h_ee_npass->SetXTitle("npass");
    h_ee_npass->Sumw2();

    TH1D *h_ee_npass_noweight = new TH1D("h_ee_npass_noweight", "", 10, 0, 10);
    h_ee_npass_noweight->SetXTitle("npass");
    h_ee_npass_noweight->Sumw2();

    /*Void Tree variable*/
    Int_t I_event;
    Int_t I_weight;
    ULong64_t I_eventID;
    Float_t f_Met;
    Float_t f_HT;

    double_t d_dileptonPt;
    double_t d_dileptonMass;
    double_t d_dileptonEta;
    double_t d_dileptonPhi;

    Int_t I_passnAK4Jets;

    vector<float> v_passJetEta;
    vector<float> v_passJetPt;
    vector<float> v_passJetCSV;
    vector<int> v_passjethadronflavor;
    vector<int> v_passJetIndex;

    TString outputfile(outputtxtFilename);

    TFile *outFile = new TFile(outputfile, "RECREATE");

    TTree *tree = new TTree("tree", "Tree");
    tree->Branch("I_event", &I_event);
    tree->Branch("I_weight", &I_weight);
    tree->Branch("I_eventID", &I_eventID);
    tree->Branch("f_Met", &f_Met);
    tree->Branch("f_HT", &f_HT);
    tree->Branch("d_dileptonPt", &d_dileptonPt);
    tree->Branch("d_dileptonMass", &d_dileptonMass);
    tree->Branch("d_dileptonEta", &d_dileptonEta);
    tree->Branch("d_dileptonPhi", &d_dileptonPhi);
    tree->Branch("I_passnAK4Jets", &I_passnAK4Jets);
    tree->Branch("v_passJetEta", &v_passJetEta);
    tree->Branch("v_passJetPt", &v_passJetPt);
    tree->Branch("v_passJetCSV", &v_passJetCSV);
    tree->Branch("v_passjethadronflavor", &v_passjethadronflavor);
    tree->Branch("v_passJetIndex", &v_passJetIndex);

    h_total_mcweight_new->Reset();

    int filenumber = 1;
    int max_filenumber = line;
    while (getline(flist, inputFile))
    {
        cout << "inputFile: " << inputFile << endl;

        if (filenumber < max_filenumber)
        {
            cout << "file number = " << filenumber << endl;
            cout << "process = "
                 << "[" << filenumber * 100.0 / (max_filenumber) << "%"
                 << "]"
                 << "\n"
                 << endl;
        }
        else
        {
            cout << "file number = " << filenumber << endl;
            cout << "finish = "
                 << "[" << filenumber * 100.0 / (max_filenumber) << "%"
                 << "]" << endl;
        }
        filenumber++;

        h_total_mcweight->Reset();

        TString file(inputFile);
        TFile *openfile = TFile::Open(file);

        h_total_mcweight = static_cast<TH1F *>(openfile->Get("h_total_mcweight"));
        h_total_mcweight_new->Add(h_total_mcweight);

        TreeReader data(inputFile.data(), "outTree");
        // Event Loop
        cout << "Total Event =" << data.GetEntriesFast() << endl;

        for (Long64_t jEntry = 0; jEntry < data.GetEntriesFast(); jEntry++)
        {
            if (jEntry % 1000 == 0)
            {
                fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
            }
            v_passJetEta.clear();
            v_passJetPt.clear();
            v_passJetCSV.clear();
            v_passjethadronflavor.clear();
            v_passJetIndex.clear();

            data.GetEntry(jEntry);

            Float_t mcWeight = data.GetFloat("mcweight");
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
            /* Get Total event number*/
            h_totevent->Fill(1, eventWeight);
            /*For inclusive sample event counter*/

            Float_t HT = data.GetFloat("st_HT");
            f_HT = HT;
            if (HT < 70)
            {
                h_HT_eventCount->Fill(1, eventWeight);
            }
            else if (HT >= 70 && HT < 100)
            {
                h_HT_eventCount->Fill(2, eventWeight);
            }
            else if (HT >= 100 && HT < 200)
            {
                h_HT_eventCount->Fill(3, eventWeight);
            }
            else if (HT >= 200 && HT < 400)
            {
                h_HT_eventCount->Fill(4, eventWeight);
            }
            else if (HT >= 400 && HT < 600)
            {
                h_HT_eventCount->Fill(5, eventWeight);
            }
            else if (HT >= 600 && HT < 800)
            {
                h_HT_eventCount->Fill(6, eventWeight);
            }
            else if (HT >= 800 && HT < 1200)
            {
                h_HT_eventCount->Fill(7, eventWeight);
            }
            else if (HT >= 1200 && HT < 2500)
            {
                h_HT_eventCount->Fill(8, eventWeight);
            }
            else if (HT >= 2500)
            {
                h_HT_eventCount->Fill(9, eventWeight);
            }

            /*1. Select electron*/
            Long64_t nEle = data.GetLong64("st_nEle");
            Float_t *elePx = data.GetPtrFloat("st_elePx");
            Float_t *elePy = data.GetPtrFloat("st_elePy");
            Float_t *elePz = data.GetPtrFloat("st_elePz");
            Float_t *eleEnergy = data.GetPtrFloat("st_eleEnergy");
            vector<bool> &eleIsPassMedium = *((vector<bool> *)data.GetPtr("st_eleIsPassMedium"));
            vector<TLorentzVector> goodElectrons;
            goodElectrons.clear();

            for (int ie = 0; ie < nEle; ie++)
            {
                // nEleBefore ++;
                TLorentzVector *myEle = new TLorentzVector(elePx[ie], elePy[ie], elePz[ie], eleEnergy[ie]);
                if (!eleIsPassMedium[ie])
                {
                    continue;
                }
                goodElectrons.push_back(*myEle);

            } // End of Ele Loop
            sort(goodElectrons.begin(), goodElectrons.end(), pt_greater);

            /*2. Select muon*/
            Long64_t nMu = data.GetLong64("st_nMu");
            Float_t *muPx = data.GetPtrFloat("st_muPx");
            Float_t *muPy = data.GetPtrFloat("st_muPy");
            Float_t *muPz = data.GetPtrFloat("st_muPz");
            Float_t *muEnergy = data.GetPtrFloat("st_muEnergy");
            vector<bool> &isSoftLooseIsoMuon = *((vector<bool> *)data.GetPtr("st_isSoftLooseIsoMuon"));
            vector<TLorentzVector> goodmuons;
            goodmuons.clear();

            for (int imu = 0; imu < nMu; imu++)
            {
                TLorentzVector *myMu = new TLorentzVector(muPx[imu], muPy[imu], muPz[imu], muEnergy[imu]);
                // if (muTrkLayers[imu] < 5) continue;
                if (!isSoftLooseIsoMuon[imu])
                {
                    continue;
                }

                goodmuons.push_back(*myMu);

            } // End of Muon Loop
            bool recoeeEvent = false;
            bool recouuEvent = false;
            if (goodmuons.size() == goodElectrons.size())
            {
                continue;
            }
            // if (goodElectrons.size() >= 2 && goodmuons.size() < 2)
            if (goodElectrons.size() == 2 && goodmuons.size() == 0)
            {
                recoeeEvent = true;
            }
            // if (goodmuons.size() >= 2 && goodElectrons.size() < 2)
            if (goodmuons.size() == 2 && goodElectrons.size() == 0)
            {
                recouuEvent = true;
            }
            if (recoeeEvent)
            {
                h_ee_npass->Fill(1, eventWeight);

                Float_t met = data.GetFloat("st_pfMetCorrPt");

                TLorentzVector DilepAfterRecoee = goodElectrons[0] + goodElectrons[1];

                // 3. Good Vertex
                Long64_t nVtx = data.GetLong64("st_nVtx");
                if (nVtx < 1)
                    continue;
                h_ee_npass->Fill(2, eventWeight);

                // 4. Tau Veto
                Long64_t nTau_DRBased_EleVeto = data.GetLong64("st_nTau_DRBased_EleVeto");
                if (nTau_DRBased_EleVeto > 0)
                    continue;
                h_ee_npass->Fill(3, eventWeight);

                // 5. Z boson
                if (goodElectrons[0].Pt() <= 25 && goodElectrons[1].Pt() <= 20)
                {
                    continue;
                }
                h_ee_npass->Fill(4, eventWeight);

                float PDGZmass = 91.1876;
                float dilepMass = DilepAfterRecoee.M();
                float deltaMass = abs(PDGZmass - dilepMass);
                if (deltaMass > 15)
                {
                    continue;
                }
                h_ee_npass->Fill(5, eventWeight);

                // 6. Thin Jet
                Long64_t nTHINjets = data.GetLong64("st_THINnJet");
                Float_t *THINjetPx = data.GetPtrFloat("st_THINjetPx");
                Float_t *THINjetPy = data.GetPtrFloat("st_THINjetPy");
                Float_t *THINjetPz = data.GetPtrFloat("st_THINjetPz");
                Float_t *THINjetEnergy = data.GetPtrFloat("st_THINjetEnergy");
                Int_t *THINjetHadronFlavor = data.GetPtrInt("st_THINjetHadronFlavor");
                vector<int> indexForPassAK4;
                indexForPassAK4.clear();

                if (nTHINjets < 1)
                    continue;
                /*
                for (int itj = 0; itj < nTHINjets; itj++)
                {
                    TLorentzVector *thisTHINjet = new TLorentzVector(THINjetPx[itj], THINjetPy[itj], THINjetPz[itj], THINjetEnergy[itj]);
                    goodTHINjets.push_back(*thisTHINjet);
                }*/
                h_ee_npass->Fill(7, eventWeight);
                
                I_event = jEntry;
                I_weight = eventWeight;

                f_Met = met;
                d_dileptonPt = DilepAfterRecoee.Pt();
                d_dileptonMass = dilepMass;
                d_dileptonEta = DilepAfterRecoee.Eta();
                d_dileptonPhi = DilepAfterRecoee.Phi();

                tree->Fill();
            } // End of rece EE
        }     // End of loop all files
    }         // End of all files loop

    outFile->cd();
    tree->Write();
    outFile->Close();
    cout << "output written to " << outputfile << endl;
} // Big Scope End
