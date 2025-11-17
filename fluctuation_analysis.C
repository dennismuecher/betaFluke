#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TBox.h>
#include <TLine.h>
#include <iostream>
#include <vector>

TH1D* GaussianSmooth(TH1D* hInput, double sigma, const char* name) {
    int nbins = hInput->GetNbinsX();
    double xmin = hInput->GetXaxis()->GetXmin();
    double xmax = hInput->GetXaxis()->GetXmax();
    double binwidth = hInput->GetBinWidth(1);
    TH1D* hSmooth = new TH1D(name, name, nbins, xmin, xmax);
    int kernelRange = (int)(5.0 * sigma / binwidth);
    for(int i = 1; i <= nbins; i++) {
        double sum = 0.0, normalization = 0.0;
        for(int j = -kernelRange; j <= kernelRange; j++) {
            int binIndex = i + j;
            if(binIndex >= 1 && binIndex <= nbins) {
                double x = j * binwidth;
                double weight = TMath::Exp(-0.5 * x * x / (sigma * sigma));
                sum += hInput->GetBinContent(binIndex) * weight;
                normalization += weight;
            }
        }
        if(normalization > 0) hSmooth->SetBinContent(i, sum / normalization);
    }
    return hSmooth;
}

TGraph* CalculateAutocorrelation(TH1D* hStationary, double Emin, double Emax, int maxShift) {
    int bin_min = hStationary->GetXaxis()->FindBin(Emin);
    int bin_max = hStationary->GetXaxis()->FindBin(Emax);
    double binwidth = hStationary->GetBinWidth(1);
    std::vector<double> epsilon_vals, C_vals;
    for(int shift = 0; shift <= maxShift; shift++) {
        double numerator = 0.0, denom1 = 0.0, denom2 = 0.0;
        int count = 0;
        for(int i = bin_min; i <= bin_max - shift; i++) {
            double d1 = hStationary->GetBinContent(i);
            double d2 = hStationary->GetBinContent(i + shift);
            numerator += d1 * d2;
            denom1 += d1;
            denom2 += d2;
            count++;
        }
        if(count > 0) {
            numerator /= count;
            denom1 /= count;
            denom2 /= count;
            double C = numerator / (denom1 * denom2);
            epsilon_vals.push_back(shift * binwidth);
            C_vals.push_back(C);
        }
    }
    TGraph* grAuto = new TGraph(epsilon_vals.size(), &epsilon_vals[0], &C_vals[0]);
    return grAuto;
}

double GetTrueLevelDensity(TH1D* hInput, double Emin, double Emax) {
    int bin_min = hInput->GetXaxis()->FindBin(Emin);
    int bin_max = hInput->GetXaxis()->FindBin(Emax);
    int nLevels = 0;
    for(int i = bin_min; i <= bin_max; i++) {
        nLevels += (int)hInput->GetBinContent(i);
    }
    return (double)nLevels / (Emax - Emin);
}

void fluctuation_analysis() {
    const char* inputFile = "output1.root";
    const char* histName = "hLevelSpectrum_J5_real0";
    double dE = 0.010;
    double sigma = 0.5 * dE;
    double sigma_ratio = 3.0;
    double sigma_large = sigma_ratio * sigma;
    double alpha = 0.7;
    double f = 1.0 - (sigma/sigma_large)*(sigma/sigma_large);
    
    const int nIntervals = 10;
    double E_start = 4.5;
    double interval_width = 0.5;
    
    TFile* file = TFile::Open(inputFile);
    if(!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open " << inputFile << std::endl;
        return;
    }
    TH1D* hInput = (TH1D*)file->Get(histName);
    if(!hInput) {
        std::cerr << "Error: Cannot find " << histName << std::endl;
        file->Close();
        return;
    }
    hInput->SetDirectory(0);
    file->Close();
    
    std::cout << "Loaded: " << histName << std::endl;
    
    TH1D* hSmooth = GaussianSmooth(hInput, sigma, "g_Ex");
    TH1D* hSmoothLarge = GaussianSmooth(hInput, sigma_large, "g_large_Ex");
    TH1D* hStationary = (TH1D*)hSmoothLarge->Clone("d_Ex");
    hStationary->Divide(hSmooth);
    
    std::vector<double> E_central, rho_fluctuation, rho_true;
    
    for(int i = 0; i < nIntervals; i++) {
        double E_min = E_start + i * interval_width;
        double E_max = E_min + interval_width;
        double E_center = (E_min + E_max) / 2.0;
        
        TGraph* grAuto = CalculateAutocorrelation(hStationary, E_min, E_max, 50);
        double C_0_minus_1 = grAuto->GetPointY(0) - 1.0;
        double mean_spacing = (C_0_minus_1 * 2.0 * sigma * TMath::Sqrt(TMath::Pi())) / (alpha * f);
        double level_density = 1.0 / mean_spacing;
        double true_density = GetTrueLevelDensity(hInput, E_min, E_max);
        
        E_central.push_back(E_center);
        rho_fluctuation.push_back(level_density);
        rho_true.push_back(true_density);
        delete grAuto;
    }
    
    TCanvas* c1 = new TCanvas("c1", "Level Density Comparison", 900, 700);
    TGraph* grFluctuation = new TGraph(nIntervals, &E_central[0], &rho_fluctuation[0]);
    grFluctuation->SetTitle("Level Density;E_{x} (MeV);#rho (states/MeV)");
    grFluctuation->SetMarkerStyle(20);
    grFluctuation->SetMarkerColor(kBlue);
    grFluctuation->SetMarkerSize(1.5);
    grFluctuation->SetLineColor(kBlue);
    grFluctuation->SetLineWidth(2);
    
    TGraph* grTrue = new TGraph(nIntervals, &E_central[0], &rho_true[0]);
    grTrue->SetMarkerStyle(21);
    grTrue->SetMarkerColor(kRed);
    grTrue->SetMarkerSize(1.5);
    grTrue->SetLineColor(kRed);
    grTrue->SetLineWidth(2);
    
    grFluctuation->Draw("APL");
    grTrue->Draw("PL SAME");
    
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(grFluctuation, "Fluctuation analysis", "lp");
    leg->AddEntry(grTrue, "True (from histogram)", "lp");
    leg->Draw();
    
    c1->SaveAs("level_density_comparison.pdf");
}