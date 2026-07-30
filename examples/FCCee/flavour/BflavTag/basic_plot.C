#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TString.h>
#include <TDirectory.h>
#include <TCut.h>
#include <TPad.h>
#include <TStyle.h>
#include <algorithm>

void ComparisonPlot(TTree* tree, 
                    TString var1_spec, 
                    TString label1, 
                    TString var2_spec, 
                    TString label2, 
                    TString xAxisTitle = "x-axis label",
                    TString savePath = "", 
                    TString cutStr1 = "", 
                    TString cutStr2 = "") 
{
    if (!tree) {
        printf("Error: Input TTree pointer is null!\n");
        return;
    }

    // Hide stats box
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // Extract histogram names dynamically ("Rec_vtx_d2PV_mag>>h1(50,0,30)" -> "h1")
    TString hname1 = "h1";
    TString hname2 = "h2";
    
    if (var1_spec.Contains(">>")) {
        hname1 = var1_spec(var1_spec.Index(">>") + 2, var1_spec.Length());
        if (hname1.Contains("(")) hname1 = hname1(0, hname1.Index("("));
    }
    if (var2_spec.Contains(">>")) {
        hname2 = var2_spec(var2_spec.Index(">>") + 2, var2_spec.Length());
        if (hname2.Contains("(")) hname2 = hname2(0, hname2.Index("("));
    }

    // Convert cuts to explicit TCut objects using .Data() to avoid Cling file-lookup bugs
    TCut cut1 = (cutStr1.IsNull() || cutStr1 == "\"\"") ? TCut("") : TCut(cutStr1.Data());
    TCut cut2 = (cutStr2.IsNull() || cutStr2 == "\"\"") ? TCut("") : TCut(cutStr2.Data());

    // Ensure active canvas exists
    TCanvas *c = (TCanvas*)gPad;
    if (!c) {
        c = new TCanvas("c_comp", "Comparison Plot", 800, 600);
    }
    c->cd();

    //Draw first variable
    tree->Draw(var1_spec.Data(), cut1);
    TH1F *h1 = (TH1F*)c->GetPrimitive(hname1.Data());
    if (!h1) {
        printf("Error: Could not retrieve histogram '%s'. Check variable name.\n", hname1.Data());
        return;
    }

    h1->SetLineColor(kBlue + 1);
    h1->GetXaxis()->SetTitle(xAxisTitle.Data());

    // Draw second variable on the SAME canvas
    tree->Draw(var2_spec.Data(), cut2, "SAME");
    TH1F *h2 = (TH1F*)c->GetPrimitive(hname2.Data());
    if (!h2) {
        printf("Error: Could not retrieve histogram '%s'. Check variable name.\n", hname2.Data());
        return;
    }

    h2->SetLineColor(kRed + 1);


    // Auto-scale Y-axis
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    h1->SetMaximum(std::max(max1, max2) * 1.25);

  
    TLegend *leg = new TLegend(0.58, 0.68, 0.78, 0.88);
    leg->SetBorderSize(0);    // Removes the border line box
    leg->SetFillStyle(0);     // Makes the legend background transparent
    leg->SetTextSize(0.035);  
    leg->AddEntry(h1, label1.Data(), "l");
    leg->AddEntry(h2, label2.Data(), "l");
    leg->Draw();

    c->Update();

    // Save output
    if (!savePath.IsNull()) {
        c->SaveAs(savePath.Data());
    }
}