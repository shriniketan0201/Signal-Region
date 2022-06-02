std::vector<int> postAnalyzer::getQcdden(double phoPtCut, double phoEtaCut, Short_t isoBit){
  std::vector<int> tmpCand;
  tmpCand.clear();
    
  for(int p=0;p<nPho;p++)    {
    Float_t uncorrectedPhoEt = ((*phoEt)[p]);
    bool passChIsoLoose = TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 6.63;
    bool passNeuIsoLoose = TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (15.57 + (0.0072* uncorrectedPhoEt) + (0.000017* pow(uncorrectedPhoEt, 2.0)));
    bool passPhoIsoLoose = TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (3.18 + (0.0019* uncorrectedPhoEt));
    
      //if(!passPFClusE || !passPFClusH || !passTrkH)
      if(!passChIsoLoose || !passNeuIsoLoose || !passPhoIsoLoose)	
	{
	  //	  bool kinematic = uncorrectedPhoEt > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
	  bool kinematic = uncorrectedPhoEt > phoPtCut  && fabs((*phoSCEta)[p])>1.566 && fabs((*phoSCEta)[p])<2.5;
	  //"Very loose" ID cuts with inverted shape cut
	  
	  bool photonID = (
                           //  ((*phoSigmaIEtaIEtaFull5x5)[p]  >  0.099 ) &&                                                                                             
                           ((*phoSigmaIEtaIEtaFull5x5)[p]  > 0.028 ) &&
                           ((*phoHoverE)[p]                <  0.0106   ) &&
                           ((*phohasPixelSeed)[p]              ==  0      ) &&
			   ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < TMath::Min(5.0*(6.63) , 0.20*uncorrectedPhoEt) )  &&
			   ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < TMath::Min(5.0*(15.57 + (0.0072 * uncorrectedPhoEt) + (0.000017* pow(uncorrectedPhoEt, 2.0))) , 0.20*uncorrectedPhoEt) )  &&
			   ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < TMath::Min(5.0*(3.18 + (0.0019 * uncorrectedPhoEt)) , 0.20*uncorrectedPhoEt) )

                           );


	  bool noncoll = fabs((*phoSeedTime)[p]) < 3.;// &&  (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001;

	  if(kinematic && photonID && noncoll)
	    {
	      tmpCand.push_back(p);
	    }
	}
    }
  return tmpCand;

}


std::vector<int> postAnalyzer::electron_veto_tightID(int pho_index, float elePtCut)
{
  std::vector<int> ele_cands;
  ele_cands.clear();
  for(int i = 0; i < nEle; i++)
    {
      if(eleIDbit->at(i)>>2&1==1)
        {
	  if(elePt->at(i) > elePtCut && abs(eleEta->at(i)) < 2.5 )
            {
	      if(dR(eleEta->at(i),elePhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5)
                {
                  ele_cands.push_back(i);
                }
            }
        }
    }
  return ele_cands;
}
std::vector<int> postAnalyzer::muon_veto_tightID(int pho_index, float muPtCut)
{
  std::vector<int> mu_cands;
  mu_cands.clear();
  for(int i = 0; i < nMu; i++)
    {
      if(muIDbit->at(i)>>3&1==1 && muIDbit->at(i)>>9&1==1)        
	{
	  if(muPt->at(i) > muPtCut && abs(muEta->at(i)) < 2.4 )        
	    {
              if(dR(muEta->at(i),muPhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5)
                {
                  mu_cands.push_back(i);
                }
            }
        }
    }
  return mu_cands;
}

std::vector<int> postAnalyzer::electron_veto_looseID(int pho_index, int mu_index, float elePtCut)
{
  std::vector<int> ele_cands;
  ele_cands.clear();
  for(int i = 0; i < nEle; i++)
    {
      if(eleIDbit->at(i)>>0&1==1)
	
	{
	  if(elePt->at(i) > elePtCut && abs(eleEta->at(i)) < 2.5 )
            {
	      
              if(dR(eleEta->at(i),elePhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5)
                {
                  ele_cands.push_back(i);
                }
            }
        }
    }
  return ele_cands;
}



std::vector<int> postAnalyzer::muon_veto_looseID(int pho_index, int ele_index, float muPtCut)
{
  std::vector<int> mu_cands;
  mu_cands.clear();
  for(int i = 0; i < nMu; i++)
    {
      if(muIDbit->at(i)>>0&1==1 && muIDbit->at(i)>>7&1==1)
	{
	  if(muPt->at(i) > muPtCut && abs(muEta->at(i)) < 2.4 )        
	    {
	      if(dR(muEta->at(i),muPhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5 && dR(muEta->at(i),muPhi->at(i),eleEta->at(ele_index),elePhi->at(ele_index)) > 0.5)
		{
		  mu_cands.push_back(i);
		}
	    }
	}
    }
  return mu_cands;
  
}
std::vector<int> postAnalyzer::JetVetoDecision(int pho_index, int ele_index) {

  bool jetVeto=true;
  std::vector<int> jetindex;
  float value =0.0;

  for(int i = 0; i < nJet; i++)
    {
      if(0.0 < abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.5)
        value =0.86;
      if(2.5 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.75)
        value =-0.10;
      if(2.75 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <3.0)
        value =-0.05;
      if(3.00 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <5.0)
        value =-0.01;
      
      double deltar = 0.0 ;
      double deltar_ele = 0.0;
      
      deltar= dR(jetEta->at(i),jetPhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index));
      deltar_ele = dR(jetEta->at(i),jetPhi->at(i),eleEta->at(ele_index),elePhi->at(ele_index));
      //std::cout<<"Delta R between photon and jet="<<dR<< "jet pt"<<pfJet_pt[i]<<std::endl; \
      
      if(deltar>0.4 && deltar_ele > 0.4 && jetPt->at(i) >30.0 && (jetID->at(i) >>1&1==1) && abs(jetEta->at(i))<5.0) 
        {
          jetindex.push_back(i);
        }

      
    }
  return jetindex;
  
}

bool postAnalyzer::OverlapWithMuon(double eta, double phi){
  
  return false;
  

}



bool postAnalyzer::OverlapWithElectron(double eta, double phi){
  return false;




}


double postAnalyzer::dR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

bool postAnalyzer::dPhiJetMET_veto(std::vector<int> jets)
{
  
  bool passes = false;
  int njetsMax = jets.size();
  
  if(njetsMax > 4)
    njetsMax = 4;
  int j = 0;
  for(; j < njetsMax; j++)
    {
      
      if(DeltaPhi(jetPhi->at(jets[j]),pfMETPhi) < 0.5)
        break;
    }
  if(j==njetsMax)
    passes = true;
  
  return passes;
}

double postAnalyzer::FakeRatePt(Float_t phoPt, int syst_number){
  //Fake rate with 2017_ID

  Float_t bin1[11] = {0.396309,0.378183,0.414314,0.415894,0.452537,0.388495,0.407523,0.396309,0.396309,0.408948,0.383671};
  Float_t bin2[11] = {0.310207,0.301206,0.306829,0.322366,0.322976,0.300976,0.318325,0.310207,0.310207,0.321899,0.298516};
  Float_t bin3[11] = {0.346647,0.337418,0.352949,0.344698,0.356272,0.318912,0.363794,0.346647,0.346647,0.369505,0.32379};
  Float_t bin4[11] = {0.221292,0.230239,0.0889309,0.227453,0.0696495,0.155854,0.304601,0.221292,0.221292,0.284499,0.158085};
  Float_t bin5[11] = {1.21786e-07,5.22148e-09,0.00628523,3.35331e-11,1.13427e-06,3.19878e-13,2.08163,1.21786e-07,1.21786e-07,0.282089,-0.282089};

  double weight = 1.0;
  if(phoPt <= 225.0)    weight = 0;
  else if(225.0 < phoPt && phoPt <= 250.0)    weight = bin1[syst_number];
  else if(250.0 < phoPt && phoPt <= 300.0)    weight = bin2[syst_number];
  else if(300.0 < phoPt && phoPt <= 400.0)    weight = bin3[syst_number];
  else if(400.0 < phoPt && phoPt <= 600.0)    weight = bin4[syst_number];
  else if(600.0 < phoPt )    weight = bin5[syst_number];
  return weight;
}
bool postAnalyzer::getMetFilter(){
  bool decision = true;
  if(metFilters>>0 & 1) decision = false;
  if(metFilters>>1 & 1) decision = false;
  if(metFilters>>2 & 1) decision = false;
  if(metFilters>>3 & 1) decision = false;
  if(metFilters>>4 & 1) decision = false;
  if(metFilters>>5 & 1) decision = false;
  if(metFilters>>7 & 1) decision = false;
  if(metFilters>>8 & 1) decision = false;
  return decision;
}

