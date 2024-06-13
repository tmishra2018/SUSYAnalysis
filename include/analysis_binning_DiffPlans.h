#ifndef ANALYSIS_BINNING
#define ANALYSIS_BINNING

#ifndef ROOT_TFile
#include "TFile.h"
#endif
#include<iostream>

  class binning{
  public:
    binning(int NBIN, float METbin1, float METbin2, float METbin3, float HTbin1, float HTbin2, float HTbin3, float PHOETbin, float PHOETBin2):
			NBIN_(NBIN),
			METbin1_(METbin1), 
			METbin2_(METbin2), 
			METbin3_(METbin3), 
			HTbin1_(HTbin1),
			HTbin2_(HTbin2),
			HTbin3_(HTbin3),
			PHOETbin_(PHOETbin),
			PHOETBin2_(PHOETBin2)
	{
		}

	~binning(){
	};
// plan5
/*   inline int findSignalBin(float MET, float HT, float Et){

                int SigBinIndex(-1);
                int halfbin = NBIN_/2;
                if(NBIN_ > 20)halfbin = NBIN_/3;

                        if(MET > 120 && MET <= METbin1_){
				if(HT  < HTbin1_)SigBinIndex = 0;
                                else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 1; 
                                else if(HT  > HTbin2_ && HT < HTbin3_)SigBinIndex = 2;
                                else if(HT >= HTbin3_)SigBinIndex = 3;
                        }
                        else if(MET > METbin1_ && MET <= METbin2_){
				if(HT  < HTbin1_)SigBinIndex = 4;
                                else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 5; 
                                else if(HT  > HTbin2_ && HT < HTbin3_)SigBinIndex = 6;
                                else if(HT >= HTbin3_)SigBinIndex = 7;
                        }
                        else if(MET > METbin2_ && MET <= METbin3_){
				if(HT  < HTbin1_)SigBinIndex = 8;
                                else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 9; 
                                else if(HT  > HTbin2_ && HT < HTbin3_)SigBinIndex = 10;
                                else if(HT >= HTbin3_)SigBinIndex = 11;
                        }
                        else if(MET > METbin3_){ 
				if(HT  < HTbin1_)SigBinIndex = 12;
                                else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 13; 
                                else if(HT  > HTbin2_ && HT < HTbin3_)SigBinIndex = 14;
                                else if(HT >= HTbin3_)SigBinIndex = 15;

                        }
                        if(Et > PHOETbin_)SigBinIndex += 16;
                return SigBinIndex;
        }
*/
// plan4

       inline int findSignalBin(float MET, float HT, float Et){

                int SigBinIndex(-1);
                int halfbin = NBIN_/2;
                if(NBIN_ > 20)halfbin = NBIN_/3;
                if(HTbin1_ < HTbin2_-10){

                        if(MET > 120 && MET <= METbin1_){
                                if(HT  < HTbin1_)SigBinIndex = 0;
                                else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 1;
                                else if(HT >= HTbin2_)SigBinIndex = 2;
                        }
                        else if(MET > METbin1_ && MET <= METbin2_){
                                if(HT  < HTbin1_)SigBinIndex = 3;
                                else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 4;
                                else if(HT >= HTbin2_)SigBinIndex = 5;
                        }
                        else if(MET > METbin2_){
                                if(HT  < HTbin1_)SigBinIndex = 6;
                                else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 7;
                                else if(HT >= HTbin2_)SigBinIndex = 8;
                        }
                }
		if(SigBinIndex >= 0 && PHOETbin_ > 10){
                        if(Et > PHOETbin_)SigBinIndex += 9;
                }
                return SigBinIndex;
        }


/*
   	// plan1 : increasing MET bins 24 bins
	inline int findSignalBin(float MET, float HT, float Et){

		int SigBinIndex(-1);
		int halfbin = NBIN_/2;
		if(NBIN_ > 20)halfbin = NBIN_/3;

			if(MET > 120 && MET <= METbin1_){
				if(HT  < HTbin1_)SigBinIndex = 0;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 1;
				else if(HT >= HTbin2_)SigBinIndex = 2; 
			}
			else if(MET > METbin1_ && MET <= METbin2_){
				if(HT  < HTbin1_)SigBinIndex = 3;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 4; 
				else if(HT >= HTbin2_)SigBinIndex = 5;
			}
			else if(MET > METbin2_ && MET <= METbin3_){
				if(HT  < HTbin1_)SigBinIndex = 6;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 7; 
				else if(HT >= HTbin2_)SigBinIndex = 8;
			}
			else if(MET > METbin3_){ 
				if(HT  < HTbin1_)SigBinIndex = 9;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 10; 
				else if(HT >= HTbin2_)SigBinIndex = 11;
			}
			if(Et > PHOETbin_)SigBinIndex += 12;
		return SigBinIndex;
	}

	// Plan2
	inline int findSignalBin(float MET, float HT, float Et){

		int SigBinIndex(-1);
		int halfbin = NBIN_/2;
		if(NBIN_ > 20)halfbin = NBIN_/3;

			if(MET > 120 && MET <= METbin1_){
				if(HT  < HTbin1_)SigBinIndex = 0;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 1;
				else if(HT  > HTbin2_ && HT < HTbin3_)SigBinIndex = 2;
				else if(HT >= HTbin3_)SigBinIndex = 3; 
			}
			else if(MET > METbin1_ && MET <= METbin2_){
				if(HT  < HTbin1_)SigBinIndex = 4;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 5; 
				else if(HT  > HTbin2_ && HT < HTbin3_)SigBinIndex = 6;
				else if(HT >= HTbin3_)SigBinIndex = 7;
			}
			else if(MET > METbin2_){ 
				if(HT  < HTbin1_)SigBinIndex = 8;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 9; 
				else if(HT  > HTbin2_ && HT < HTbin3_)SigBinIndex = 10;
				else if(HT >= HTbin3_)SigBinIndex = 11;
			}
			if(Et > PHOETbin_)SigBinIndex += 12;
		return SigBinIndex;
	}





	//Plan 3  : 27 bins
	inline int findSignalBin(float MET, float HT, float Et){

		int SigBinIndex(-1);
		int halfbin;
		if(NBIN_ > 20)halfbin = NBIN_/3;

			if(MET > 120 && MET <= METbin1_){
				if(HT  < HTbin1_)SigBinIndex = 0;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 1;
				else if(HT >= HTbin2_)SigBinIndex = 2; 
			}
			else if(MET > METbin1_ && MET <= METbin2_){
				if(HT  < HTbin1_)SigBinIndex = 3;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 4; 
				else if(HT >= HTbin2_)SigBinIndex = 5;
			}
			else if(MET > METbin2_){ 
				if(HT  < HTbin1_)SigBinIndex = 6;
				else if(HT  > HTbin1_ && HT < HTbin2_)SigBinIndex = 7; 
				else if(HT >= HTbin2_)SigBinIndex = 8;
			}
			if(SigBinIndex >= 0 && Et > PHOETbin_ && Et < PHOETBin2_) SigBinIndex += 9;
			if(SigBinIndex >= 0 && Et >= PHOETBin2_) SigBinIndex += 18;
		return SigBinIndex;
	}


*/
	private:
		int	NBIN_;
		float	METbin1_;
		float	METbin2_;
		float	METbin3_;
		float	HTbin1_;
		float	HTbin2_;
		float	HTbin3_;
		float	PHOETbin_;
		float	PHOETBin2_;
	};

#endif
