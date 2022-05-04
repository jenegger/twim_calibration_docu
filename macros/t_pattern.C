//insert all headers
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include <string>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"

int t_pattern(R3BEventHeader* tpat_tclone){
	Int_t tpatbin = 0;
	int trigger_pattern = 0;
	if (tpat_tclone->GetTpat() > 0){
		for (Int_t i = 0; i < 16; i++){
			tpatbin = (tpat_tclone->GetTpat() & (1 << i));
			if (tpatbin != 0){
				trigger_pattern = i + 1;
				}
		}
	}
return trigger_pattern;
}

uint64_t wrt_ts (R3BEventHeader* tpat_tclone){
	uint64_t timestamp = tpat_tclone->GetTimeStamp();
	return timestamp;
	}
