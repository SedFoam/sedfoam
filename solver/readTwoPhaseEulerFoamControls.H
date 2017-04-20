#include "readTimeControls.H"

int nAlphaCorr(readInt(pimple.dict().lookup("nAlphaCorr")));
Switch correctAlpha(pimple.dict().lookup("correctAlpha"));
