## Using this will make things run faster...

numMaxProcesses = 4

## This will take all of the .root files in inDir
## And make friend trees in files in outDir
## If inDir and outDir are the same, then the files will be updated

inDir  = "/afs/cern.ch/work/d/dabercro/public/Winter15/flatTreesSkimmedV5/"
outDir = "/afs/cern.ch/work/d/dabercro/public/Winter15/flatTreesTest/"

## If the inTreeName and outTreeName are the same
## with inDir and outDir the same, the trees will just get an added branch
## Recommend backing up files before though

inTreeName  = "events"
outTreeName = "test"

## The inputs are TTreeFormula, in case you store Pt in a TLorentzVector #
## This is for the footprint correction

PhotonPtExpression  = "photonPtRaw"

OutputPhotonPt      = "photonFoot"
OutputPhotonSysUp   = "photonFootUp"
OutputPhotonSysDown = "photonFootDown"

## This is for the recoil smearing

GenBosonPtExpression    = "genBos_pt"
GenBosonPdgIdExpression = "genBos_PdgId"
EventNumExpression      = "eventNum"

OutputPerpName       = "uPerpSmeared"
OutputParaName       = "uParaSmeared"
OutputRecoilPtName   = "uMagSmeared"
