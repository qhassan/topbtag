import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000))

from Topbjets.Eff.btag.TTJets_cfi import source as events_source
#from test_fi.Testfi.btag.SingleMu_2017B_cfi import source as events_source

process.source=events_source

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#tfileservice
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("TTJets.root")
#        fileName = cms.string("SingleMu_2017B.root")
)

#running sequence
#process.load('UserCode.TopAnalysis.myChargedPFJets_cfi')
#process.p = cms.Path(process.myChargedPFJets*process.demo)
process.load('Topbjets.Eff.miniAnalyzer_cfi')
process.p = cms.Path(process.demo)


