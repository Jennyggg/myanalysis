import uproot
import pandas as pd
import numpy as np

files2016mcbkg_dyjetsll=["/scratch/mratti/NEWSnTtrees/2016/dyjetsll_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2016/dyjetsll_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2016/dyjetsll_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2016/dyjetsll_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2016/dyjetsll_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2016/dyjetsll_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2016/dyjetsll_ht800to1200.root"]
#/scratch/mratti/NEWSnTtrees/2016/dyjetsll_m50_incl.root

files2016mcbkg_gjets=["/scratch/mratti/NEWSnTtrees/2016/gjets_dr0p05_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2016/gjets_dr0p05_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2016/gjets_dr0p05_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2016/gjets_dr0p05_ht40to100.root",
"/scratch/mratti/NEWSnTtrees/2016/gjets_dr0p05_ht600toInf.root"]

files2016mcbkg_qcd=["/scratch/mratti/NEWSnTtrees/2016/qcd_ht1000to1500.root",
"/scratch/mratti/NEWSnTtrees/2016/qcd_ht1500to2000.root",
"/scratch/mratti/NEWSnTtrees/2016/qcd_ht2000toInf.root",
"/scratch/mratti/NEWSnTtrees/2016/qcd_ht300to500.root",
"/scratch/mratti/NEWSnTtrees/2016/qcd_ht500to700.root",
"/scratch/mratti/NEWSnTtrees/2016/qcd_ht700to1000.root",
"/scratch/mratti/NEWSnTtrees/2016/singletop_schan.root"]

files2016mcbkg_singletop=["/scratch/mratti/NEWSnTtrees/2016/singletop_tW_tbar.root",
"/scratch/mratti/NEWSnTtrees/2016/singletop_tW_top.root",
"/scratch/mratti/NEWSnTtrees/2016/singletop_tchan_tbar.root",
"/scratch/mratti/NEWSnTtrees/2016/singletop_tchan_top.root"]

files2016mcbkg_tt=["/scratch/mratti/NEWSnTtrees/2016/ttdl.root",
"/scratch/mratti/NEWSnTtrees/2016/ttsl_tbar.root",
"/scratch/mratti/NEWSnTtrees/2016/ttsl_top.root",
"/scratch/mratti/NEWSnTtrees/2016/ttz_mg_lo.root"]

files2016mcbkg_tt_negligible=["/scratch/mratti/NEWSnTtrees/2016/ttg_amcatnlo.root",
"/scratch/mratti/NEWSnTtrees/2016/tttt.root",
"/scratch/mratti/NEWSnTtrees/2016/ttw_lnu_amcatnlo.root",
"/scratch/mratti/NEWSnTtrees/2016/ttw_qq_amcatnlo.root"]

files2016mcbkg_wjets=["/scratch/mratti/NEWSnTtrees/2016/wjets_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2016/wjets_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2016/wjets_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2016/wjets_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2016/wjets_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2016/wjets_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2016/wjets_ht800to1200.root"]

#/scratch/mratti/NEWSnTtrees/2016/wjets_incl.root

files2016mcbkg_ww_wz=["/scratch/mratti/NEWSnTtrees/2016/ww.root",
"/scratch/mratti/NEWSnTtrees/2016/wz.root"]

files2016mcbkg_zinv=["/scratch/mratti/NEWSnTtrees/2016/zinv_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2016/zinv_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2016/zinv_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2016/zinv_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2016/zinv_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2016/zinv_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2016/zinv_ht800to1200.root"]

files2017mcbkg_dyjetsll=["/scratch/mratti/NEWSnTtrees/2017/dyjetsll_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2017/dyjetsll_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2017/dyjetsll_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2017/dyjetsll_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2017/dyjetsll_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2017/dyjetsll_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2017/dyjetsll_ht800to1200.root"]

#/scratch/mratti/NEWSnTtrees/2017/dyjetsll_m50_incl.root
files2017mcbkg_gjets=["/scratch/mratti/NEWSnTtrees/2017/gjets_dr0p05_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2017/gjets_dr0p05_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2017/gjets_dr0p05_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2017/gjets_dr0p05_ht40to100.root",
"/scratch/mratti/NEWSnTtrees/2017/gjets_dr0p05_ht600toInf.root"]

files2017mcbkg_qcd=["/scratch/mratti/NEWSnTtrees/2017/qcd_ht1000to1500.root",
"/scratch/mratti/NEWSnTtrees/2017/qcd_ht1500to2000.root",
"/scratch/mratti/NEWSnTtrees/2017/qcd_ht2000toInf.root",
"/scratch/mratti/NEWSnTtrees/2017/qcd_ht300to500.root",
"/scratch/mratti/NEWSnTtrees/2017/qcd_ht500to700.root",
"/scratch/mratti/NEWSnTtrees/2017/qcd_ht700to1000.root"]

files2017mcbkg_singletop=["/scratch/mratti/NEWSnTtrees/2017/singletop_schan.root",
"/scratch/mratti/NEWSnTtrees/2017/singletop_tW_tbar.root",
"/scratch/mratti/NEWSnTtrees/2017/singletop_tW_top.root",
"/scratch/mratti/NEWSnTtrees/2017/singletop_tchan_tbar.root",
"/scratch/mratti/NEWSnTtrees/2017/singletop_tchan_top.root"]

files2017mcbkg_tt=["/scratch/mratti/NEWSnTtrees/2017/ttdl_mg.root",
"/scratch/mratti/NEWSnTtrees/2017/ttsl_tbar_mg.root",
"/scratch/mratti/NEWSnTtrees/2017/ttsl_top_mg.root",
"/scratch/mratti/NEWSnTtrees/2017/ttz_mg.root"]

files2017mcbkg_tt_negligible=["/scratch/mratti/NEWSnTtrees/2017/ttg_amcatnlo.root",
"/scratch/mratti/NEWSnTtrees/2017/tttt.root",
"/scratch/mratti/NEWSnTtrees/2017/ttw_mg.root"]

files2017mcbkg_wjets=["/scratch/mratti/NEWSnTtrees/2017/wjets_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2017/wjets_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2017/wjets_ht1200to2500_1.root",
"/scratch/mratti/NEWSnTtrees/2017/wjets_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2017/wjets_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2017/wjets_ht2500toInf_1.root",
"/scratch/mratti/NEWSnTtrees/2017/wjets_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2017/wjets_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2017/wjets_ht800to1200.root"]
#/scratch/mratti/NEWSnTtrees/2017/wjets_incl.root
files2017mcbkg_ww_wz=["/scratch/mratti/NEWSnTtrees/2017/ww.root",
"/scratch/mratti/NEWSnTtrees/2017/wz.root"]

files2017mcbkg_zinv=["/scratch/mratti/NEWSnTtrees/2017/zinv_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2017/zinv_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2017/zinv_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2017/zinv_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2017/zinv_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2017/zinv_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2017/zinv_ht800to1200.root"]

files2018mcbkg_dyjetsll=["/scratch/mratti/NEWSnTtrees/2018/dyjetsll_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2018/dyjetsll_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2018/dyjetsll_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2018/dyjetsll_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2018/dyjetsll_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2018/dyjetsll_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2018/dyjetsll_ht800to1200.root"]
#/scratch/mratti/NEWSnTtrees/2018/dyjetsll_m50_incl.root

files2018mcbkg_gjets=["/scratch/mratti/NEWSnTtrees/2018/gjets_dr0p05_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2018/gjets_dr0p05_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2018/gjets_dr0p05_ht40to100.root",
"/scratch/mratti/NEWSnTtrees/2018/gjets_dr0p05_ht600toInf_ext1.root"]

files2018mcbkg_qcd=["/scratch/mratti/NEWSnTtrees/2018/qcd_ht1000to1500.root",
"/scratch/mratti/NEWSnTtrees/2018/qcd_ht1500to2000.root",
"/scratch/mratti/NEWSnTtrees/2018/qcd_ht2000toInf.root",
"/scratch/mratti/NEWSnTtrees/2018/qcd_ht300to500.root",
"/scratch/mratti/NEWSnTtrees/2018/qcd_ht500to700.root",
"/scratch/mratti/NEWSnTtrees/2018/qcd_ht700to1000.root"]

files2018mcbkg_singletop=["/scratch/mratti/NEWSnTtrees/2018/singletop_schan.root",
"/scratch/mratti/NEWSnTtrees/2018/singletop_tW_tbar.root",
"/scratch/mratti/NEWSnTtrees/2018/singletop_tW_top.root",
"/scratch/mratti/NEWSnTtrees/2018/singletop_tchan_tbar.root",
"/scratch/mratti/NEWSnTtrees/2018/singletop_tchan_top.root"]

files2018mcbkg_tt=["/scratch/mratti/NEWSnTtrees/2018/ttdl_mg.root",
"/scratch/mratti/NEWSnTtrees/2018/ttsl_tbar_mg.root",
"/scratch/mratti/NEWSnTtrees/2018/ttsl_top_mg.root",
"/scratch/mratti/NEWSnTtrees/2018/ttz_mg.root"]

files2018mcbkg_tt_negligible=["/scratch/mratti/NEWSnTtrees/2018/tttt.root",
"/scratch/mratti/NEWSnTtrees/2018/ttw_mg.root"]

files2018mcbkg_wjets=["/scratch/mratti/NEWSnTtrees/2018/wjets_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2018/wjets_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2018/wjets_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2018/wjets_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2018/wjets_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2018/wjets_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2018/wjets_ht800to1200.root"]
#/scratch/mratti/NEWSnTtrees/2018/wjets_incl.root

files2018mcbkg_ww_wz=["/scratch/mratti/NEWSnTtrees/2018/ww.root",
"/scratch/mratti/NEWSnTtrees/2018/wz.root"]

files2018mcbkg_zinv=["/scratch/mratti/NEWSnTtrees/2018/zinv_ht100to200.root",
"/scratch/mratti/NEWSnTtrees/2018/zinv_ht1200to2500.root",
"/scratch/mratti/NEWSnTtrees/2018/zinv_ht200to400.root",
"/scratch/mratti/NEWSnTtrees/2018/zinv_ht2500toInf.root",
"/scratch/mratti/NEWSnTtrees/2018/zinv_ht400to600.root",
"/scratch/mratti/NEWSnTtrees/2018/zinv_ht600to800.root",
"/scratch/mratti/NEWSnTtrees/2018/zinv_ht800to1200.root"]

files2016signalT1bbbb=["/scratch/mratti/NEWSnTtrees/2016/signal/T1bbbb_94x_1.root",
"/scratch/mratti/NEWSnTtrees/2016/signal/T1bbbb_94x_2.root",
"/scratch/mratti/NEWSnTtrees/2016/signal/T1bbbb_94x.root"]

files2016signalT1qqqq=["/scratch/mratti/NEWSnTtrees/2016/signal/T1qqqq_94x_1.root",
"/scratch/mratti/NEWSnTtrees/2016/signal/T1qqqq_94x_2.root",
"/scratch/mratti/NEWSnTtrees/2016/signal/T1qqqq_94x.root"]

files2016signalT2bb=["/scratch/mratti/NEWSnTtrees/2016/signal/T2bb_1.root",
"/scratch/mratti/NEWSnTtrees/2016/signal/T2bb_mSbot1650to2600.root",
"/scratch/mratti/NEWSnTtrees/2016/signal/T2bb.root"]

files2016signalT2qq=["/scratch/mratti/NEWSnTtrees/2016/signal/T2qq_1.root",
"/scratch/mratti/NEWSnTtrees/2016/signal/T2qq_mSq1850to2600.root",
"/scratch/mratti/NEWSnTtrees/2016/signal/T2qq.root"]


files2017signalT1bbbb=["/scratch/mratti/NEWSnTtrees/2017/signal/T1bbbb_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1bbbb_2.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1bbbb_ext1_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1bbbb_ext1_2.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1bbbb_ext1_3.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1bbbb_ext1_4.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1bbbb_ext1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1bbbb.root"]

files2017signalT1qqqq=["/scratch/mratti/NEWSnTtrees/2017/signal/T1qqqq_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1qqqq_2.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1qqqq_ext1_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1qqqq_ext1_2.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1qqqq_ext1_3.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1qqqq_ext1_4.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1qqqq_ext1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T1qqqq.root"]

files2017signalT2bb=["/scratch/mratti/NEWSnTtrees/2017/signal/T2bb_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2bb_ext1_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2bb_ext1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2bb_mSbot1650to2600_ext1_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2bb_mSbot1650to2600_ext1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2bb_mSbot1650to2600.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2bb.root"]

files2017signalT2qq=["/scratch/mratti/NEWSnTtrees/2017/signal/T2qq_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2qq_2.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2qq_ext1_1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2qq_ext1_2.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2qq_ext1_3.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2qq_ext1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2qq_mSq_1850to2600_ext1.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2qq_mSq_1850to2600.root",
"/scratch/mratti/NEWSnTtrees/2017/signal/T2qq.root"]

filessignal=files2016signalT1bbbb+files2016signalT1qqqq+files2016signalT2bb+files2016signalT2qq \
+files2017signalT1bbbb+files2017signalT1qqqq+files2017signalT2bb+files2017signalT2qq


year2016=[files2016mcbkg_dyjetsll,files2016mcbkg_gjets,files2016mcbkg_qcd,files2016mcbkg_singletop, \
files2016mcbkg_tt,files2016mcbkg_tt_negligible,files2016mcbkg_wjets,files2016mcbkg_ww_wz,files2016mcbkg_zinv, \
files2016signalT1bbbb,files2016signalT1qqqq, \
files2016signalT2bb,files2016signalT2qq]

year2017=[files2017mcbkg_dyjetsll,files2017mcbkg_gjets,files2017mcbkg_qcd,files2017mcbkg_singletop, \
files2017mcbkg_tt,files2017mcbkg_tt_negligible,files2017mcbkg_wjets,files2017mcbkg_ww_wz,files2017mcbkg_zinv, \
files2017signalT1bbbb,files2017signalT1qqqq,files2017signalT2bb,files2017signalT2qq]

bkgfile=[files2016mcbkg_dyjetsll,files2016mcbkg_gjets,files2016mcbkg_qcd,files2016mcbkg_singletop, \
files2016mcbkg_tt,files2016mcbkg_tt_negligible,files2016mcbkg_wjets,files2016mcbkg_ww_wz,files2016mcbkg_zinv, \
files2017mcbkg_dyjetsll,files2017mcbkg_gjets,files2017mcbkg_qcd,files2017mcbkg_singletop, \
files2017mcbkg_tt,files2017mcbkg_tt_negligible,files2017mcbkg_wjets,files2017mcbkg_ww_wz,files2017mcbkg_zinv, \
files2018mcbkg_dyjetsll,files2018mcbkg_gjets,files2018mcbkg_qcd,files2018mcbkg_singletop, \
files2018mcbkg_tt,files2018mcbkg_tt_negligible,files2018mcbkg_wjets,files2018mcbkg_ww_wz,files2018mcbkg_zinv]

features_preselection1=["nJet30","nJet40","nBJet20","nBJet30","nBJet40", \
                        "deltaPhiMin","diffMetMht","ht","mt2", \
                        "jet1_pt","jet2_pt","mht_pt","mht_phi","met_pt","met_phi", \
                        "GenSusyMScan1","GenSusyMScan2", \
                        "njet","jet_pt","jet_eta","jet_phi","jet_mass","jet_btagDeepCSV","evt_scale1fb"]
#features_preselection2=["jet_chf","jet_nhf","jet_cemf","jet_nemf", \
#                        "ntau","tau_pt","tau_eta","tau_phi","tau_mass","tau_charge","tau_dxy","tau_dz"]
features_preselection3=['nMuons10','nElectrons10','nPFLep5LowMT','nPFHad10LowMT']
filecount=0
for filelist in bkgfile:
    m=0
    arraysfromlist={}
    for fn in filelist:
        tree=uproot.open(fn)["mt2"]
        arrays=tree.arrays(features_preselection1+features_preselection3)
        selectindex=(arrays[b'njet']>0) & (arrays[b'ht']>250) & \
        (((arrays[b'ht']<1200) & (arrays[b'met_pt']>250) )| ((arrays[b'ht']>=1200) & (arrays[b'met_pt']>30))) \
        & (((arrays[b'ht']<1500) & (arrays[b'mt2']>200)) | ((arrays[b'ht']>=1500) & (arrays[b'mt2']>400))) \
        & (arrays[b'deltaPhiMin']>0.3) \
        & ((arrays[b'mht_pt']<(1.5*arrays[b'met_pt'])) | (arrays[b'mht_pt']>(0.5*arrays[b'met_pt']))) \
        & (arrays[b'nMuons10']==0) & (arrays[b'nElectrons10']==0) \
        & (arrays[b'nPFLep5LowMT']==0) & (arrays[b'nPFHad10LowMT']==0)
        for i in arrays.keys():
            if m:
                try:
                    arraysfromlist[i]=np.concatenate((arrays[i][selectindex],arraysfromlist[i]),axis=0)
                except ValueError:
                    arraysselectindex=np.ndarray(shape=(arrays[i][selectindex].size,), dtype="object")
                    for mm in range(arrays[i][selectindex].size):
                        arraysselectindex[mm]=np.array([],dtype="float32")
                    arraysfromlist[i]=np.concatenate((arraysselectindex,arraysfromlist[i]),axis=0)
            else:
                arraysfromlist[i]=arrays[i][selectindex]
        m=m+1
        print("load samples from"+fn)
    if filelist in year2016:
        year='2016'
    elif filelist in year2017:
        year='2017'
    else:
        year='2018'

    if filelist in [files2016mcbkg_dyjetsll,files2017mcbkg_dyjetsll,files2018mcbkg_dyjetsll]:
        cati='dyjetsll'
    elif filelist in [files2016mcbkg_gjets,files2017mcbkg_gjets,files2018mcbkg_gjets]:
        cati='gjets'
    elif filelist in [files2016mcbkg_qcd,files2017mcbkg_qcd,files2018mcbkg_qcd]:
        cati='qcd'
    elif filelist in [files2016mcbkg_singletop,files2017mcbkg_singletop,files2018mcbkg_singletop]:
        cati='singletop'
    elif filelist in [files2016mcbkg_tt,files2017mcbkg_tt,files2018mcbkg_tt]:
        cati='tt'
    elif filelist in [files2016mcbkg_tt_negligible,files2017mcbkg_tt_negligible,files2018mcbkg_tt_negligible]:
        cati='tt_negligible'
    elif filelist in [files2016mcbkg_wjets,files2017mcbkg_wjets,files2018mcbkg_wjets]:
        cati='wjets'
    elif filelist in [files2016mcbkg_ww_wz,files2017mcbkg_ww_wz,files2018mcbkg_ww_wz]:
        cati='ww_wz'
    elif filelist in [files2017mcbkg_zinv,files2017mcbkg_zinv,files2018mcbkg_zinv]:
        cati='zinv'
#    else:
#        if filelist in T1bbbbfile:
#            cati='sig_T1bbbb'
#        elif filelist in T1qqqqfile:
#            cati='sig_T1qqqq'
#        elif filelist in T2bbfile:
#            cati='sig_T2bb'
#        else:
 #           cati='sig_T2qq'

    del arraysfromlist[b'nMuons10']
    del arraysfromlist[b'nElectrons10']
    del arraysfromlist[b'nPFLep5LowMT']
    del arraysfromlist[b'nPFHad10LowMT']
    np.save('/scratch/wjin/featurereduced2/bkg/preselection_'+year+'_'+cati+str(filecount)+'_slim.npy', arraysfromlist)
    print("produced file /scratch/wjin/featurereduced1/bkg/preselection_"+year+'_'+cati+str(filecount)+"_slim.npy")
    filecount=filecount+1



for fn in filessignal:
    arraysfromlist={}
    tree=uproot.open(fn)["mt2"]
    arrays=tree.arrays(features_preselection1+features_preselection3)
    selectindex=(arrays[b'njet']>0) & (arrays[b'ht']>250) & \
    (((arrays[b'ht']<1200) & (arrays[b'met_pt']>250) )| ((arrays[b'ht']>=1200) & (arrays[b'met_pt']>30))) \
    & (((arrays[b'ht']<1500) & (arrays[b'mt2']>200)) | ((arrays[b'ht']>=1500) & (arrays[b'mt2']>400))) \
    & (arrays[b'deltaPhiMin']>0.3) \
    & ((arrays[b'mht_pt']<(1.5*arrays[b'met_pt'])) | (arrays[b'mht_pt']>(0.5*arrays[b'met_pt']))) \
    & (arrays[b'nMuons10']==0) & (arrays[b'nElectrons10']==0) \
    & (arrays[b'nPFLep5LowMT']==0) & (arrays[b'nPFHad10LowMT']==0)
    for i in arrays.keys():
        arraysfromlist[i]=arrays[i][selectindex]
    print("load samples from"+fn)
    del arraysfromlist[b'nMuons10']
    del arraysfromlist[b'nElectrons10']
    del arraysfromlist[b'nPFLep5LowMT']
    del arraysfromlist[b'nPFHad10LowMT']
    np.save('/scratch/wjin/featurereduced2/'+fn.replace("/scratch/mratti/NEWSnTtrees/",'').replace(".root",'')+'_slim.npy', arraysfromlist)
    print("produced file /scratch/wjin/featurereduced2/"+fn.replace("/scratch/mratti/NEWSnTtrees/",'').replace(".root",'')+"_slim.npy")

