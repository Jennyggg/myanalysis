import uproot
import pandas as pd
import numpy as np
filesbkg=['/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_dyjetsll.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_gjets.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_qcd.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_singletop.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_tt.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_tt_negligible.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_wjets.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_ww_wz.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2016_zinv.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_dyjetsll.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_gjets.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_qcd.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_singletop.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_tt.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_tt_negligible.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_wjets.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_ww_wz.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2017_zinv.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_dyjetsll.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_gjets.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_qcd.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_singletop.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_tt.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_tt_negligible.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_wjets.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_ww_wz.root',
'/scratch/wjin/featurereduced2/rootfiles/bkg/preselection_2018_zinv.root']

filessig=['/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T1bbbb1.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T1bbbb2.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T1bbbb3.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T1qqqq1.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T1qqqq2.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T1qqqq3.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T2bb1.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T2bb2.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T2bb3.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T2qq1.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T2qq2.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2016_T2qq3.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1bbbb1.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1bbbb2.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1bbbb3.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1bbbb4.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1bbbb5.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1bbbb6.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1bbbb7.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1bbbb8.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1qqqq1.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1qqqq2.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1qqqq3.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1qqqq4.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1qqqq5.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1qqqq6.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1qqqq7.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T1qqqq8.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2bb1.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2bb2.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2bb3.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2bb4.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2bb5.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2bb6.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2bb7.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq1.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq2.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq3.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq4.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq5.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq6.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq7.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq8.root',
'/scratch/wjin/featurereduced2/rootfiles/sig/preselection_2017_T2qq9.root']


for fn in filesbkg:
    tree=uproot.open(fn)["mt2"]
    arrays=tree.arrays()
    print("load samples from"+fn)
    np.save('/scratch/wjin/featurereduced3/bkg/'
            +fn.replace("/scratch/wjin/featurereduced2/rootfiles/bkg/",'').replace(".root",'')+'_slim.npy', arrays)


for fn in filessig:
    tree=uproot.open(fn)["mt2"]
    arrays=tree.arrays()
    print("load samples from"+fn)
    np.save('/scratch/wjin/featurereduced3/sig/'
            +fn.replace("/scratch/wjin/featurereduced2/rootfiles/sig/",'').replace(".root",'')+'_slim.npy', arrays)
