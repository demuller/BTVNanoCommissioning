import awkward as ak


def close_muon_tag(jets, muons, muon_dr=0.4, jet_dr=1.5, btag=0.0614, nminjets=2, nminmuons=1, jet_pt=30., muon_pt=5.):
    """_summary_

    Args:
        jets (PtEtaPhiMCandidate): Jet collection
        muons (PtEtaPhiMCandidate): Muon collections 
        muon_dr (float, optional): Maxiumum distance between muon to tag jet. Defaults to 0.4.
        jet_dr (float, optional): Minimum distance between tag to probe jet. Defaults to 1.4.
        btag (float, optional):  btag value for tag-jet. Defaults to 0.0614 (DeepJet loose WP Summer22EE).
        nminjets (int, optional):  Minimum number of jets per event. Defaults to 2.
        nminmuons (int, optional):  Minimum number of muons per event. Defaults to 1.
        jet_pt (float, optional): Minimum pt of jets. Defaults to 30.
        muon_pt (float, optional): Minimum pt of muon. Defaults to 5.

    Returns:
        PtEtaPhiMCandidate Pairs: Pairs of jet with fields "tag" and "probe". The probe jets
                                  have to be probed by the user afterwards 
                                  The tag-jet needs to have a muon close by
    """        
    jets["muon"] = jets.nearest(muons, threshold=muon_dr)
    # only work with jets that have a close muon
    jets = jets[~ak.is_none(jets.muon)]
    # kinematic and object cuts
    tnp_jet_mask =  (ak.num(jets.pt > jet_pt ) >= nminjets) & (jets.pt > jet_pt) & (abs(jets.eta) <= 2.5) & (jets.jetId >= 5)
    tnp_muon_mask = (ak.num(muons.pt) >= nminmuons) & (ak.any(muons.pt > muon_pt, axis=-1)) & (ak.any(abs(muons.eta) < 2.4, axis=-1))# & (ak.any(muons.tightId > 0.5, axis=-1)) & (ak.any(muons.pfRelIso04_all > 0.2, axis=-1))

    jets = jets[tnp_jet_mask & tnp_muon_mask]

    tnp_jet_pairs = ak.combinations(jets, 2, fields=["tag", "probe"])
    tnp_jet_pairs = tnp_jet_pairs[tnp_jet_pairs.tag.delta_r(tnp_jet_pairs.probe) > jet_dr]
    tnp_muon_mask = tnp_jet_pairs.tag.delta_r(tnp_jet_pairs.tag.muon) < muon_dr
    tnp_btag_mask = tnp_jet_pairs.tag.btagDeepFlavB > btag 
    tnp_jet_pairs = tnp_jet_pairs[tnp_muon_mask]
    tnp_jet_pairs = tnp_jet_pairs[tnp_btag_mask]
    #tnp_jet_pairs = tnp_jet_pairs[tnp_muon_mask & tnp_btag_mask]
    #tnp_jet_pairs = ak.flatten(tnp_jet_pairs)
    print("tag pt:",tnp_jet_pairs.tag.pt)
    print("probe pt:",tnp_jet_pairs.probe.pt)
    '''
    pnt_jet_pairs = ak.combinations(jets, 2, fields=["probe", "tag"])
    pnt_jet_pairs = pnt_jet_pairs[pnt_jet_pairs.tag.delta_r(pnt_jet_pairs.probe) > jet_dr]
    pnt_muon_mask = pnt_jet_pairs.tag.delta_r(pnt_jet_pairs.tag.muon) < muon_dr 
    pnt_btag_mask = pnt_jet_pairs.tag.btagDeepFlavB > btag
    pnt_jet_pairs = pnt_jet_pairs[pnt_muon_mask & pnt_btag_mask]
    #pnt_jet_pairs = ak.flatten(pnt_jet_pairs)
    print("tag pt 2:",tnp_jet_pairs.tag.pt)
    print("probe pt 2:",tnp_jet_pairs.probe.pt)

    all_jets = ak.concatenate((tnp_jet_pairs, pnt_jet_pairs), axis=0)
    print("all jets tag pt before ak any",all_jets.tag.pt)
    all_jets = all_jets[~ak.is_none(all_jets)]
    print(all_jets.tag.pt)
    print(all_jets.probe.pt)

    return all_jets
    '''
    return tnp_jet_pairs

def close_muon_probe(jets, muons, muon_dr=0.4, jet_dr=1.5, btag=0.0614, nminjets=2, nminmuons=1, jet_pt=30., muon_pt=5.):
    """_summary_

    Args:
        jets (PtEtaPhiMCandidate): Jet collection
        muons (PtEtaPhiMCandidate): Muon collections 
        muon_dr (float, optional): Maxiumum distance between muon to probe jet. Defaults to 0.4.
        jet_dr (float, optional): Minimum distance between tag to probe jet. Defaults to 1.4.
        btag (float, optional):  btag value for tag-jet. Defaults to 0.0614 (DeepJet loose WP Summer22EE).
        nminjets (int, optional):  Minimum number of jets per event. Defaults to 2.
        nminmuons (int, optional):  Minimum number of muons per event. Defaults to 1.
        jet_pt (float, optional): Minimum pt of jets. Defaults to 30.
        muon_pt (float, optional): Minimum pt of muon. Defaults to 5.

    Returns:
        PtEtaPhiMCandidate Pairs: Pairs of jet with fields "tag" and "probe". The probe jets
                                  have to be probed by the user afterwards.
                                  The probe-jet needs to have a muon close by
    """        
    jets["muon"] = jets.nearest(muons, threshold=muon_dr)
    # only work with jets that have a close muon
    jets = jets[~ak.is_none(jets.muon)]
    # kinnematic and object cuts
    tnp_jet_mask =  (ak.num(jets.pt) >= nminjets) & (ak.num(jets.pt > jet_pt ) >= nminjets) & (jets.pt > jet_pt)
    tnp_muon_mask = (ak.num(muons.pt) >= nminmuons) & (ak.any(muons.pt > muon_pt, axis=-1))
    #jets = jets[tnp_jet_mask & tnp_muon_mask]
    print("jet pt:",jets.pt)
    print("muon pt:",jets.muon.pt)
    tnp_jet_pairs = ak.combinations(jets, 2, fields=["tag", "probe"])
    tnp_jet_dR_mask = tnp_jet_pairs.tag.delta_r(tnp_jet_pairs.probe) > jet_dr
    #tnp_jet_pairs = tnp_jet_pairs[tnp_jet_pairs.tag.delta_r(tnp_jet_pairs.probe) > jet_dr]
    print(tnp_jet_pairs[3].to_list())
    tnp_muon_dR_mask = tnp_jet_pairs.probe.delta_r(tnp_jet_pairs.probe.muon) < muon_dr 
    print("tag jet pt",tnp_jet_pairs[2].tag.pt.to_list(), "probe muon pt",tnp_jet_pairs[2].probe.muon.pt.to_list())
    tnp_btag_mask = tnp_jet_pairs.tag.btagDeepFlavB > btag 
    print("tnp_muon_dR_mask & tnp_btag_mask",(tnp_muon_dR_mask & tnp_btag_mask)[:10])
    #tnp_jet_pairs = tnp_jet_pairs[tnp_btag_mask & tnp_muon_dR_mask]
    #tnp_jet_pairs = ak.flatten(tnp_jet_pairs)
    print("tag pt:",tnp_jet_pairs.tag.pt)
    print("probe pt:",tnp_jet_pairs.probe.pt)
    '''
    pnt_jet_pairs = ak.combinations(jets, 2, fields=["probe", "tag"])
    pnt_jet_pairs = pnt_jet_pairs[pnt_jet_pairs.tag.delta_r(pnt_jet_pairs.probe) > jet_dr]
    pnt_muon_mask = pnt_jet_pairs.probe.delta_r(pnt_jet_pairs.probe.muon) < muon_dr 
    pnt_btag_mask = pnt_jet_pairs.tag.btagDeepFlavB > btag
    pnt_jet_pairs = pnt_jet_pairs[pnt_muon_mask & pnt_btag_mask]
    pnt_jet_pairs = ak.flatten(pnt_jet_pairs)

    all_jets = ak.concatenate((tnp_jet_pairs, pnt_jet_pairs), axis=0)
    all_jets = all_jets[~ak.is_none(all_jets)]
    print("tag pt:",all_jets.tag.pt)
    print("probe pt:",all_jets.probe.pt)

    return all_jets
    '''
    print(ak.sum(tnp_jet_dR_mask & tnp_muon_dR_mask & tnp_btag_mask, axis=-1))
    print(ak.type(ak.sum(tnp_jet_mask,axis=-1)), ak.type(tnp_muon_mask))
    print(ak.type(ak.sum(tnp_jet_dR_mask & tnp_muon_dR_mask & tnp_btag_mask, axis=-1)))
    print(ak.sum(tnp_jet_dR_mask & tnp_muon_dR_mask & tnp_btag_mask, axis=-1))
    print(ak.sum(tnp_jet_mask,axis=-1))
    evt_lvl = ( (ak.sum(tnp_jet_mask, axis=-1) > 1) & tnp_muon_mask  & (ak.sum(tnp_jet_dR_mask & tnp_muon_dR_mask & tnp_btag_mask, axis=-1) == 1 ))
    '''
    tnp_jet_pairs.tag = tnp_jet_pairs.tag[ak.fill_none(((ak.sum(tnp_jet_mask, axis=-1) > 1) & tnp_muon_mask),False,axis=-1)]
    tnp_jet_pairs.probe = tnp_jet_pairs.probe[ak.fill_none(((ak.sum(tnp_jet_mask, axis=-1) > 1) & tnp_muon_mask),False,axis=-1)]

    tnp_jet_pairs = tnp_jet_pairs[ak.fill_none((tnp_jet_dR_mask & tnp_btag_mask & tnp_muon_dR_mask),False,axis=-1)] 
    tnp_jet_pairs = tnp_jet_pairs[ak.count(tnp_jet_pairs.tag.pt,axis=-1)==1] 
    print("here",ak.sum(evt_lvl),ak.type(tnp_jet_pairs)) 
    '''
    jets = jets[tnp_jet_mask] # & tnp_muon_mask] 
    jets = jets[ak.sum(tnp_jet_mask, axis=-1) > 1]
    tnp_jet_pairs = ak.combinations(jets, 2, fields=["tag", "probe"])
    #tnp_jet_pairs = tnp_jet_pairs[tnp_jet_pairs.tag.delta_r(tnp_jet_pairs.probe) > jet_dr]
    tnp_jet_dR_mask = tnp_jet_pairs.tag.delta_r(tnp_jet_pairs.probe) > jet_dr
    tnp_btag_mask = tnp_jet_pairs.tag.btagDeepFlavB > btag
    tnp_muon_dR_mask = tnp_jet_pairs.probe.delta_r(tnp_jet_pairs.probe.muon) < muon_dr
    #tnp_jet_pairs = tnp_jet_pairs[tnp_btag_mask & tnp_muon_dR_mask]
    tnp_jet_pairs = tnp_jet_pairs[ak.fill_none((tnp_jet_dR_mask & tnp_btag_mask & tnp_muon_dR_mask),False,axis=-1)] 
    tnp_jet_pairs = tnp_jet_pairs[ak.count(tnp_jet_pairs.tag.pt,axis=-1)==1] 

    print("here",ak.sum(evt_lvl),ak.type(tnp_jet_pairs))

    return tnp_jet_pairs, evt_lvl
