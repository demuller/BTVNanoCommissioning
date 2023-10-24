import collections, awkward as ak, numpy as np
import os
import uproot
from coffea import processor
from coffea.analysis_tools import Weights

# functions to load SFs, corrections
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    eleSFs,
    muSFs,
    puwei,
    btagSFs,
    JME_shifts,
    Roccor_shifts,
)

# user helper function
from BTVNanoCommissioning.helpers.func import (
    flatten,
    update,
    uproot_writeable,
    dump_lumi,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch

from BTVNanoCommissioning.helpers.TnP_helper import (
    close_muon_tag,
    close_muon_probe,
)



## load histograms & selctions for this workflow
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, mu_idiso, ele_cuttightid


class NanoProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        ## Load corrections
        self.SF_map = load_SF(self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    ## Apply corrections on momentum/mass on MET, Jet, Muon
    def process(self, events):
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        events = missing_branch(events)
        shifts = []
        if "JME" in self.SF_map.keys():
            syst_JERC = True if self.isSyst != None else False
            if self.isSyst == "JERC_split":
                syst_JERC = "split"  # JEC splitted into 11 sources instead of JES_total
            shifts = JME_shifts(
                shifts, self.SF_map, events, self._campaign, isRealData, syst_JERC
            )
        else:
            if "Run3" not in self._campaign:
                shifts = [
                    ({"Jet": events.Jet, "MET": events.MET, "Muon": events.Muon}, None)
                ]
            else:
                shifts = [
                    (
                        {
                            "Jet": events.Jet,
                            "MET": events.PuppiMET,
                            "Muon": events.Muon,
                        },
                        None,
                    )
                ]
        if "roccor" in self.SF_map.keys():
            shifts = Roccor_shifts(shifts, self.SF_map, events, isRealData, False)
        else:
            shifts[0][0]["Muon"] = events.Muon

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    ## Processed events per-chunk, made selections, filled histogram, stored root files
    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        ######################
        #  Create histogram  # : Get the histogram dict from `histogrammer`
        ######################
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, "TnP")
        )

        output = {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)
        ####################
        #    Selections    #
        ####################
        ## Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        # only dump for nominal case
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        ## HLT
        triggers = [
            "BTagMu_AK4DiJet110_Mu5",
            "BTagMu_AK4DiJet170_Mu5",
            "BTagMu_AK4DiJet20_Mu5",
            "BTagMu_AK4DiJet40_Mu5",
            "BTagMu_AK4DiJet70_Mu5",
            "BTagMu_AK4Jet300_Mu5"
        ]
        checkHLT = ak.Array([hasattr(events.HLT, _trig) for _trig in triggers])
        if ak.all(checkHLT == False):
            raise ValueError("HLT paths:", triggers, " are all invalid in", dataset)
        elif ak.any(checkHLT == False):
            print(np.array(triggers)[~checkHLT], " not exist in", dataset)
        trig_arrs = [
            events.HLT[_trig] for _trig in triggers if hasattr(events.HLT, _trig)
        ]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ##### Add some selections
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2

        ## Jet cuts

        ## Other cuts

        ## Apply all selections
        event_level = (
            req_trig & req_lumi  # & req_muon & req_ele & req_jets & req_opposite_charge
        )
        #event_level = ak.fill_none(event_level, False)
        ## Skip empty events
        #if len(events[event_level]) == 0:
        #    return {dataset: output}

        ####################
        # Selected objects # : Pruned objects with reduced event_level
        ####################
        tnp = close_muon_probe(events.Jet[event_level],events.Muon[event_level])
        sjets_tag = tnp[0].tag
        #sjets_tag = close_muon_probe(sjets,smu).tag
        print(ak.type(sjets_tag))
        sjets_probe = tnp[0].probe
        evt_lvl = tnp[1]
        print(evt_lvl)
        event_level = evt_lvl

        event_level = ak.fill_none(event_level, False)
        # Skip empty events
        if len(events[event_level]) == 0:
            return {dataset: output}

        smu = tnp[0].probe.muon
        #sjets_probe = close_muon_probe(sjets,smu).probe
        print("n tag jets:",ak.count(sjets_tag.pt,axis=1))
        print("n probe jets:",ak.count(sjets_probe.pt,axis=1))
        print(len(sjets_tag.pt),len(sjets_probe.pt))
        print(len(events[event_level]))
        print("sum n tags / probes",sum(ak.count(sjets_tag.pt,axis=1)),sum(ak.count(sjets_probe.pt,axis=1)))
        print(ak.flatten(sjets_tag.pt),len(ak.flatten(sjets_tag.pt)))
        print(ak.flatten(sjets_probe.pt), len(ak.flatten(sjets_probe.pt)))
        sjets = ak.concatenate([sjets_tag, sjets_probe])
        #sjets = np.vstack([sjets_tag,sjets_probe])
        print(ak.type(sjets.pt))
        for b in sjets_tag.fields:
            print(b, np.shape(sjets[b].to_numpy()), ak.type(sjets[b]))

        sjets= ak.zip(
        {
        b : ak.Array(np.reshape(sjets[b].to_numpy(), (len(sjets_tag[b]), 2)))
        for b in sjets_tag.fields if "muon" not in b and "IdxG" not in b
        }
        ) 
        print(ak.type(sjets))
        ####################
        # Weight & Geninfo # : Add weight to selected events
        ####################
        # create Weights object to save individual weights
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            # Load SFs
            if len(self.SF_map.keys()) > 0:
                syst_wei = (
                    True if self.isSyst != None else False
                )  # load systematic flag
                if "PU" in self.SF_map.keys():
                    puwei(
                        events[event_level].Pileup.nTrueInt,
                        self.SF_map,
                        weights,
                        syst_wei,
                    )
                if "MUO" in self.SF_map.keys():
                    muSFs(
                        smu, self.SF_map, weights, syst_wei, False
                    )  # input selected muon
                if "EGM" in self.SF_map.keys():
                    eleSFs(sele, self.SF_map, weights, syst_wei, False)
                if "BTV" in self.SF_map.keys():
                    # For BTV weight, you need to specify type
                    btagSFs(sjets, self.SF_map, weights, "DeepJetC", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepJetB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVC", syst_wei)
        else:
            genflavor = ak.zeros_like(sjets.pt)
        print("gen:",ak.any(ak.count(genflavor, axis=-1) == 0))
        # Systematics information (add name of systematics)
        if shift_name is None:  # weight variations
            systematics = ["nominal"] + list(weights.variations)
        else:  # resolution/ scale variation would use the shift_name
            systematics = [shift_name]
        exclude_btv = [
            "DeepCSVC",
            "DeepCSVB",
            "DeepJetB",
            "DeepJetB",
        ]  # exclude b-tag SFs for btag inputs

        ####################
        #  Fill histogram  #
        ####################
        for syst in systematics:
            if self.isSyst == None and syst != "nominal":
                break
            if self.noHist:
                break
            weight = (
                weights.weight()
                if syst == "nominal" or syst == shift_name
                else weights.weight(modifier=syst)
            )  # shift up/down for weight systematics
            print("hello")
            for histname, h in output.items():
                print("histname:",histname)
                if "jet_tag" in histname and "dr" not in histname and "njet" != histname:
                    sel_tag = sjets_tag
                    print("sjets_tag:",sjets_tag)
                    print("hist:",histname)
                    print(ak.type(genflavor))
                    print(genflavor[0])
                    print("len gen:",len(flatten(genflavor[:,0])),"len sel_tag:",len(flatten(sel_tag[histname.replace(f"jet_tag_","")])),len(weight))
                    h.fill(syst, flatten(genflavor[:,0]), flatten(sel_tag[histname.replace(f"jet_tag_","")]), weight=weight)
            #output["jet_tag_pt"].fill(syst, seljet_tag.pt, weight=weight)
            #output["jet_probe_pt"].fill(syst, seljet_probe.pt, weight=weight)
            output["njet"].fill(syst, ak.count(sjets_tag.pt,axis=1), weight=weight)


        #######################
        #  Create root files  # : Save arrays in to root file, keep axis structure
        #######################
        if self.isArray:
            # Keep the structure of events and pruned the object size
            pruned_ev = events[event_level]  # pruned events
            pruned_ev.Electron = (
                sele  # replace electron collections with selected electron
            )
            # Add custom variables
            if not isRealData:
                pruned_ev["weight"] = weights.weight()
                for ind_wei in weights.weightStatistics.keys():
                    pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                        include=[ind_wei]
                    )
            # Create a list of variables want to store. For objects from the PFNano file, specify as {object}_{variable}, wildcard option only accepted at the end of the string
            out_branch = np.setdiff1d(
                np.array(pruned_ev.fields), np.array(events.fields)
            )  # stored customed variables
            out_branch = np.append(out_branch, ["Jet_btagDeep*", "Electron_pt"])
            # write to root files
            os.system(f"mkdir -p {self.name}/{dataset}")
            with uproot.recreate(
                f"{self.name}/{dataset}/f{events.metadata['filename'].split('_')[-1].replace('.root','')}_{systematics[0]}_{int(events.metadata['entrystop']/self.chunksize)}.root"
            ) as fout:
                fout["Events"] = uproot_writeable(pruned_ev, include=out_branch)
        return {dataset: output}

    ## post process, return the accumulator, compressed
    def postprocess(self, accumulator):
        return accumulator
