import numpy as np
from proteingraph import read_pdb
from proteingraph.pin import filter_dataframe,compute_chain_pos_aa_mapping,compute_rgroup_dataframe,compute_distmat,get_interacting_atoms
import torch
from os.path import exists
from torch_geometric.data import Data
from proteingraph.pin import pdb2df
from proteingraph import read_pdb
## proteingraph version "0.3.0"



def number_of_nodes(pdb_path,pdb_id):
    G = read_pdb(str(pdb_path)+'/'+str(pdb_id)+'.pdb')
    return G.number_of_nodes()

def number_of_edges(pdb_path,pdb_id):
    G = read_pdb(str(pdb_path)+'/'+str(pdb_id)+'.pdb')
    return G.number_of_edges()

def edge_result(G,res2node):
    start_aa=[]
    target_aa=[]
    for edge in G.edges():
        residue_start=str(edge[0][:-3])
        residue_target=str(edge[1][:-3])
        node_start=res2node[residue_start]
        node_target=res2node[residue_target]
        start_aa.append(int(node_start))
        target_aa.append(int(node_target))
    result_edges=[start_aa,target_aa]
    return result_edges

def residueID2nodeID(pdb_df):
    chain_pos_aa=compute_chain_pos_aa_mapping(pdb_df)
    residue2node={}
    node_index=0

    for i in chain_pos_aa.keys():
        my_chain=chain_pos_aa[i]
        for item in my_chain:
            residue2node[str(i)+str(item)]=str(node_index)
            node_index+=1
    return residue2node

def node_result(G,res2node):
    sequence=['N']*G.number_of_nodes()
    for node in G.nodes():

        residueID=str(node[:-3])
        aa_name=str(node[-3:])
        nodeID=res2node[residueID]
        sequence[int(nodeID)]=aa_name
    return sequence


def edge_feature(G):
    edge_type=[]
    for _,_,d in G.edges(data=True):
        edge_type.append(d['kind'])
    return edge_type

def add_pi_pi_interactions(G,pdb_df):
    """
    Find all (non-aromatic) pi-pi interactions.
    Performs searches between the following residues:
    "TYR","PHE","TRP","HIS","GLN","ASN","GLU","ASP","ARG","GLY","SER","THR","PRO"
    Criteria: R-group residues are within 6A distance.
    """
    non_aromatic_pi_RESIS=["TYR","PHE","TRP","HIS","GLN","ASN","GLU","ASP","ARG","GLY","SER","THR","PRO"]
    rgroup_df = compute_rgroup_dataframe(pdb_df)
    non_aromatic_pi_df = filter_dataframe(rgroup_df, "residue_name", non_aromatic_pi_RESIS, True)
    distmat = compute_distmat(non_aromatic_pi_df)
    interacting_atoms = get_interacting_atoms(6, distmat)

    # this function is NOT heritaged from proteingraph.pin, its definition is below
    add_interacting_resis_additional(G, interacting_atoms, non_aromatic_pi_df, ["pipi"])

def add_interacting_resis_additional(G, interacting_atoms, dataframe, kind): 
    # modify this function from pin.py in order to add new bond types
    resi1 = dataframe.loc[interacting_atoms[0]]["node_id"].values
    resi2 = dataframe.loc[interacting_atoms[1]]["node_id"].values

    interacting_resis = set(list(zip(resi1, resi2)))
    for i1, i2 in interacting_resis:
        if i1 != i2:
            if G.has_edge(i1, i2):
                for k in kind:
                    if k not in G.edges[i1, i2]["kind"]:
                        G.edges[i1, i2]["kind"].append(k)
            else:
                G.add_edge(i1, i2, kind=list(kind))



def convert_pdb2graph(input):
 
    pdb_path=input[0]
    my_protein=input[1]
    featureData=input[2]
    graph_labels=input[3]
    protein_index=input[4]
    type_input=input[5]
    if type_input=='AlphaFold':
        path=str(pdb_path)+'/'+'AF-'+str(my_protein)+'-F1-model_v3.pdb'
    elif type_input=='PDB':
        path=str(pdb_path)+'/'+str(my_protein)+'.pdb'
    else:
        print('Unkown input type, must be AlphaFold or PDB')
        return 
    if exists(path):
        try:
            pdb_df=pdb2df(path)
            res2node=residueID2nodeID(pdb_df)
            if type_input=='AlphaFold':
                G=read_pdb(str(pdb_path)+'/'+'AF-'+str(my_protein)+'-F1-model_v3.pdb')
            elif type_input=='PDB':
                G=read_pdb(str(pdb_path)+'/'+str(my_protein)+'.pdb')
            else:
                print('Unkown input type, must be AlphaFold or PDB')
                return
            add_pi_pi_interactions(G,pdb_df)

            node=node_result(G,res2node)
            edge=edge_result(G,res2node)
            print("Loaded ",str(my_protein))
        except (IndexError, ValueError):
            print("Can't load ",str(my_protein))
            return

        
    ### readin feature list of all amino acids
        try:
            complete_list_feature=[]
            for aminoAcid in featureData.buildProtein(node):
                my_features=[float(tmp) for tmp in list(aminoAcid)]
                complete_list_feature.append(my_features)

            nodes_features=torch.tensor(complete_list_feature,dtype=torch.float) # feature vector of all nodes
            edge_index = torch.tensor(edge, dtype=torch.long) # edges, 1st list: index of the source nodes, 2nd list: index of target nodes.
            my_label=torch.tensor(graph_labels[protein_index], dtype=torch.long)
            g = Data(x=nodes_features, edge_index=edge_index,y=my_label)

            return g
        except KeyError:
            print("Failed loading aminoacids info for ",str(my_protein))
            return

