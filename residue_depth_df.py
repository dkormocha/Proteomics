from Bio.PDB.ResidueDepth import ResidueDepth
from Bio.PDB.PDBParser import PDBParser
import pandas as pd 
import os 


def residue_depth(file): #This function calculates the residue depth of a given PDB file and makes a dataframe with the results
    parser = PDBParser() #Creates a parser object
    structure = parser.get_structure(f'{file}', file)   #Creates a structure object from the PDB file
    model = structure[0] #Creates a model object from the structure object
    rd = ResidueDepth(model)    #Creates a ResidueDepth object from the model object

    appender = []   #Creates a list to append the results to
    for k in rd.property_keys:  #Iterates through the keys in the ResidueDepth object
        x = rd.property_dict[k] #Creates a variable to store the value of the key
        resdepth = x[0] #Creates a variable to store the value of the first element of the list
        cadepth = x[1]  
        appender.append((resdepth, cadepth))    #Appends the results to the list

    df_rd = pd.DataFrame.from_records(appender, columns=['res_depth', 'ca_depth'])  #Creates a dataframe from the list of results
    df_rd['index'] = range(1, len(df_rd) + 1)   #Adds an index column to the dataframe
    return df_rd


def get_protein_name(df): #Function uses proteins in data frame to find all the PDB files in the directory and run the residue depth function on them
    proteinID= [] #List to store the protein IDs
    for p in df_dssp['ProteinID']:  #Loop through the proteins in the data frame
        proteinID.append(p) #Append the protein ID to the list

    proteins= set(proteinID)   #Convert the list to a set to remove duplicates
    for k in proteins:  #Loop through the set of proteins
        print(k)
        filename= f'AF-{k}-F1-model_v2.pdb' #Create the filename
        if filename in os.listdir('human_proteome/'):   #Check if the file is in the directory
            print(filename)
            df_rd= residue_depth(f'human_proteome/{filename}')  #Run the residue depth function on the file
            df_rd.to_csv(f'rd/{filename}.csv')  #Save the dataframe to a csv file




def add_values(df): #Function to add the values of the residue depth and cadepth to the dataframe

    res_depth = [] #List to store the residue depth values
    ca_depth = []   
    for i,j in df_dssp.iterrows():  #Loop through the dataframe
        file= f'AF-{j["ProteinID"]}-F1-model_v2.pdb.csv'    #Create the filename

        if file in os.listdir('rd/'):   #Check if the file is in the directory

            df_rd= pd.read_csv(f'rd/{file}')    #Read the csv file
            if j['LastAA_position'] in df_rd['index']:  #Check if the last amino acid position is in the dataframe
                res_depth.append(df_rd.loc[df_rd['index'] == j['LastAA_position']]['res_depth'].values[0])  #Append the residue depth value to the list
                ca_depth.append(df_rd.loc[df_rd['index'] == j['LastAA_position']]['ca_depth'].values[0])    #Append the cadepth value to the list
            else:
                res_depth.append('') #Append a blank value to the list
                ca_depth.append('')  
        else:
            res_depth.append('')    #Append a blank value to the list
            ca_depth.append('') 

    print(res_depth)

    df_dssp['res_depth'] = res_depth    #Add the values to the dataframe
    df_dssp['ca_depth'] = ca_depth   #Add the values to the dataframe



df_dssp = pd.read_csv('df_dssp.csv')

add_values(df_dssp)

df_dssp.to_csv('df_dssp1.csv')