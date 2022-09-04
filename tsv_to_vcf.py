import numpy as np
import os
import shutil
import pandas as pd
from tqdm import tqdm

def filter_tsv(filename):
    tsv_info = pd.read_csv(filename, sep='\t', usecols= [0,2,3], header=0)
    # Extract relevant information from the tsv file
    tsv_info['Chromosome'] = tsv_info['DNA Change'].str.split(':').str[0]
    tsv_info['Pos'] = tsv_info['DNA Change'].str.split(':').str[1].str.split('>').str[0]
    tsv_info['Gene ID'] = tsv_info['Consequences'].str.split().str[1]
    tsv_info['Cons_extra'] = tsv_info['Consequences'].str.split().str[2]
    tsv_info['Pos'] = tsv_info['Pos'].str.replace('g.', '')
    tsv_info['Position'] = tsv_info['Pos'].str[:-1]
    tsv_info['Reference'] = tsv_info['Pos'].str[-1]
    tsv_info['Alteration'] = tsv_info['DNA Change'].str.split('>').str[1]
    tsv_info['fraction'] = tsv_info['# Affected Cases in Cohort'].str.split(',').str[1].str[:-1].astype(float)/100
    # create a new dataframe with specific columns from tsv_info
    relevant = tsv_info[['Chromosome', 'Position', 'Gene ID', 'Reference', 'Alteration', 'fraction']].copy()
    relevant['Position'] = relevant['Position'].astype(int)
    relevant['fraction'] = relevant['fraction'].astype(float)
    relevant['Chromosome'] = relevant['Chromosome'].str.split('r').str[1]
    relevant['Chromosome'] = relevant['Chromosome'].str.replace('X', '30')
    #relevant['Chromosome'] = relevant['Chromosome'].str.replace('Y', '31')
    relevant['Chromosome'] = relevant['Chromosome'].astype(int)
    relevant.sort_values(by=['Chromosome','Position'],inplace=True)
    relevant['Chromosome'] = relevant['Chromosome'].astype(str)
    relevant['Chromosome'] = relevant['Chromosome'].str.replace('30', 'X')
    #relevant['Chromosome'] = relevant['Chromosome'].str.replace('31', 'Y')
    relevant = relevant.dropna(axis=0)
    relevant.reset_index(drop=True, inplace=True)
    # save relevant information to a csv file
    relevant.to_csv("relevant_info.csv",index=False)
    return relevant

def generate_variants(relevant,num_samples,category):
    new_seq = np.empty((relevant.shape[0],num_samples)).astype(str)
    for i in tqdm(range(relevant.shape[0])):
        rng = np.random.default_rng()
        newp = np.array(["0|0"]*num_samples)
        if relevant['Gene ID'][i] == 'ERBB2':
            newp[:] = "1|1"
        elif relevant['Gene ID'][i] == 'PIK3CA':
            if category == 'all':
                newp[:] = "1|1"
            elif category == 'before':
                mut = np.arange(0,num_samples,5)
                newp[mut] = "1|1"
            elif category == 'none':
                pass
                #newp[[j for j in range(num_samples) if j%2==0]] = "1|1"
        else:
            prob = (relevant.iloc[i,-1]*num_samples).astype(np.int32)
            mut = rng.choice(num_samples,prob,replace=False)
            assert mut.shape[0] == len(set(mut))        
            newp[mut] = "1|1"
        new_seq[i,:] = newp
    simulated_seq = pd.DataFrame(new_seq,columns=range(1,num_samples+1))
    return simulated_seq

def create_vcf(subject_index,relevant,sim,startText,num_samples,category):    
    vcf_data = ""
    # add header to vcf_data
    ll = len(str(num_samples))
    if category == 'none':
        subject_ID = str(subject_index+1).zfill(ll)
    else:
        subject_ID = str(subject_index+(num_samples//2 + 1)).zfill(ll)

    if category == 'before':
        if subject_index%5 == 0:
            vcf_data = startText + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+subject_ID+"_10\n"
        else:                
            vcf_data = startText + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+subject_ID+"_40\n"
    else:
        vcf_data = startText + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+subject_ID+"\n"

    relevant['Combined'] = relevant['Chromosome'] + "\t" + relevant['Position'].astype(str) + \
    "\t" + relevant['Gene ID'] + "\t" + relevant['Reference'] + "\t" + relevant['Alteration'] + \
    "\t.\t.\t.\tGT\t" + sim[subject_index+1] #sim[str(subject_index+1)]

    vcf_data += "\n".join(relevant['Combined'])
    # save vcf file
    if category == 'none':
        with open("all_files/no_pik3ca/subject_"+subject_ID+".vcf", "w") as text_file:
            text_file.write(vcf_data)
    if category == 'all':
        with open("all_files/both_genes/subject_"+subject_ID+".vcf", "w") as text_file:
            text_file.write(vcf_data)
    elif category == 'before':
        if subject_index%5 == 0:
            with open("all_files/ten_percent/subject_"+subject_ID+"_10.vcf", "w") as text_file:
                text_file.write(vcf_data)
        else:
            with open("all_files/forty_percent/subject_"+subject_ID+"_40.vcf", "w") as text_file:
                text_file.write(vcf_data)

def tsv_2_vcf(filename,startfile,num_samples):
    # read relevant from csv file. If it doesn't exist, create it.
    if os.path.isfile("relevant_info.csv"):
        relevant = pd.read_csv("relevant_info.csv")
    else:
        relevant = filter_tsv(filename)
    print("Filtering complete")
    r = relevant.copy()

    # read header from startfile and save it to a variable
    with open(startfile, 'r') as f:
        startText = f.read()

    os.mkdir("all_files") # create a directory to store all vcf files
    os.mkdir("all_files/no_pik3ca") # throws error if folder already exists
    os.mkdir("all_files/both_genes") # throws error if folder already exists
    os.mkdir("all_files/ten_percent") # throws error if folder already exists
    os.mkdir("all_files/forty_percent") # throws error if folder already exists

    sim_none = generate_variants(relevant,num_samples//2,'none')
    print("Creating files with no PIK3CA mutations")
    for i in tqdm(range(num_samples//2)):
        create_vcf(i,r,sim_none,startText,num_samples,'none')

    sim_all = generate_variants(relevant,num_samples//2,'all')
    print("Creating files with all PIK3CA mutations")
    for i in tqdm(range(num_samples//2)):
            create_vcf(i,r,sim_all,startText,num_samples,'all')

    sim_before = generate_variants(relevant,num_samples//2,'before')
    print("Creating files where some have PIK3CA mutations and some don't")
    for i in tqdm(range(num_samples//2)):
            create_vcf(i,r,sim_before,startText,num_samples,'before')
    
def reorganize(num_samples,max_files,main_folder):
    os.mkdir(main_folder) # throws error if folder already exists
    a = int(max_files//14)*14
    b = int(a*0.4/1.4)
    c = (a-b)//2
    d = round(c*2.8)
    assert d == a
    if num_samples <= c*2:
        num_folders = 1
    else:
        num_folders = num_samples//(c*2)
    # Create num_folders number of sub-folders in data_for_babyships folder
    for i in range(num_folders):
        folder_name = "set"+str(i+1)+"_simulated_gastric_guardant_panel"
        os.mkdir(main_folder+"/"+folder_name)

    # Moving all files to the correct folder
    for i in range(num_folders):
        folder_name = "set"+str(i+1)+"_simulated_gastric_guardant_panel"

        # Copy every c files from no_pik3ca folder to the new folders
        for file in os.listdir("all_files/no_pik3ca"):
            if file != ".DS_Store":
                if i*c < int(file.split('.')[0].split('_')[1]) <= (i+1)*c:
                    shutil.copy("all_files/no_pik3ca/"+file,main_folder+"/"+folder_name+"/"+file)
        
        # Copy every c files from both_genes folder to the new folders
        for file in os.listdir("all_files/both_genes"):
            if file != ".DS_Store":
                if i*c < int(file.split('.')[0].split('_')[1])-num_samples//2 <= (i+1)*c:
                    shutil.copy("all_files/both_genes/"+file,main_folder+"/"+folder_name+"/"+file)

        # Repeat the above steps for the 40 percent folder    
        for file in os.listdir("all_files/forty_percent"):
            if file != ".DS_Store":
                if i*c < int(file.split('_')[1])-num_samples//2 <= (i+1)*c:
                    shutil.copy("all_files/forty_percent/"+file,main_folder+"/"+folder_name+"/"+file)

def main():
    #Create vcf files from the given tsv file
    tsv_2_vcf(filename='frequent-mutations_GDC_GastricCancer_SNPs.tsv', startfile='startText.txt', num_samples=100000)
    #Reorganize the files into the correct format
    reorganize(num_samples=100000, max_files=4000, main_folder = "data_for_babyships")

if __name__ == "__main__":
    main()