import tkinter as tk
from tkinter import *

import matplotlib

matplotlib.use('TkAgg',force=True)
from matplotlib import pyplot as plt
print("Switched to:",matplotlib.get_backend())
import numpy as np
from pymed import PubMed

root = tk.Tk()
root.geometry("600x600")

#To select the datafile used
Files = sorted(["hg18","hg19","hg38"])
selected_file = StringVar()
selected_file.set("hg18")
dropdown = OptionMenu(root, selected_file, *Files)
dropdown.place(x=150,y=200)

#User selected data file type
def UserSelectedFile():
    FileChoice = selected_file.get()
    FileInAction = ""
    if FileChoice == "hg18":
        FileInAction = "NCBI36_hg18_variants_2020-02-25.txt"
    elif FileChoice == "hg19":
        FileInAction = "GRCh37_hg19_variants_2020-02-25.txt"
    else:
        FileInAction = "GRCh38_hg38_variants_2020-02-25.txt"
    return FileInAction




#Gene dict
#Why Hg38 is not working???
FileInAction=UserSelectedFile()
geneDict={}
with open(FileInAction,"r") as dataFile:
    for line in dataFile:
        if ("variantaccession" not in line) and "CNV" in line:
            line = line.split("\t")
            if "," not in line[-2]:
                geneDict[line[-2]]=[]

            else:
                genes=line[-2].split(",")
                for gene in genes:
                    geneDict[gene]=[]


with open(FileInAction,"r") as dataFile:
    for line in dataFile:
        if ("variantaccession" not in line) and "CNV" in line:
            line = line.split("\t")
            if "," not in line[-2]:
                geneDict[line[-2]].append(line[0])

            else:
                genes=line[-2].split(",")
                for gene in genes:
                    geneDict[gene].append(line[0])

#CNV sizes------------------------------------------------------------------------------------------------
def cnvSize():
    FileInAction=UserSelectedFile()

    seqLength = []
    with open(FileInAction, "r") as dataFile:
        for line in dataFile:
            if ("variantaccession" not in line) and "CNV" in line:
                line = line.split("\t")
                seqLength.append(int(line[3]) - int(line[2]))

    print(min(seqLength))
    print(max(seqLength))

    majoritySeqLength = []


    for length in seqLength:
        if int(text_from.get(1.0, "end-1c"))<length<int(text_to.get(1.0, "end-1c")):
            majoritySeqLength.append(length)


    plt.hist(majoritySeqLength, bins=int(text_bins.get(1.0, "end-1c")), edgecolor='black')
    plt.xlabel('Length')
    plt.ylabel('Frequency')
    plt.title('Histogram for Lengths of CNVs between '+str(int(text_from.get(1.0, "end-1c")))+" and "+str(int(text_to.get(1.0, "end-1c"))))
    plt.show()


#Possible CNVs for each gene
def cnvforGene():
    text_box.delete("1.0", "end")

    FileInAction=UserSelectedFile()

    # Gene dict
    geneDict = {}
    with open(FileInAction, "r") as dataFile:
        for line in dataFile:
            if ("variantaccession" not in line) and "CNV" in line:
                line = line.split("\t")
                if "," not in line[-2]:
                    geneDict[line[-2]] = []

                else:
                    genes = line[-2].split(",")
                    for gene in genes:
                        geneDict[gene] = []

    with open(FileInAction, "r") as dataFile:
        for line in dataFile:
            if ("variantaccession" not in line) and "CNV" in line:
                line = line.split("\t")
                if "," not in line[-2]:
                    geneDict[line[-2]].append(line[0])

                else:
                    genes = line[-2].split(",")
                    for gene in genes:
                        geneDict[gene].append(line[0])

    gene=selected_gene.get()
    for k,v in geneDict.items():
        if k==gene:
            text_box.insert(tk.END, "CNVs present in gene "+k + "\n")
            for cnv in v:
                text_box.insert(tk.END,cnv)
                text_box.insert(tk.END,"\n")


#Single function for chr/var/subvar data
def UserSelectedDataOption():
    userneeddata=selected_DataOption.get()
    colValue=0
    if userneeddata=="CNVs by chromosome":
        colValue=1
    elif userneeddata=="Variable Types of Data":
        colValue=4
    else:
        colValue=5

    text_box.delete("1.0", "end")

    FileInAction = UserSelectedFile()

    selected = selected_radio.get()
    userDict = {}
    key = ""
    val = 0
    with open(FileInAction, "r") as dataFile:
        for line in dataFile:
            if ("variantaccession" not in line):
                line = line.split("\t")
                key = line[colValue]
                userDict[key] = 0

    with open(FileInAction, "r") as dataFile:
        for line in dataFile:
            if ("variantaccession" not in line):
                line = line.split("\t")
                key = line[colValue]
            for k, v in userDict.items():
                if k == key:
                    v += 1
                    userDict[k] = v


    labels = list(userDict.keys())
    sizes = list(userDict.values())

    if selected == "Pie Chart":
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
        plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.title('Chromosome Data Availability')
        plt.show()

        for k, v in userDict.items():
            text_box.insert(tk.END, k + ":" + str(round((100 * v / sum(userDict.values())), 2)) + "%\n")

    else:
        fig = plt.figure(figsize=(20, 8))

        # creating the bar plot
        plt.bar(labels, sizes, color='maroon', width=0.2)

        plt.xlabel("Chr Number")
        plt.ylabel("Amount of data")
        plt.xticks(rotation=90)
        plt.title("Chromosome data availability in the database")
        plt.show()

        sumvar = 0
        for k, v in userDict.items():
            sumvar += v
            text_box.insert(tk.END, k + ":" + str(v) + "\n")
        text_box.insert(tk.END, "Sum:" + str(sumvar))




#-----------------------------------------------------------------------------------------------------------
#
# def UserSelectedFile():
#     FileChoice = selected_file.get()
#     FileInAction = ""
#     if FileChoice == "hg18":
#         FileInAction = "NCBI36_hg18_variants_2020-02-25.txt"
#     elif FileChoice == "hg19":
#         FileInAction = "GRCh37_hg19_variants_2020-02-25.txt"
#     else:
#         FileInAction = "GRCh38_hg38_variants_2020-02-25.txt"
#     return FileInAction




#----------------------------------------------------------------------------------------------------------

label0=tk.Label(root,text="Select Data").place(x=10,y=200)
label1 = tk.Label(root,text="Select graph you want \nto retrieve").place(x=10,y=250)
# label2 = tk.Label(root,text="Variable Types of Data").place(x=10,y=224)
# label3 = tk.Label(root,text="Subvariable Types of Data").place(x=10,y=276)
label4 = tk.Label(root,text="Number of CNVs\n between the lengths").place(x=10,y=328)
label5 = tk.Label(root,text="From").place(x=150,y=328)
label6 = tk.Label(root,text="To").place(x=235,y=328)
label7 = tk.Label(root,text="Bins").place(x=310,y=328)
label8 = tk.Label(root,text="CNVs by gene").place(x=10,y=400)



text_box = tk.Text(root, bg='white', bd=3, cursor='hand2', height=10.4,width=65)
text_box.place(x=10,y=10)

scrollbar = tk.Scrollbar(root, orient="vertical", command=text_box.yview)
scrollbar.grid(row=0, column=1, sticky="ns")

text_box.config(yscrollcommand=scrollbar.set)


#Radio Buttons for chromosomes
selected_radio = tk.StringVar()
radio_bar = tk.Radiobutton(root, text='Count', value='Bar Graph', variable=selected_radio,command=UserSelectedDataOption).place(x=350,y=250)
radio_pie = tk.Radiobutton(root, text='Percentage', value='Pie Chart', variable=selected_radio,command=UserSelectedDataOption).place(x=425,y=250)

# #Radio Buttons for variations
# selected_var=tk.StringVar()
# radio_var_bar = tk.Radiobutton(root, text='Count', value='Bar Graph', variable=selected_var,command=varGraphs).place(x=150,y=224)
# radio_var_pie = tk.Radiobutton(root, text='Percentage', value='Pie Chart', variable=selected_var,command=varGraphs).place(x=250,y=224)
#
# #Radio Buttons for sub-variations
# selected_subvar=tk.StringVar()
# radio_subvar_bar = tk.Radiobutton(root, text='Count', value='Bar Graph', variable=selected_subvar,command=subvarGraphs).place(x=150,y=276)
# radio_subvar_pie = tk.Radiobutton(root, text='Percentage', value='Pie Chart', variable=selected_subvar,command=subvarGraphs).place(x=250,y=276)

#Widgets for CNV lengths
text_from=tk.Text(root, bg='white', bd=3, cursor='hand2', height=1,width=5)
text_from.place(x=185,y=328)
text_to=tk.Text(root, bg='white', bd=3, cursor='hand2', height=1,width=5)
text_to.place(x=260,y=328)
text_bins=tk.Text(root, bg='white', bd=3, cursor='hand2', height=1,width=5)
text_bins.place(x=335,y=328)
length_button=tk.Button(root,text="Create\nHistogram",command=cnvSize).place(x=400,y=328)


#Widgets for cnvs for genes
options = sorted(list(geneDict.keys()))  # or list(options_dict.values())

selected_gene = StringVar()

dropdown = OptionMenu(root, selected_gene, *options)
dropdown.place(x=150,y=400)

drop_button=tk.Button(root,text="Select gene",command=cnvforGene).place(x=270,y=400)

#Widgets for chromosomes, variant types and sub variant types
DataOptions = sorted(["CNVs by chromosome","Variable Types of Data","Sub-variable Types of Data"])  # or list(options_dict.values())

selected_DataOption = StringVar()
selected_DataOption.set("CNVs by chromosome")
dropdown = OptionMenu(root, selected_DataOption, *DataOptions)
dropdown.place(x=150,y=250)















root.mainloop()





