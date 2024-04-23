

# for find out the number of bases
# returns the four integer type output
def no_of_code(nucleotide_seq):
    # gets the parameter as DNA sequence not "aplied for RNA" 
    (no_A, no_T, no_G, no_C) = (0, 0, 0, 0)
    for i in nucleotide_seq:
        if i == 'A':
            no_A += 1
        if i == 'T':
            no_T += 1
        if i == 'G':
            no_G += 1
        if i == 'C':
            no_C += 1
    return no_A, no_T, no_G, no_C  


# returns aproxmatly weight of protein molecules 
# returns integer type output value   
def peptide_Molecular_weight(protein_seq):
    no_of_amino_acid = len(protein_seq)
        # multiply with mean value of amino acid in kiloDalton 
    mw = no_of_amino_acid * 0.11
    return mw     # answer is given in "kilo Dalton"


# for finding the value of moleculer weight of single standerd of DNA
# returns the integer type of single value output   
def ssDNA_MoleculerWeight(nucleotide_seq):   # parameter for DNA
    (no_A, no_T, no_G, no_C) = no_of_code(nucleotide_seq)
    MW_ssDNA = (no_A * 313.2) + (no_T * 304.2) + (no_G * 329.2) + (no_C * 289.2) + 79.0
    return MW_ssDNA     # returns moleculer weight in "gram/mol"


# short function for finding temperature on which melt the perticular DNA sequence 
# returns single output of integer type
def Melting_Temperature(DNA_sequence):
    # input are checked by as a DNA sequence "not for RNA"
    a, t, g, c = no_of_code(DNA_sequence)
    # eqution  returns °C answer
    return 2*a+t + 4*g+c   # returns answer in "°C"

def ssRNA_MoleculerWeight(nucleotide_seq):
    # takes single input as mRNA sequence parameter
    (no_A, no_U, no_G, no_C) = no_of_code(nucleotide_seq)
    MW_ssRNA = (no_A * 329.2) + (no_U * 306.2) + (no_G * 345.2) + (no_C * 305.2) + 159.0
    return MW_ssRNA     # answer will apear in "gram/mole"


# problem is there conversion in temperature
# using simple formula of temp. conversion
class Temperature_Conversion():
    def celsius_to_kelvin(self, celsius_):
        # return answer in K
        return celsius_ + 273.15
    
    def kelvin_to_celsius(self, kelvin_):
        # return answer in °C
        return kelvin_ - 273.15

    def fahrenheit_to_celsius(self, fahrenheit_):
        # # return answer in °C
        return fahrenheit_ - 32 * 0.55
    
    def celsius_to_fahrenheit(self, celsius_):
        # return answer in °F
        return celsius_ * 0.55 + 32

    def fahrenheit_to_kelvin(self, fahrenheit_):
        # # return answer in K
        return fahrenheit_ - 32 * 0.55 + 273.15
    
    def kelvin_to_fahrenheit(self, kelvin_):
        # return answer in °F
        return kelvin_ - 273.15 * 0.55 + 32


# for read emprical formula of a chemical structure, only compatible on organic compound
# to get input as a emprical formula of drug and gives output as a two list:-
# first is number of element which is contain in formula, and second one is name of the element  
def read_empirical(empirical_formula):
    import re
    # re is the regular expression for finding numeric in emprical
    # empirical like this = 'C10H16N2O2'
    element = []
    # finding the organic element into emprical 
    for i in empirical_formula:
        if i == 'C':            # find carbon
            element.append(i)
        elif i == 'H':          # find hydorgen
            element.append(i)
        elif i == 'N':          # find nitrogen
            element.append(i)
        elif i == 'O':          # find oxygen
            element.append(i)
        elif i == 'F':          # find florine
            element.append(i)
        elif i == 'I':          # find iodine
            element.append(i)
    # sepratly find bromine, clorine  
    if 'Br' in empirical_formula:
        element.append('Br')
    if 'Cl' in empirical_formula:
        element.append('Cl')
    # find how many element present in a particular compound 
    no_contant_list = re.findall(r"\d+", empirical_formula)
    # return number of element and name of element
    return no_contant_list, element


# block of code which find molecular weight using emprical formula('read_empirical' function above)
# it's take input as return function of read_empirical and you gets output a integer type mol. weight 
def Organic_MolecularWeight(no_contant_list, element):
    no_cont = 0
    no_element = 0
    # count number of containt and element 
    # for genrating a condition, is both are equal -> so find the mol. weight at the line 134.
    for a in no_contant_list:
        no_cont += 1
    for b in element:
        no_element += 1
    if no_cont == no_element:
        periodic = {'H': 1.007, 'C': 12.01, 'N': 14.006, 'O': 15.999, 'S': 32.065,
                    'F': 18.998, 'Cl': 35.453, 'Br': 79.904, 'I': 126.90}
        # there is the mol. weight of a perticular organic element
        Mw = 0
        # initially mol. weight is zero
        for ele, c in zip(element, no_contant_list):
            if ele in periodic:
                # multiply element, if it's present in periodic dictionary 
                # and perticular number of element
                Mw += periodic[ele]*int(c)
        # return mol. weight
        return Mw
    # if both are not equal -> so add 1 on number of element list type
    # means only single element is present, not multiple form
    else:
        diff = no_element - no_cont
        for n in range(0, diff):
            no_contant_list.append(1)
        # at the last recall all over function(by recursion)
        return Organic_MolecularWeight(no_contant_list, element)


def Chem_MolecularWeight(emprical_formula):
    no_contant_list, element = read_empirical(emprical_formula)
    return Organic_MolecularWeight(no_contant_list, element)

