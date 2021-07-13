import pandas as pd
from pyteomics import mass
import copy

SEQUENCE = 'KQLATKAAR'

DICT_MODIFICATIONS = {
    'ac12' : mass.calculate_mass(formula='C[12]2H3O') - mass.calculate_mass(formula='H'),
    'ac13' : mass.calculate_mass(formula='C[13]2H3O') - mass.calculate_mass(formula='H'),
    'ac*' : mass.calculate_mass(formula='C[13]2H[2]3O') - mass.calculate_mass(formula='H')
}

ROUND_POSITION = 11

def nextCombination(listOfElements, combination, length):
    """Gets the next combination of a given length and elements
    
    Parameters
    ----------
    listOfElements : list
        A list of elements allowed in the combination
    combination : list
        A list of elements contained in lisOfElements of length length
    length : int
        The length of the combination
        Necessary to perform recursion

    Returns
    -------
    list
        a list of the next combination
    """

    indexChecked = length - 1 
    if combination[indexChecked] != listOfElements[-1]:        
        indexElement = listOfElements.index(combination[indexChecked])
        combination[indexChecked] = listOfElements[indexElement + 1]      
        return combination

    combination[indexChecked] = listOfElements[0]    
    return nextCombination(listOfElements, combination, length-1)

def combinations(listOfElements, length):
    """Generates all combinations of given elements with given length 

    Parameters
    ----------
    listOfElements : list
        A list of elements allowed in the combinations
    length : int
        The length of each combination

    Returns
    -------
    list
        Nested list. Returns a list of each possible combination, which in turn are also a list
    """
    listOfElements = list(listOfElements)
    size = len(listOfElements) ** length
    #generate first combination
    start = [[]]
    for i in range(length):
        start[0].append(listOfElements[0])
    
    #generate the next combination of the last element and append it to the return list
    for j in range(size-1):
        next = nextCombination(listOfElements, copy.deepcopy(start)[-1], length)
        start.append(next)

    return start

def precursorMZ(sequence, modifications, charge):
    """Calculates the m/z of a given peptide with a list of modifications

    Parameters
    ----------
    sequence : str
        The sequence of a peptide as one-letter code and in uppercase
    modifications : list
        a list of modifications as strings. Modifications must be included in DICT_MODIFICATIONS
    charge : int
        charge of the peptide (to be automated in future versions to calculate it automatically)

    Returns
    -------
    float
        m/z of peptide
    """

    precursorMass = mass.calculate_mass(sequence = sequence) + charge * mass.calculate_mass(formula = 'H+')
    for mod in modifications:
        precursorMass += DICT_MODIFICATIONS[mod]

    return precursorMass / charge

def modifiedSequence(sequence, modifications, aminoacid, ion_type):
    """modifies string of the sequence

    Modifies the string of a given sequence that it contains the modifications
    of a given aminoacid in brackets

    Parameters
    ----------
    sequence : str
        sequence of the peptide or fragment (one-letter code, uppercase)
    modifications : list
        list of strings of the modifications in the correct order
    aminoacid : str
        modified aminoacid as one-letter code
    ion-type : str
        ion type of the fragment (currently b/y, to be added: M, a/z, c/x)
        
    Returns
    -------
    str
        string of the modified sequence        
    """

    split = sequence.split(aminoacid)
    modified = ''    
    mod = copy.deepcopy(modifications)
    if ion_type == 'y':
        if sequence.count(aminoacid) < len(mod):
            difference = len(mod) - sequence.count(aminoacid)            
            for i in range(difference):
                mod.pop(0)   

        for i in range(len(split)):
            if i == len(split) - 1:
                modified = modified + split[i]
            else:
                modified = modified + split[i]+ aminoacid + '(' + mod[i] + ')'
    else:
        for i in range(len(split)):
            if i == len(split) - 1:
                modified = modified + split[i]
            else:
                modified = modified + split[i]+ aminoacid + '(' + mod[i] + ')'
    return modified
    
def ionMass(sequence, modifications, aminoacid, ion_type):
    """Calculates mass of a modified fragment

    Currently just charge +1. Automatic charge calculation to be added.

    Parameters
    ----------
    sequence : str
        sequence of fragment (one-letter code, uppercase)
    modifications : list
        list of strings of the modifications in the correct order
    aminoacid : str
        modified aminoacid as one-letter code
    ion-type : str
        ion type of the fragment (currently b/y, to be added: M, a/z, c/x)

    Returns
    -------
    float
        m/z of modified fragment
    """

    count = sequence.count(aminoacid)
    fragmentMass = mass.calculate_mass(sequence = sequence, ion_type = ion_type, charge = 1)
    if ion_type == 'y':        
        for i in range(count):
            fragmentMass += DICT_MODIFICATIONS[modifications[-(i+1)]]
    else:
        for i in range(count):
            fragmentMass += DICT_MODIFICATIONS[modifications[i]]

    return fragmentMass

def fragmentIons(sequence, modifications, aminoacid):
    """Generate a pandas dataframe of all possible fragment ions, given a list of modifications of an aminoacid

    Parameters
    ----------
    sequence : str
        sequence of the peptide (one-letter code, uppercase)
    modifications : list
        list of strings of the allowed modifications
    aminoacid : str
        modified aminoacid (one-letter code, uppercase)

    Returns
    -------
    object
        returns a pandas dataframe with all possible fragment ions
        columns: peptide, precursor sequence, precursor m/z, fragment ion, fragment sequence, fragment m/z
    """

    fragmentIons = []
    modCombinations = combinations(modifications, sequence.count(aminoacid))    
    for modification in modCombinations:            
        for i in range(len(sequence)):
            ion = []
            fragmentSequence = sequence[i+1:]
            if fragmentSequence == '':
                break
            ion.append(sequence) #peptide/precursor sequence
            ion.append(modifiedSequence(sequence, modification, aminoacid, 'y')) #precursor sequence modified
            ion.append(precursorMZ(sequence, modification, 2)) #precursor m/z
            ion.append('y{number}'.format(number = len(sequence)-i-1)) #y-ion                    
            ion.append(modifiedSequence(fragmentSequence, modification, aminoacid, 'y')) #y ion sequence
            ion.append(ionMass(fragmentSequence, modification, aminoacid,ion_type = 'y')) #fragment m/z            
            ion.append(modification)
            fragmentIons.append(ion)
            
           
        for i in range(len(sequence)-1):
            ion = []
            ion.append(sequence) #peptide/precursor sequence
            ion.append(modifiedSequence(sequence, modification, aminoacid, 'b')) #precursor sequence modified
            ion.append(precursorMZ(sequence, modification, 2)) #precursor m/z
            fragmentSequence = sequence[:i+1]#b ions
            ion.append('b{number}'.format(number = i+1)) #b-ion                    
            ion.append(modifiedSequence(fragmentSequence, modification, aminoacid, 'b')) #b ion sequence            
            ion.append(ionMass(fragmentSequence, modification, aminoacid,'b')) #fragment m/z
            ion.append(modification)
            fragmentIons.append(ion)
            
          
    return pd.DataFrame(fragmentIons, columns=['peptide', 'precursor sequence', 'precursor m/z', 'fragment_ion', 'fragment sequence', 'fragment m/z', 'modification'])
   
def uniqueFragment(dataframe, row):
    """Auxiliary function for precurserSpecificIons()

    This function is intended to be applied on a pandas dataframe.

    Parameters
    ----------
    dataframe : object
        pandas dataframe of fragment ions
    row : object
        row of a pandas dataframe

    Returns
    -------
    bool
        Returns weather a fragment of a given row is unique or not        
    """
      
    return len(dataframe.loc[dataframe['fragment m/z'].round(ROUND_POSITION) == round(row['fragment m/z'], ROUND_POSITION)]['precursor sequence'].value_counts()) == 1

def precursorSpecificIons(dataframe):
    """Detects all fragment ions that are unique for one sequence for each precurser m/z

    Parameters
    ----------
    dataframe : object
        pandas dataframe containing at least the columns:
        precurser m/z, precurser sequence, fragment m/z
    
    Returns
    -------
    object
        pandas dataframe of all fragment ions thar are unique for one sequence for each precurser m/z
    """

    specific = pd.DataFrame()
    for mass in dataframe['precursor m/z'].unique():
        #generate dataframe of desired m/z
        df = dataframe.loc[dataframe['precursor m/z'].round(ROUND_POSITION) == round(mass,ROUND_POSITION)]

        df = df.loc[df.apply(lambda row: uniqueFragment(df, row), axis=1)]        
        specific = pd.concat([specific, df])        
    return specific


if __name__ == "__main__":
    df = fragmentIons(SEQUENCE, DICT_MODIFICATIONS, 'K')
    df.to_csv('{sequence}_Fragmentation.csv'.format(sequence = SEQUENCE), index=False)