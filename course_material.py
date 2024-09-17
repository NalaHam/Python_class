#review of python course material
----------------------------------------------------data structure-----------------------------------------------------------------------------
#https://colab.research.google.com/drive/1kMwDW2c2oG7s6X0cjQqvBvQI8wmO7-e4?usp=sharing#scrollTo=VpxAdelKZyaJ
#Create a dictionary with the following data: {"Homo sapiens": 10, "Pan troglodytes": 5, "Mus musculus": 15}.
#Write code to add a new species to the dictionary ("Danio rerio": 8).
#Print the total number of individuals observed across all species.

# Solution
species_counts = {"Homo sapiens": 10, "Pan troglodytes": 5, "Mus musculus": 15}

# Add new species
species_counts["Danio rerio"] = 8

# Calculate total number of individuals
total_individuals = sum(species_counts.values())

print("Total individuals observed:", total_individuals)

#Create a list of DNA sequences with some duplicates: ["ATCG", "GCTA", "ATCG", "CGTA", "GCTA"].
#Convert the list to a set to remove duplicates.
#Print the set of unique sequences.

# Solution
dna_sequences = ["ATCG", "GCTA", "ATCG", "CGTA", "GCTA"]

# Convert to set to remove duplicates
unique_sequences = set(dna_sequences)

print("Unique DNA sequences:", unique_sequences)

#Create two sets of species observed by each research team:
#team_1 = {"Homo sapiens", "Pan troglodytes", "Mus musculus"}
#team_2 = {"Pan troglodytes", "Drosophila melanogaster", "Mus musculus"}
#Find and print the species common to both sets.

# Solution
team_1 = {"Homo sapiens", "Pan troglodytes", "Mus musculus"}
team_2 = {"Pan troglodytes", "Drosophila melanogaster", "Mus musculus"}

# Find common species
common_species = team_1.intersection(team_2)

print("Species common to both studies:", common_species)

-------------------------------------------AI_troubleshooting---------------------------------------------------------------
#https://colab.research.google.com/drive/1u4SLSxppWV5rXeC84Zz94_cP2ToiQIo9?usp=sharing

-------------------------------------------If Statements---------------------------------------------------------------
#https://colab.research.google.com/drive/1LUmI02xUpd12OyrWwJ2fr_I7hCbD-7sU?usp=sharing#scrollTo=H61pTrmEivuE
#We may also want to test if the enzyme is inactive below a certain threshold, and handle both cases.
# if-else example
temperature = 35  # in degrees Celsius

if temperature > 37:
    print("The enzyme is active!")
else:
    print("The enzyme is inactive.")


# if-elif-else example
temperature = 25  # in degrees Celsius

if temperature > 37:
    print("The enzyme is highly active!")
elif 25 <= temperature <= 37:
    print("The enzyme is moderately active.")
else:
    print("The enzyme is inactive.")

# Using pass to skip a condition
hemoglobin = 13  # g/dL

if hemoglobin < 10:
    print("Anemia detected.")
elif 10 <= hemoglobin < 12:
    pass  # Ignore these cases for now
else:
    print("Normal hemoglobin levels.")


# Nested if example
gene_expression = 150  # arbitrary expression level units
temperature = 39  # in degrees Celsius

if gene_expression > 100:
    print("Gene is highly expressed.")
    if temperature > 37:
        print("And the enzyme is highly active at this temperature.")
    else:
        print("But the enzyme is inactive due to low temperature.")
else:
    print("Gene expression is low.")

#Summary:
#if condition: Evaluates a condition (e.g., enzyme activity or gene expression) and runs the block if true.
#elif: Optional, allows for checking additional conditions (e.g., moderate vs high activity).
#else: Optional, runs if none of the previous conditions are true.
#pass: Skips specific conditions when no action is required.
#Only one else: You can only have one else block at the end of the chain.
#Order matters: Conditions are checked in order, and only the first true condition’s block will execute.

-------------------------------------------Loops---------------------------------------------------------------------------
#Let's start by using a for loop to iterate over a list of gene names, which could represent genes we are studying.
# List of gene names
genes = ["BRCA1", "TP53", "EGFR", "MYC"]

# Using a for loop to iterate over the list
for gene in genes:
    print("Studying gene:", gene)

#Next, let's write a for loop to count the number of occurrences of each nucleotide (A, T, C, G) in a DNA sequence.
# DNA sequence
dna_sequence = "ATGCTAGCCTGA"

# Initialize a dictionary to store nucleotide counts
nucleotide_counts = {"A": 0, "T": 0, "C": 0, "G": 0}

# Iterate through the DNA sequence and count each nucleotide
for nucleotide in dna_sequence:
    if nucleotide in nucleotide_counts:
        nucleotide_counts[nucleotide] += 1

print("Nucleotide counts:", nucleotide_counts)

#compare two DNA sequences of equal length and print whether each position is a match or mismatch.
# Two DNA sequences of the same length
sequence_1 = "ATGCTA"
sequence_2 = "ATGCGT"

# Iterate through both sequences and compare nucleotides at each position
for i in range(len(sequence_1)):
    if sequence_1[i] == sequence_2[i]:
        print(f"Position {i+1}: Match ({sequence_1[i]})")
    else:
        print(f"Position {i+1}: Mismatch ({sequence_1[i]} vs {sequence_2[i]})"

#A while loop runs as long as a specified condition is true. Let' say we are monitoring gene expression levels, 
#and we want to keep running an experiment until the expression reaches a certain threshold.

# Starting gene expression level
gene_expression = 0

# Target expression level
target_expression = 100

# Simulating gene expression growth
while gene_expression < target_expression:
    gene_expression += 20  # Incrementing expression
    print("Current gene expression:", gene_expression)

print("Target expression reached!")

#The break keyword is used to exit the loop once the expression reaches or exceeds 90, even though the loop's condition allows it to continue up to 100.
# Starting gene expression level
gene_expression = 0

# Break the loop if expression reaches 90
while gene_expression < 100:
    gene_expression += 10
    print("Gene expression:", gene_expression)

    if gene_expression >= 90:
        print("Expression level high enough, stopping experiment.")
        break  # Exit the loop when expression reaches 90

#Summary
#for loops- Used to iterate over a sequence (like a list or string).
#Useful for tasks like counting nucleotides, processing gene lists, and comparing DNA sequences.
#while loops- Run as long as a condition is true.
#Ideal for tasks like monitoring gene expression or other experimental conditions.
#Loop control keywords - break: Exit a loop early if a specific condition is met.
#continue: Skip the rest of the current iteration and move to the next one.
#pass: Do nothing and move on.

#Exercises:

#Write a Python program to classify a patient's blood pressure into one of 
#three categories: "normal," "elevated," or "hypertension." Use the systolic and diastolic values as inputs.

#Normal: Systolic < 120 and Diastolic < 80
#Elevated: 120 ≤ Systolic < 130 and Diastolic < 80
#Hypertension: Systolic ≥ 130 or Diastolic ≥ 80

# Solution
systolic = 125
diastolic = 78

if systolic < 120 and diastolic < 80:
    print("Normal blood pressure.")
elif 120 <= systolic < 130 and diastolic < 80:
    print("Elevated blood pressure.")
else:
    print("Hypertension.")

#Write a program that iterates through a list of gene names and prints only those with a length greater than 5 characters.
# Solution
genes = ["BRCA1", "TP53", "MYC", "EGFR", "APC"]

for gene in genes:
    if len(gene) > 5:
        print("Gene with long name:", gene)

#Simulate a pH monitoring system where you start with a neutral pH of 7.0, 
#and with each iteration of a loop the pH decreases by 0.5. 
#Stop the loop when the pH drops below 5.0 and print "Acidic environment detected."

# Solution
ph = 7.0

while ph >= 5.0:
    print("Current pH:", ph)
    ph -= 0.5

print("Acidic environment detected.")

#Write a program to count the number of each nucleotide (A, T, C, G) in a given DNA sequence. 
#If an unknown character is found, print "Invalid nucleotide found" and skip that character.
# Solution
dna_sequence = "ATGCTAAGCTNTGC"

# Initialize nucleotide counts
nucleotide_counts = {"A": 0, "T": 0, "C": 0, "G": 0}

for nucleotide in dna_sequence:
    if nucleotide in nucleotide_counts:
        nucleotide_counts[nucleotide] += 1
    else:
        print("Invalid nucleotide found:", nucleotide)

print("Nucleotide counts:", nucleotide_counts)

#Write a program to compare two DNA sequences of the same length 
#and print out whether each position is a match or mismatch. 
#Keep track of how many matches and mismatches there are.

# Solution
sequence_1 = "ATGCTA"
sequence_2 = "ATGCGT"

matches = 0
mismatches = 0

for i in range(len(sequence_1)):
    if sequence_1[i] == sequence_2[i]:
        matches += 1
        print(f"Position {i+1}: Match ({sequence_1[i]})")
    else:
        mismatches += 1
        print(f"Position {i+1}: Mismatch ({sequence_1[i]} vs {sequence_2[i]})")

print(f"Total matches: {matches}")
print(f"Total mismatches: {mismatches}")

#Use the genetic code table to translate an RNA sequence into an amino acid sequence. 
#The RNA sequence is read in triplets (codons), and the translation stops if a "stop" codon (UAA, UAG, UGA) is encountered. 
#Use a dictionary for the genetic code and a loop to process the RNA sequence three characters at a time.

# Genetic code dictionary
genetic_code = {
    "AUG": "Methionine", "UUU": "Phenylalanine", "UUC": "Phenylalanine",
    "UUA": "Leucine", "UUG": "Leucine", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop"
    # Add more codons as needed
}

rna_sequence = "AUGUUUUUAAUGA"

# Initialize an empty list to store the amino acids
protein = []

# Iterate through the RNA sequence in steps of 3 (codon size)
for i in range(0, len(rna_sequence), 3):
    codon = rna_sequence[i:i+3]

    if len(codon) < 3:
        break  # Stop if the last codon is incomplete

    # Check if the codon is a stop codon
    if genetic_code[codon] == "Stop":
        print("Stop codon encountered. Translation terminated.")
        break

    # Append the corresponding amino acid to the protein sequence
    protein.append(genetic_code[codon])

# Output the translated amino acid sequence
print("Protein sequence:", protein)


-------------------------------------------AI_troubleshooting---------------------------------------------------------------

































