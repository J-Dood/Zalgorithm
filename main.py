# Z-Algorithm
# Jordan Dood 9/3/2021
# This program takes in a single FASTA file to search and prompts the user for a query sequence.
# This program returns a list of full matches to the command line, along with their start index.


def z_algorithm(sequence, query):
    """
    :argument sequence is the string to be searched by
    :argument query is the string searched for in the sequence
    :rtype: list
    :returns a list of the start indexes of full matches of query in sequence
    """
    # Construct the single string
    combined_string = query + "$" + sequence
    # Create empty list for matches
    match_list = []
    # Create empty list for all Z values
    z_list = [-1]
    # Initiate the variables for the z-algorithm
    right_edge = 0
    left_edge = 0
    k_prime = 0
    beta = 0

    # Loop over each letter in the combined string
    for k in range(1, len(combined_string)+1):
        temp_z = 0
        increment = 0
        # Case 1: Not in Z-box
        if k > right_edge:
            # Match letters until mismatch found
            while ((len(combined_string) - 1 >= k + increment)
                    and (combined_string[increment] == combined_string[k + increment])):
                temp_z += 1
                increment += 1
            if temp_z > 0:
                right_edge = k + temp_z - 1
                left_edge = k
            z_list.append(temp_z)
        # Case 2: In Z-box
        else:
            k_prime = k - left_edge + 1
            beta = right_edge - k + 1
            # Case 2: Sub-case 1: Matching does not extend past Z-box end
            if z_list[k_prime] < beta:
                temp_z = z_list[k_prime]
            # Case 2: Sub-case 2: Matching extends past Z-box end
            else:
                # Match letters past Z-box end until mismatch found
                while ((len(combined_string) - 1 >= right_edge + 1 + increment)
                       and (combined_string[beta + 1 + increment] == combined_string[right_edge + 1 + increment])):
                    temp_z += 1
                    increment += 1
                temp_z = beta + increment - 1
                # Update Z-box edges
                right_edge = increment - 1
                left_edge = k
            z_list.append(temp_z)

    # Use Z scores to compose list of indexes where perfect matches were found
    for x in range(len(z_list) -1):
        if z_list[x] == len(query):
            match_list.append(x - len(query))
    return match_list


def main():
    """handels the input, calls the z-algorithm, and prints the output"""
    # Handel the file opening of the FASTA file to be searched
    file = open(r"C:\Users\doodw\PycharmProjects\Zalgorithm\venv\input\testfasta.fasta", 'r')
    # Create an empty string to hold text sequence
    sequence = ""
    # Flag to ensure only one sequence is read from the Fasta
    single_sequence_flag = False
    # Read in and concatenate sequence lines from Fasta file
    for line in file:
        if line.find('>') != -1:    # Header ignored
            # Ignores multiple sequences in Fasta file
            if single_sequence_flag:
                break
            continue
        else:
            line = line.strip().lower()
            sequence += line

    # Prompt the user for query
    query = input("Please enter a DNA sequence to search for: ").strip().lower()

    # Call the z-algorithm
    output = z_algorithm(sequence, query)

    # Print the output form the z-algorithm
    for index in output:
        print("A perfect match found at: " + str(index))


# Call main
main()