def readData (fileName):
    """
    Opening file, sifting through data, returning similarities.
    """

    file = open(fileName, 'r')
    file = file.readlines()

    time1 = round(float(file[0]), 2) # time for max score, indices

    score = file[1] # maxscore

    word1 = file[2] # first sequences
    word2 = file[3] # second sequence

    time2 = round(float(file[5]), 2) # time to run Hirschberg on sequence

    gap = 0
    match = 0
    mismatch = 0

    for i in range(len(word1)): # comparing species to human protein
        if (word1[i] == '-' or word2[i] == '-'):
            gap += 1
        elif (word1[i] == word2[i]):
            match += 1
        elif (word1[i] != word2[i]):
            mismatch += 1

    return time1, time2, score, match, mismatch, gap, (match/len(word1)), len(word1)

readData ('catout.txt')
